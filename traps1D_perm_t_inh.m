%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                     1D spherical trap                         %%%%%
%%%%%                   By Fernando Leon-Cazares                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extension (simplified) from multitrap_TDS_spheretraps_4
%
% Spherical 1D H diffusion equation in 2 phases
% 
% Solved via an explicit finite differences method
%
%
% Inputs:
% xs:       distance modelled [m]
% xi:       distance where phase1 ends [m]
%        --> ** given as a close value to the desired size but updated to match the mesh
% phi:      heating rate [K/s]
% cs:       H surface concentration [mol/m3]
% c#_ini:   H initial concentration [mol/m3]
% coord:    Coordinate system (1->Cartesian, 2->Cylindrical, 3->Spherical)
% l:        Sample dimensions (Cartesian [L,w], Cylindrical [L,~] or Spherical [~])
%
% Isothermal case: phi -> tend, dT -> dt
%
% Nomenclature:
% H:        Total H in the sample [mol]
% cs:       Concentration of H at the free surface
% cl:       Concentration of H in free sites
% clr:      cl without the interface values
% civ:      Concentrations at the [left side of the interface, trap, right side of the interface]
% Ht:       H in phase1 [mol]
% Js:       H flux at the surface of phase 2 [mol/s]
% Ji:       H flux at the interface [mol/s]

% Updates:
% spheretraps_1D_TDS
% _2: Simplified flux calculation
% _3: Flux at the interface given by modified Ilin's model 
% traps_1D_TDS_1: coord -> For both spherical and cartesian coordinates
%                 r changed for x
% _2: Flux at the interface given by gradient of mu
%     Introducing fds
% _perm: fixed concentration at both ends (it doesn't make sense for cylindrical and spherical coordinates but whatever)
% _3: Added intdiff = 5 for the local heterogeneous model to generalise energy barriers
% _Iel: Now the interface is in the middle of one element
% _4: Interface at any point
% _5: New homogenisation step: Equal fluxes Dleq-Dreq. Choosing method via variable 'redist'
% _6: New homogenisation step: Equal fluxes D1-Db-D2
% _t: With a trap of thickness a at the interface
% 7: Changed trap homogeneisation scheme to phi values
% _inheq: Inhomogeneous interface assuming a constant thermodynamic equilibrium (different lattice sites)

%% Function
function [xr,x,t,ts,clr,cl,ct,civ,Js1,Js2,Ji,H,k] = traps1D_perm_t_inh(xs,xi,T0,dx,dt,tend,S1,D01,Q1,Eb1,S2,D02,Q2,Eb2,St,DeltaE,DeltaEt,D0t,a,cs1,cs2,c1_ini,c2_ini,coord,l,nsvs)

if numel(unique([length(St),length(DeltaEt),length(D0t),length(Eb1),length(Eb2)])) ~= 1
   error('Trap parameters must have the same length')
end

%% Geometry
x = 0:dx:xs;
if x(end) ~= xs
    error('xs must be a multiple of dx')
end

t = 0:dt:tend;
T = zeros(size(t)) + T0;
disp(['Inhomogeneous: tend = ',num2str(tend),'s',', dt = ',num2str(dt)])

if nsvs > length(t)
    js = 2:length(t);
else
    js = round(length(t).*(1:nsvs)/nsvs);   % Indices where data is saved
end
js = unique(js(js>1));                  % In case nsvs>=length(t)
ts = [0,t(js)];    % Recorded data
js = [js,0];

% To calculate the dimensions of each element
x_l = ([x(1),x(1:end-1)] + x)/2;            % Boundaries of each element given by [x_l,x_r]
x_r = (x + [x(2:end),x(end)])/2;
I = find((xi>=x_l)&(xi<x_r),1);             % Index of the element with the trap
k = (xi-x_l(I))/(x_r(I)-x_l(I));            % Fraction of the element that corresponds to material 1
It = 1;                                     % Index of the trap site to reference all equations to

if k ~= 0.5
    error('k must be 0.5, i.e. the interface must fall at the centre of a cell')
end
phi1 = 1 - a/dx;
phib1 = a/dx;
phib2 = a/dx;
phi2 = 1- a/dx;

if coord == 1
    w = l(2);
    L = l(1);
    Ax_l = L*w*ones(size(x_l))';            % Left boundary area for each element [m^2]
    Ax_r = Ax_l;                            % Right boundary area for each element [m^2]
    Vx = L*w*(x_r-x_l)';                    % Volume of each element [m^3]
    As1 = L*w;                              % Area of the surface [m^2]
    As2 = L*w;
    Ai = [L*w,L*w];                         % Areas of the interfaces [m^2]
    Vi = [L*w*((xi-a/2)-x_l(I)),L*w*a,L*w*(x_r(I)-(xi+a/2))];   % Volume of the left/trap/right regions [m^3]
elseif coord == 2
    L = l(1);
    Ax_l = (2*pi*x_l*L)';
    Ax_r = (2*pi*x_r*L)';
    Vx = pi*(x_r.^2-x_l.^2)'*L;
    As1 = 0;
    As2 = 2*pi*xs*L;
    Ai = [2*pi*(xi-a/2)*L,2*pi*(xi+a/2)*L];
    Vi = [pi*((xi-a/2)^2-x_l(I)^2)*L,pi*((xi+a/2)^2-(xi-a/2)^2)*L,pi*(x_r(I)^2-(xi+a/2)^2)*L];   
elseif coord == 3
    Ax_l = (4*pi*x_l.^2)';
    Ax_r = (4*pi*x_r.^2)';
    Vx = 4/3*pi*(x_r.^3-x_l.^3)';
    As1 = 0;
    As2 = 4*pi*xs^2;
    Ai = [4*pi*(xi-a/2)^2,4*pi*(xi+a/2)^2];
    Vi = [4/3*pi*((xi-a/2)^3-x_l(I)^3),4/3*pi*((xi+a/2)^3-(xi-a/2)^3),4/3*pi*(x_r(I)^3-(xi-a/2)^3)];
end

%% State variable initialisation and paramters to use
cl = zeros(length(x),length(ts));           % Lattice H concentration
cl(:,1) = c1_ini;                           % Boundary conditions
cl(I:end,1) = c2_ini;
cl(1,:) = cs1;
cl(end,:) = cs2;
ct = zeros(length(DeltaEt),length(ts));
HIv = zeros(1,length(ts));                  % Hydrogen in the regions of element I
civ = zeros(3,length(ts));                  % Concentrations at the interface (ci1,ct,ci2)
Ji = zeros(2,length(ts));
D_Ts = zeros(size(cl));

% Energetics and trap layout
R = 8.31446261815324;                       % Gas constant [J/molK]
G = @(Q,T) exp(-Q/(R*T));                   % Fraction of succesful jumps [-]
Ebt1 = DeltaEt+Q1+Eb1;
Ebt2 = DeltaEt+DeltaE+Q2+Eb2;
DeltaE1t = DeltaEt;
DeltaE2t = DeltaEt+DeltaE;

%% Stability: mk1 < 1/2
mk1 = max([D01.*exp(-Q1./(R*max(T))),D02.*exp(-Q2./(R*max(T)))])*dt/dx^2;
disp(['DeltaE = ',num2str(DeltaE),' J'])
disp(['Maximum mesh size = ',num2str(mk1)])
if mk1 >= 1/2
   error('Unstable solution: mk1<1/2 required') 
end

%% Calculation - step 1
wb = waitbar(0,'Progress:');
clj = cl(:,1);                          % cl at that point in time ...
D1 = D01*G(Q1,T(1));
D2 = D02*G(Q2,T(1));
Db1 = 0;
Db2 = 0;
for k = 1:length(DeltaEt)
    Db1 = Db1+2*D01*G(Q1+Eb1(k),T(1))*D0t(k)*G(Ebt1(k),T(1))*S1*St(k)/(D01*G(Q1+Eb1(k),T(1))*S1+D0t(k)*G(Ebt1(k)-DeltaE1t(k),T(1))*St(k));
    Db2 = Db2+2*D0t(k)*G(Ebt2(k),T(1))*D02*G(Q2+Eb2(k),T(1))*St(k)*S2/(D0t(k)*G(Ebt2(k),T(1))*St(k)+D02*G(Q2+Eb2(k)+DeltaE2t(k),T(1))*S2);
end
[ctj(1),HI,civ(1,1),~,civ(3,1)] = redist(c1_ini*Vi(1)+c1_ini*Vi(2)+c2_ini*Vi(3),DeltaE1t,DeltaE2t,clj(I-1),clj(I+1),dx,x(I-1),l,D1,D2,Db1,Db2,phi1,phi2,phib1,phib2,S1,S2,St(It),St,Vi(2),R,T(1),coord);
ctj = ctj(1)*St/St(1).*exp(-(DeltaEt(1)-DeltaEt)/(R*T(1)));
clj(I) = sum(ctj);
civ(2,1) = clj(I);
ct(:,1) = ctj;   
HIv(1) = HI;
cljr = clj(2:end-1);                    % clj reduced to ignore both boundaries
Jx_l = zeros(size(clj));                % Flux at each r_l element boundary (Jr_l(1)=0 always due to symmetry)
jss = max(js(1),2);
ji = 1;
%% Calculation - for loop
for j = 1:(length(t)-1)                 % T loop
    Tj = T(j);
    D1 = D01*G(Q1,Tj);
    D2 = D02*G(Q2,Tj);
    Db1 = 0;
    Db2 = 0;
    for k = 1:length(DeltaEt)
        Db1 = Db1+2*D01*G(Q1+Eb1(k),Tj)*D0t(k)*G(Ebt1(k),Tj)*S1*St(k)/(D01*G(Q1+Eb1(k),Tj)*S1+D0t(k)*G(Ebt1(k)-DeltaE1t(k),Tj)*St(k));
        Db2 = Db2+2*D0t(k)*G(Ebt2(k),Tj)*D02*G(Q2+Eb2(k),Tj)*St(k)*S2/(D0t(k)*G(Ebt2(k),Tj)*St(k)+D02*G(Q2+Eb2(k)+DeltaE2t(k),Tj)*S2);
    end
    D_T = zeros(size(clj)) + D1;      % Diffusion coefficient for this temperature [m2/s] - right boundary of each element
    D_T(I:end) = D2;                    
    
    %%% Fluxes
    Jx_l(2:end) = -D_T(1:end-1) .* (clj(2:end) - clj(1:end-1))/dx;          % Fick's 1st law - the result does not apply at the boundaries of element I
    
                            % Interface fluxes
    Jx_l(I)   = -D1*Db1*S1/(Db1*phi1 + D1*S1*phib1*G(DeltaE1t(1),Tj)) * (ctj(1)/St(1)*G(DeltaE1t(1),Tj) - clj(I-1)/S1)/dx;
    Jx_l(I+1) = -D2*Db2*S2/(Db2*phi2 + D2*S2*phib2)                   * (clj(I+1)/S2 - ctj(1)/St(1)*G(DeltaE2t(1),Tj))/dx;
    
    Jx0 = Jx_l(2:end-1).*Ax_l(2:end-1) - Jx_l(3:end).*Ax_r(2:end-1);        % Flux into each element * the boundary area [mol/s]

    % Fick's 2nd law    
    cljr = Jx0./Vx(2:end-1) * dt + cljr;    % dcl/dt
    clj(2:end-1) = cljr;                    % Updating cl
    
    %%% Redistribution of solute in the heterogeneous element
    [ctj(1),HI,ci(1),~,ci(3)] = redist(HI+Jx0(I-1)*dt,DeltaE1t,DeltaE2t,clj(I-1),clj(I+1),dx,x(I-1),l,D1,D2,Db1,Db2,phi1,phi2,phib1,phib2,S1,S2,St(1),St,Vi(2),R,Tj,coord);
    ctj = ctj(1)*St/St(1).*exp(-(DeltaEt(1)-DeltaEt)/(R*Tj));
    clj(I) = sum(ctj);
    ci(2) = clj(I);
    
    % Updating cl and cr
    if j == 1
        D_Ts(:,1) = D_T;
        Ji(:,1) = Jx_l([I,I+1]).*Ai';
    end
    if j == jss-1
        waitbar(j/(length(t)-1),wb)
        ji = ji+1;
        jss = js(ji);
        cl(:,ji) = clj;
        ct(:,ji) = ctj;
        HIv(ji) = HI;
        civ(:,ji) = ci;
        D_Ts(:,ji) = D_T;
        Ji(:,ji) = Jx_l([I,I+1]).*Ai';
    end
end
close(wb);

%% Final calculations 
% To calculate H and Js
H = sum(cl([1:I-1,I+1:end],:).*Vx([1:I-1,I+1:end])) + HIv;

Js1 = -D_Ts(1,:) .* (cl(2,:)-cl(1,:))/dx * As1;   % Flux at the surfaces [mol/s]
Js2 = -D_Ts(end,:) .* (cl(end,:)-cl(end-1,:))/dx * As2;

% To incorporate two concentration values at the interface (#r = variables without the interface values)
xr = x;
clr = cl;
x = [x(1:I-1),xi-a,xi,xi+a,x(I+1:end)];
cl = [cl(1:I-1,:);civ;cl(I+1:end,:)];

end

%% Solute redistribution functions
function [c2,HI,ci1,ct,ci2] = redist(H,DeltaE1t,DeltaE2t,c1,c2,dx,x1,l,D1,D2,Db1,Db2,phi1,phi2,phib1,phib2,S1,S2,SIt,St,Vt,R,T,coord)   % Redistribution of solute in an interface element - equal fluxes D1-Db1-Db2-D2 - integral in terms of phi
    
    if coord == 1
        A = l(1)*l(2);
        m1 = (A*Db1*dx*exp(-DeltaE1t(1)/(R*T))*S1*(4*phi1^2-1)) / (8*SIt *(Db1*phi1+D1*S1*phib1));
        b1 = (A*c1*dx*(2*phi1-1)*(Db1*(2*phi1-1)+4*D1*S1*phib1)) / (8*(Db1*phi1+D1*S1*phib1));
        m2 = (A*Db2*dx*exp(-DeltaE2t(1)/(R*T))*S2*(4*phi2^2-1)) / (8*SIt *(Db2*phi2+D2*S2*phib2));
        b2 = (A*c2*dx*(2*phi2-1)*(Db2*(2*phi2-1)+4*D2*S2*phib2)) / (8*(Db2*phi2+D2*S2*phib2));
    elseif coord == 2
        L = l(1);
        m1 = (Db1*dx*exp(-DeltaE1t(1)/(R*T))*L*pi*S1*(3*x1*(4*phi1^2-1)+dx*(8*phi1^3-1))) / (12*SIt *(Db1*phi1+D1*S1*phib1));
        b1 = (c1*dx*L*pi*(2*phi1-1)*(Db1*(2*phi1-1)*(dx+3*x1*dx*phi1)+3*D1*S1*(dx+4*x1+2*dx*phi1)*phib1)) / (12*(Db1*phi1+D1*S1*phib1));
        m2 = (Db2*dx*exp(-DeltaE2t(1)/(R*T))*L*pi*S2*(3*x1*(4*phi1^2-1)+dx*(-5-8*(phi2-3)*phi2^2))) / (12*SIt *(Db2*phi2+D2*S2*phib2));
        b2 = (c2*dx*L*pi*(2*phi2-1)*(Db2*(3*x1-dx*(phi2-5))*(2*phi2-1)+3*D2*S2*(4*x1+dx*(7-2*phi2))*phib2)) / (12*(Db2*phi2+D2*S2*phib2));
    elseif coord == 3
        m1 = (Db1*dx*exp(-DeltaE1t(1)/(R*T))*pi*S1*(24*x1^2*(4*phi1^2-1)+16*dx*x1*(8*phi1^3-1)+dx^2*(48*phi1^4-3))) / (48*SIt *(Db1*phi1+D1*S1*phib1));
        b1 = (c1*dx*pi*(2*phi1-1)*(Db1*(2*phi1-1)*(24*x1^2+16*dx*x1*(phi1+1)+dx^2*(3+4*phi1*(1+phi1)))+8*D1*S1*(12*x1^2+6*dx*(x1+2*x1*phi1)+dx^2*(1+2*phi1+4*phi1^2))*phib1)) / (48*(Db1*phi1+D1*S1*phib1));
        m2 = (exp(-DeltaE2t(1)/(R*T))*pi*(1/16*Db2*S2*(3*dx+2*x1)^3*(7*dx+2*x1)-Db2*S2*(x1-dx*(phi2-2))^3*(x1+dx*(3*phi2+2)))) / (3*dx*SIt *(Db2*phi2+D2*S2*phib2));
        b2 = (c2*dx*pi*(2*phi2-1)*(Db2*(2*phi2-1)*(24*x1^2-16*dx*x1*(phi2-5)+dx^2*(67+4*(-7+phi2)*phi2))+8*D2*S2*(12*x1^2+6*dx*x1*(7-2*phi2)+dx^2*(37-22*phi2+4*phi2^2))*phib2)) / (48*(Db2*phi2+D2*S2*phib2));
    end
    mt = sum(Vt.*St/SIt .*exp(-(DeltaE1t(1)-DeltaE1t)/(R*T)));
    c2 = (H-b1-b2)/(mt+m1+m2);
    ct = c2;
    ci1 = S1*(ct*Db1*phi1*exp(-DeltaE1t(1)/(R*T))+c1*D1*SIt*phib1)/(SIt *(Db1*phi1+D1*S1*phib1));
    ci2 = S2*(ct*Db2*phi2*exp(-DeltaE2t(1)/(R*T))+c2*D2*SIt*phib2)/(SIt *(Db2*phi2+D2*S2*phib2));
    HI = H;
end
