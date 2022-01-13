%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                   By Fernando Leon-Cazares                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Monolayer interfacial trap
% Solved via an explicit finite differences method
% Zero flux (left) and constant concentration (right)
%
% Inputs:
% xs:       distance modelled [m]
% xi:       distance where phase1 ends [m]
% T0:       temperature [K]
% dx:       length step size [m]
% dt:       time step size [s]
% tend:     final time simulated [s]
% S:        lattice site number density [mol/m3]
% D0:       Diffusivity prefactor [m2/s]
% Q:        Energy barrier in the bulk [J/mol]
% Eb:       Interfacial energy barrier [J/mol]
% DeltaE:   Energy difference between phases [J/mol]
% DeltaEt:  Energy difference between crystal 1 and the interfacial trap [J/mol]
% a:        Layer spacing [m]
% c#_ini:   H initial concentration [mol/m3]
% coord:    Coordinate system (1->Cartesian, 2->Cylindrical, 3->Spherical)
% l:        Sample dimensions (Cartesian [L,w], Cylindrical [L,~] or Spherical [~])
%
% Outputs:
% xr:       Spatial grid [m]
% x:        Spatial grid + two points at the side of the interface [m]
% t:        Time points simulated [s]
% ts:       Time points saved [s]
% clr:      Concentration in the lattice sites [mol/m3]
% cl:       Concentration at all x points [mol/m3]
% civ:      Concentration at the [left of interface; interface site; right of interface] [mol/m3]
% Js1:      Flux through left interface [mol/s]
% Js2:      Flux through right interface [mol/s]
% Ji:       Flux at the interface [mol/s]
% H:        Total amount of solute in the sample [mol]

 %% Function
function [xr,x,t,ts,clr,cl,civ,Js1,Js2,Ji,H] = traps1D_open_t(xs,xi,T0,dx,dt,tend,S1,D01,Q1,Eb1,S2,D02,Q2,Eb2,St,DeltaE,DeltaEt,D0t,a,cs2,c1_ini,c2_ini,coord,l,nsvs)

%% Geometry
x = 0:dx:xs;
if x(end) ~= xs
    error('xs must be a multiple of dx')
end

t = 0:dt:tend;
T = zeros(size(t)) + T0;
disp(['Monolayer interfacial trap: tend = ',num2str(tend),'s',', dt = ',num2str(dt)])

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
if k ~= 0.5
    error('The interface must be situated at the centre of a cell')
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
cl(end,:) = cs2;
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
Db1 = 2*D01*G(Q1+Eb1,T(1))*D0t*G(Ebt1,T(1))*S1*St/(D01*G(Q1+Eb1,T(1))*S1+D0t*G(Ebt1-DeltaE1t,T(1))*St);
Db2 = 2*D0t*G(Ebt2,T(1))*D02*G(Q2+Eb2,T(1))*St*S2/(D0t*G(Ebt2,T(1))*St+D02*G(Q2+Eb2+DeltaE2t,T(1))*S2);
[clj(I),HI,civ(1,1),civ(2,1),civ(3,1)] = redist(c1_ini*Vi(1)+c1_ini*Vi(2)+c2_ini*Vi(3),DeltaE1t,DeltaE2t,clj(I-1),clj(I+1),dx,x(I-1),l,Vi(2),D1,D2,Db1,Db2,phi1,phi2,phib1,phib2,S1,S2,St,R,T(1),coord);
HIv(1) = HI;
cljr = clj(1:end-1);                    % clj reduced to ignore both boundaries
Jx_l = zeros(size(clj));                % Flux at each r_l element boundary (Jr_l(1)=0 always due to symmetry)
jss = max(js(1),2);
ji = 1;
%% Calculation - for loop
for j = 1:(length(t)-1)                 % T loop
    Tj = T(j);
    D1 = D01*G(Q1,Tj);
    D2 = D02*G(Q2,Tj);
    Db1 = 2*D01*G(Q1+Eb1,Tj)*D0t*G(Ebt1,Tj)*S1*St/(D01*G(Q1+Eb1,Tj)*S1+D0t*G(Ebt1-DeltaE1t,Tj)*St);
    Db2 = 2*D0t*G(Ebt2,Tj)*D02*G(Q2+Eb2,Tj)*St*S2/(D0t*G(Ebt2,Tj)*St+D02*G(Q2+Eb2+DeltaE2t,Tj)*S2);
    D_T = zeros(size(clj)) + D1;      % Diffusion coefficient for each temperature [m2/s] - right boundary of each element
    D_T(I:end) = D2;                    
    
    %%% Fluxes
    Jx_l(2:end) = -D_T(1:end-1) .* (clj(2:end) - clj(1:end-1))/dx;          % Fick's 1st law - the result does not apply at the boundaries of element I
    
                            % Interface fluxes
    Jx_l(I)   = -D1*Db1*S1/(Db1*phi1 + D1*S1*phib1*G(DeltaE1t,Tj)) * (clj(I)/St*G(DeltaE1t,Tj) - clj(I-1)/S1)/dx;
    Jx_l(I+1) = -D2*Db2*S2/(Db2*phi2 + D2*S2*phib2)                * (clj(I+1)/S2 - clj(I)/St*G(DeltaE2t,Tj))/dx;

    Jx0 = Jx_l(1:end-1).*Ax_l(1:end-1) - Jx_l(2:end).*Ax_r(1:end-1);        % Flux into each element * the boundary area [mol/s]
    
    % Fick's 2nd law    
    cljr = Jx0./Vx(1:end-1) * dt + cljr;    % dcl/dt
    clj(1:end-1) = cljr;                    % Updating cl
    
    %%% Redistribution of solute in the heterogeneous element
    [clj(I),HI,ci(1),ci(2),ci(3)] = redist(HI+Jx0(I)*dt,DeltaE1t,DeltaE2t,clj(I-1),clj(I+1),dx,x(I-1),l,Vi(2),D1,D2,Db1,Db2,phi1,phi2,phib1,phib2,S1,S2,St,R,Tj,coord);

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

Js1 = 0;   % Flux at the surfaces [mol/s]
Js2 = -D_Ts(end,:) .* (cl(end,:)-cl(end-1,:))/dx * As2;

% To incorporate two concentration values at the interface (#r = variables without the interface values)
xr = x;
clr = cl;
x = [x(1:I-1),xi-a,xi,xi+a,x(I+1:end)];
cl = [cl(1:I-1,:);civ;cl(I+1:end,:)];

end

%% Solute redistribution functions
function [c2,HI,ci1,ct,ci2] = redist(H,DeltaE1t,DeltaE2t,c1,c2,dx,x1,l,Vt,D1,D2,Db1,Db2,phi1,phi2,phib1,phib2,S1,S2,St,R,T,coord)   % Redistribution of solute in an interface element - equal fluxes D1-Db1-Db2-D2 - integral in terms of phi
    
    if coord == 1
        A = l(1)*l(2);
        m1 = (A*Db1*dx*exp(-DeltaE1t/(R*T))*S1*(4*phi1^2-1)) / (8*St*(Db1*phi1+D1*S1*phib1));
        b1 = (A*c1*dx*(2*phi1-1)*(Db1*(2*phi1-1)+4*D1*S1*phib1)) / (8*(Db1*phi1+D1*S1*phib1));
        m2 = (A*Db2*dx*exp(-DeltaE2t/(R*T))*S2*(4*phi2^2-1)) / (8*St*(Db2*phi2+D2*S2*phib2));
        b2 = (A*c2*dx*(2*phi2-1)*(Db2*(2*phi2-1)+4*D2*S2*phib2)) / (8*(Db2*phi2+D2*S2*phib2));
    elseif coord == 2
        L = l(1);
        m1 = (Db1*dx*exp(-DeltaE1t/(R*T))*L*pi*S1*(3*x1*(4*phi1^2-1)+dx*(8*phi1^3-1))) / (12*St*(Db1*phi1+D1*S1*phib1));
        b1 = (c1*dx*L*pi*(2*phi1-1)*(Db1*(2*phi1-1)*(dx+3*x1*dx*phi1)+3*D1*S1*(dx+4*x1+2*dx*phi1)*phib1)) / (12*(Db1*phi1+D1*S1*phib1));
        m2 = (Db2*dx*exp(-DeltaE2t/(R*T))*L*pi*S2*(3*x1*(4*phi1^2-1)+dx*(-5-8*(phi2-3)*phi2^2))) / (12*St*(Db2*phi2+D2*S2*phib2));
        b2 = (c2*dx*L*pi*(2*phi2-1)*(Db2*(3*x1-dx*(phi2-5))*(2*phi2-1)+3*D2*S2*(4*x1+dx*(7-2*phi2))*phib2)) / (12*(Db2*phi2+D2*S2*phib2));
    elseif coord == 3
        m1 = (Db1*dx*exp(-DeltaE1t/(R*T))*pi*S1*(24*x1^2*(4*phi1^2-1)+16*dx*x1*(8*phi1^3-1)+dx^2*(48*phi1^4-3))) / (48*St*(Db1*phi1+D1*S1*phib1));
        b1 = (c1*dx*pi*(2*phi1-1)*(Db1*(2*phi1-1)*(24*x1^2+16*dx*x1*(phi1+1)+dx^2*(3+4*phi1*(1+phi1)))+8*D1*S1*(12*x1^2+6*dx*(x1+2*x1*phi1)+dx^2*(1+2*phi1+4*phi1^2))*phib1)) / (48*(Db1*phi1+D1*S1*phib1));
        m2 = (exp(-DeltaE2t/(R*T))*pi*(1/16*Db2*S2*(3*dx+2*x1)^3*(7*dx+2*x1)-Db2*S2*(x1-dx*(phi2-2))^3*(x1+dx*(3*phi2+2)))) / (3*dx*St*(Db2*phi2+D2*S2*phib2));
        b2 = (c2*dx*pi*(2*phi2-1)*(Db2*(2*phi2-1)*(24*x1^2-16*dx*x1*(phi2-5)+dx^2*(67+4*(-7+phi2)*phi2))+8*D2*S2*(12*x1^2+6*dx*x1*(7-2*phi2)+dx^2*(37-22*phi2+4*phi2^2))*phib2)) / (48*(Db2*phi2+D2*S2*phib2));
    end
    c2 = (H-b1-b2)/(Vt+m1+m2);
    ct = c2;
    ci1 = S1*(ct*Db1*phi1*exp(-DeltaE1t/(R*T))+c1*D1*St*phib1)/(St*(Db1*phi1+D1*S1*phib1));
    ci2 = S2*(ct*Db2*phi2*exp(-DeltaE2t/(R*T))+c2*D2*St*phib2)/(St*(Db2*phi2+D2*S2*phib2));
    HI = H;
end