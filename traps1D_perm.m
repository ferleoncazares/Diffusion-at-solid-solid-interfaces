%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                   By Fernando Leon-Cazares                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interfacial energy barrier
% Solved via an explicit finite differences method
% Constant concentration at the surfaces
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
% a:        Layer spacing [m]
% cs:       Surface concentration [mol/m3]
% c#_ini:   H initial concentration [mol/m3]
% coord:    Coordinate system (1->Cartesian, 2->Cylindrical, 3->Spherical)
% L:        Sample dimensions (Cartesian [L,w], Cylindrical [L,~] or Spherical [~])
%
% Outputs:
% x:        Spatial grid [m]
% t:        Time points simulated [s]
% ts:       Time points saved [s]
% cl:       Concentration in the lattice sites [mol/m3]
% Js1:      Flux through left interface [mol/s]
% Js2:      Flux through right interface [mol/s]
% Ji:       Flux at the interface (left and right) [mol/s]
% H:        Total amount of solute in the sample [mol]
% DeltaE:   Energy difference between phases [J/mol]

function [x,t,ts,cl,Js1,Js2,Ji,H,DeltaE] = traps1D_perm(xs,xi,T0,dx,dt,tend,S1,D01,Q1,Eb1,S2,D02,Q2,Eb2,a,cs1,cs2,c1_ini,c2_ini,coord,L,nsvs)

%% Geometry
x = 0:dx:xs;
if x(end) ~= xs
    error('xs must be a multiple of dx')
end

t = 0:dt:tend;
T = zeros(size(t)) + T0;
disp(['Energy barrier: tend = ',num2str(tend),'s',', dt = ',num2str(dt)])

if nsvs > length(t)
    js = 2:length(t);
else
    js = round(length(t).*(1:nsvs)/nsvs);   % Indices where data is saved
end
ts = [0,t(js)];    % Recorded data
js = [js,0];

% To calculate the dimensions of each element
x_l = ([x(1),x(1:end-1)] + x)/2;            % Boundaries of each element given by [x_l,x_r]
x_r = (x + [x(2:end),x(end)])/2;
I = find((xi>=x_l)&(xi<x_r),1);             % Index of the first element with phase2
k = (xi-x_l(I))/(x_r(I)-x_l(I));            % Fraction of the element that corresponds to material 1
if k ~= 0
    error('The interface must be situated in the middle of two cells')
end
h = 1/2 + k;                                % Fraction of the Ji region to the left of the interface
phi1 = h-a/(2*dx);
phib = a/dx;
phi2 = 1-h-a/(2*dx);
if coord == 1
    w = L(2);
    L = L(1);
    Ax_l = L*w*ones(size(x_l))';            % Left boundary area for each element [m^2]
    Ax_r = Ax_l;                            % Right boundary area for each element [m^2]
    Vx = L*w*(x_r-x_l)';                    % Volume of each element [m^3]
    As1 = L*w;                              % Area of the surface [m^2]
    As2 = L*w;
    Ai = L*w;                               % Area of the interface [m^2]
elseif coord == 2
    L = L(1);
    Ax_l = (2*pi*x_l*L)';
    Ax_r = (2*pi*x_r*L)';
    Vx = pi*(x_r.^2-x_l.^2)'*L;
    As1 = 0;
    As2 = 2*pi*xs*L;
    Ai = 2*pi*xi*L;
elseif coord == 3
    Ax_l = (4*pi*x_l.^2)';
    Ax_r = (4*pi*x_r.^2)';
    Vx = 4/3*pi*(x_r.^3-x_l.^3)';
    As1 = 0;
    As2 = 4*pi*xs^2;
    Ai = 4*pi*xi^2;
end

%% State variable initialisation
cl = zeros(length(x),length(ts));           % Lattice H concentration
cl(:,1) = c1_ini;                           % Boundary conditions
cl(I:end,1) = c2_ini;
cl(1,:) = cs1;
cl(end,:) = cs2;
Ji = zeros(size(ts));
D_Ts = zeros(size(cl));

% Energetics and trap layout
R = 8.31446261815324;                       % Gas constant [J/molK]
G = @(Q,T) exp(-Q/(R*T));                   % Fraction of succesful jumps [-]
DeltaE = Q1 + Eb1 - Q2 - Eb2;               % Difference in energy state E2-E1 [J]

%% Stability: mk1 < 1/2
mk1 = max([D01.*exp(-Q1./(R*max(T))),D02.*exp(-Q2./(R*max(T)))])*dt/dx^2;
disp(['DeltaE = ',num2str(DeltaE),' J'])
disp(['Maximum mesh size = ',num2str(mk1)])
if mk1 >= 1/2
   error('Unstable solution: mk1<1/2 required') 
end

%% Calculation - Step 1
wb = waitbar(0,'Progress:');
clj = cl(:,1);                          % cl at that point in time ...

cljr = clj(2:end-1);                    % clj reduced to ignore both boundaries
Jx_l = zeros(size(clj));                % Flux at each r_l element boundary (Jr_l(1)=0 always due to symmetry)
jss = js(1);
ji = 1;
%% Calculation - for loop
for j = 1:(length(t)-1)                 % T loop
    Tj = T(j);
    D1 = D01*G(Q1,Tj);
    D2 = D02*G(Q2,Tj);
    Db = 2*D01*G(Q1+Eb1,Tj)*D02*G(Q2+Eb2,Tj)*S1*S2/(D01*G(Q1+Eb1,Tj)*S1+D02*G(Q2+Eb2+DeltaE,Tj)*S2);
    D_T = zeros(size(clj)) + D1;      % Diffusion coefficient for each temperature [m2/s] - right boundary of each element
    D_T(I:end) = D2;                    
    
    %%% Fluxes
    Jx_l(2:end) = -D_T(1:end-1) .* (clj(2:end) - clj(1:end-1))/dx;          % Fick's 1st law - the result does not apply at the boundaries of element I
    
                            % Interface flux
    Jx_l(I)   = -D1*D2*Db*S1*S2/(D1*Db*S1*phi2 + D2*S2*(D1*S1*phib + Db*phi1*G(DeltaE,Tj))) * (clj(I)/S2 - clj(I-1)/S1*G(DeltaE,Tj))/dx; % Flux at the interface
      
    Jx0 = Jx_l(2:end-1).*Ax_l(2:end-1) - Jx_l(3:end).*Ax_r(2:end-1);        % Flux into each element * the boundary area [mol/s]

    % Fick's 2nd law    
    cljr = Jx0./Vx(2:end-1) * dt + cljr;    % dcl/dt
    clj(2:end-1) = cljr;                    % Updating cl
    
    % Updating cl and cr
    if j == 1
        D_Ts(:,1) = D_T;
        Ji(1) = Jx_l(I)*Ai;                 % Flux at the interface [mol/s]     
    end
    if j == jss-1
        waitbar(j/(length(t)-1),wb)
        ji = ji+1;
        jss = js(ji);
        cl(:,ji) = clj;
        D_Ts(:,ji) = D_T;
        Ji(ji) = Jx_l(I)*Ai;
    end
end
close(wb);

%% Final calculations
% To calculate H and Js
H = sum(cl.*Vx,1);

Js1 = -D_Ts(1,:) .* (cl(2,:)-cl(1,:))/dx * As1;   % Flux at the surfaces [mol/s]
Js2 = -D_Ts(end,:) .* (cl(end,:)-cl(end-1,:))/dx * As2;

end

