%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples for the multiple scenarios of solid-solid interfaces developed
% in the paper:

    % F.D. León-Cázares, E.I. Galindo-Nava, General model for the 
    % kinetics of solute diffusion at solid-solid interfaces.
    % Physical Review Materials, 5 (2021) 123802.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%%%%% Drawing energy landscapes
%{

% Atomic jump
n = 4;
dx = pi/(2*n)/1000;
Q1 = 3;
Q2 = 1;
Eb2 = 0;
DeltaE = 3.5;
Eb1 = DeltaE+Q2+Eb2-Q1;
x1 = (2*n-1.2)/(2*n)*pi:dx:(2*n-1)/(2*n)*pi;
x2 = (2*n+1)/(2*n)*pi:dx:(2*n+1.3)/(2*n)*pi;
xI1 = (2*n-1)/(2*n)*pi:dx:pi;
xI2 = pi:dx:(2*n+1)/(2*n)*pi;
y1 = Q1/2*cos(2*n*x1)+Q1/2;
y2 = Q2/2*cos(2*n*x2)+(Q1+Eb1)-(Q2+Eb2)+(Q2)/2;
yI1 = (Q1+Eb1)/2*cos(2*n*xI1)+(Q1+Eb1)/2;
yI2 = (Q2+Eb2)/2*cos(2*n*xI2)+(Q1+Eb1)-(Q2+Eb2)/2;
figure;
plot(x1,y1,x2,y2,xI1,yI1,xI2,yI2,'Color','k')
xlim([-1,2*pi+1])
ylim([min([yI1,yI2])-1,max([yI1,yI2])+1])
set(gca,'visible','off')
text(0,0,'Atomic jump')

% Homogeneous
Q1 = 3;
Q2 = 3;
Eb2 = 0;
DeltaE = 0;
Eb1 = DeltaE+Q2+Eb2-Q1;
energy_landscape_plot(4,Q1,Q2,DeltaE,Eb1,[],[]);
text(0,-0.4,'Homogeneous')

% Point trap
Q1 = 3;
Q2 = 3;
Eb1 = 0;
Eb2 = 0;
DeltaE_21 = 0;
DeltaE_1t = 3;
energy_landscape_plot(3,Q1,Q2,DeltaE_21,Eb1,DeltaE_1t,Eb2);
text(0,-0.6,'Point trap')

% Simplest bicrystal
Q1 = 3;
Q2 = 1;
Eb2 = 0;
DeltaE = 3.5;
Eb1 = DeltaE+Q2+Eb2-Q1;
energy_landscape_plot(4,Q1,Q2,DeltaE,Eb1,[],[]);
text(0,-0.4,'Simplest bicrystal')

% Bicrystal - Eb
Q1 = 3;
Q2 = 1;
Eb2 = 1.5;
DeltaE = 3.5;
Eb1 = DeltaE+Q2+Eb2-Q1;
energy_landscape_plot(4,Q1,Q2,DeltaE,Eb1,[],[]);
text(0,-0.4,'Bicrystal - energy barrier')

% Trap - no Eb
Q1 = 3;
Q2 = 1;
Eb1 = 0;
Eb2 = 0;
DeltaE_21 = 3.5;
DeltaE_1t = 2;
energy_landscape_plot(4,Q1,Q2,DeltaE_21,Eb1,DeltaE_1t,Eb2);
text(-0.2,-0.4,'Bicrystal - monolayer trap')

% Trap - with Ebs
Q1 = 3;
Q2 = 1;
Eb1 = 1.5;
Eb2 = 1;
DeltaE_21 = 3.5;
DeltaE_1t = 2;
energy_landscape_plot(4,Q1,Q2,DeltaE_21,Eb1,DeltaE_1t,Eb2);
text(-0.2,-2.4,'Bicrystal - monolayer trap with energy barriers')

% Negative trap - with Ebs
Q1 = 3;
Q2 = 1;
Eb1 = 3;
Eb2 = 2;
DeltaE_21 = 3.5;
DeltaE_1t = -5;
energy_landscape_plot(4,Q1,Q2,DeltaE_21,Eb1,DeltaE_1t,Eb2);
text(-0.2,-0.4,'Bicrystal - negative trap')

% Diffuse interface
Q1 = 3;
Q2 = 1;
DeltaE_21 = 3.5;
                        % Trap parameters
DeltaEt1max1 = -5;
DeltaEt1max2 = -5;
DeltaEt1 = [[0.05,0.15,0.4,0.65,1]*DeltaEt1max1,DeltaE_21+[0.6,0.2]*DeltaEt1max2];
DeltaEbt1 = [Q1,DeltaEt1(1)+Q1,DeltaEt1(2)+Q1,DeltaEt1(3)+Q1,DeltaEt1(4)+Q1,...
             DeltaEt1(6)+Q2,DeltaEt1(7)+Q2,DeltaE_21+Q2];

energy_landscape_plot_diff(5,1,1,Q1,Q2,DeltaE_21,ones(size(DeltaEt1)),DeltaEt1,DeltaEbt1);
text(-0.2,-0.4,'Bicrystal - diffuse trap')

% Inhomogeneous trap
Q1 = 3;
Q2 = 1;
Eb1 = [0,0];
Eb2 = [0,0];
DeltaE_21 = 3.5;
DeltaE_1t = [2,1.5];
energy_landscape_plot(4,Q1,Q2,DeltaE_21,Eb1(1),DeltaE_1t(1),Eb2(1));
text(-0.2,-0.4,'Bicrystal - monolayer trap')
[f,ff] = deal(cell(1,length(DeltaE_1t))); % Energy landscape
for i = 1:length(DeltaE_1t)
    f{i} = energy_landscape_plot(4,Q1,Q2,DeltaE_21,Eb1(i),DeltaE_1t(i),Eb2(i));
    if i ~= 1
        ff{i} = findobj(f{i},'type','line');
        copyobj(ff{i},findobj(f{1},'type','axes'));
        close(f{i});
    end
end
ylim([-max(DeltaE_1t),max(max(DeltaE_21+Q2+Eb2),max(Q1+Eb1))])
text(-0.2,-2.4,'Inhomogeneous monolayer trap')


%}

%%%%% Boundary conditions: constant concentrations at the surfaces (Section III B.)
%%%%% 
%%% 1. Energy barrier
% Ni - fcc-octahedral site - Sigma3(111)[-110] GB
%{

c1_ini = 0;                             % Initial concentration in the lattice [mol/m3]
c2_ini = 0;
cs1 = 1;                                % Concentration at free surface [mol/m3]
cs2 = 0;

aspacing = 3.53e-10;                    % Lattice parameter [m]
DeltaE_21 = 0;                          % Energy difference [J/mol]
D01 = 7e-7;                             % {Louthan1975} Coefficient D0 [m2/s]
Q1 = 0.37*1.60218e-19*6.022e23;         % Enery barrier [J/mol]
S1 = 4/aspacing^3 / 6.022e23;           % Site density [mol/m^3]    
D02 = D01;          
Q2 = Q1;
S2 = S1;
Eb2 = (0.55-0.37)*1.60218e-19*6.022e23; % INTERFACIAL ENERGY BARRIER [J/mol]
Eb1 = DeltaE_21-Q1+Q2+Eb2;

a = aspacing/2;                         % Layer spacing [m] (1/4*lattice parameter)
xs = 5001*a;                            % Size of the simulated space [m]
xi = xs/2;                              % Interface position [m]
dx = xs/21;

T0 = 300 + 273;                         % Temperature  [K]
tend = 1.5e-3;                          % Time [s]
dt1 = 0.45*a^2/(D01*exp(-Q1/(8.31446261815324*T0)));
dt = dt1*(dx/a).^2;

coord = 1;              % Coordinate system (1 -> Cartesian, 2 -> Cylindrical, 3 -> Spherical)
L = [10e-6,10e-6];      % Sample dimensions (Cartesian [L,w], Cylindrical [L,~] or Spherical [~])
nsvs = 200;             % Number of saving points

% Simulation
[x,t,ts,cl,Js1,Js2,Ji,H,DeltaE] = traps1D_perm(xs,xi,T0,dx,dt,tend,S1,D01,Q1,Eb1,S2,D02,Q2,Eb2,a,cs1,cs2,c1_ini,c2_ini,coord,L,nsvs);

% Plots
energy_landscape_plot(4,Q1,Q2,DeltaE_21,Eb1,[],[]);

tI = [1,5,21,201];
figure
plot(x*1e6,cl(:,tI(1)),'k-o','LineWidth',1.5)
hold on
plot(x*1e6,cl(:,tI(2)),'r-o','LineWidth',1.5)
plot(x*1e6,cl(:,tI(3)),'b-o','LineWidth',1.5)
plot(x*1e6,cl(:,tI(4)),'g-o','LineWidth',1.5)
set(gca,'fontsize',16,'FontName','Times')
xlim([0,x(end)*1e6])
ylim([0,1.2*max(max(cl))])
xlabel('$x$ \textrm{[$\mu$m]}','Interpreter','latex')
ylabel('$c(x,t)$ \textrm{[mol/m$^3$]}','Interpreter','latex')
legend({[num2str(ts(tI(1))*1000),' ms'],[num2str(ts(tI(2))*1000),' ms'],[num2str(ts(tI(3))*1000),' ms'],[num2str(ts(tI(4))*1000),' ms']},'FontSize',16,'Location','Northeast')
grid off
box on

figure
plot(ts,[Js1;Ji;Js2]/(L(1)*L(2)),'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Flux [mol/m^2s]')
legend({'Left boundary','Interface','Right boundary'},'FontSize',16,'Location','Northeast')
set(gca,'fontsize',16)
grid off
box on



%}

%%% 2. Monolayer trap
% Du 2011
% Fe - bcc-tetrahedral site - Sigma3 GB
%{

c1_ini = 0;                             % Initial concentration in the lattice [mol/m3]
c2_ini = 0;
cs1 = 1;                                % H concentration at free surface [mol/m3]
cs2 = 0;

aspacing = 2.85e-10;                    % Lattice parameter [m]
a = aspacing/4;                         % Layer spacing [m] (1/4*lattice parameter)
DeltaE_21 = 0;                          % Energy difference [J/mol]
D01 = 1.379e-8;                         % {He2017} Coefficient D0 [m2/s]
Q1 = 0.088*1.60218e-19*6.022e23;        % [Du 2012] Enery barrier [J/mol]
S1 = 12/aspacing^3 / 6.022e23;          % Site density [mol/m^3]    
Eb1 = 0;                                % ASSUMED - Interfacial energy barrier [J/mol]
D02 = D01;          
Q2 = Q1;
S2 = S1;
Eb2 = 0;
DeltaE_1t = 0.43*1.60218e-19*6.022e23;  % Trap parameters
St = 1/(2.44e-10*4e-10*a) / 6.022e23;       % From the area of the supercall and thickness a
D0t = D01;
Ebt1 = DeltaE_1t+Q1+Eb1;
Ebt2 = DeltaE_1t+DeltaE_21+Q2+Eb2;

xs = 5000*a;                            % Size of the simulated space [m]
xi = xs/2;                              % Interface position [m]
dx = xs/20;

T0 = 300 + 273;                         % Temperature  [K]
tend = 2e-4;                            % Time [s]
dt1 = 0.45*a^2/(D01*exp(-Q1/(8.31446261815324*T0)));
dt = dt1*(dx/a).^2;

coord = 1;              % Coordinate system (1 -> Cartesian, 2 -> Cylindrical, 3 -> Spherical)
L = [10e-6,10e-6];      % Sample dimensions (Cartesian [L,w], Cylindrical [L,~] or Spherical [~])
nsvs = 200;             % Number of saving points

% Simulation
[xr,x,t,ts,clr,cl,civ,Js1,Js2,Ji,H] = traps1D_perm_t(xs,xi,T0,dx,dt,tend,S1,D01,Q1,Eb1,S2,D02,Q2,Eb2,St,DeltaE_21,DeltaE_1t,D0t,a,cs1,cs2,c1_ini,c2_ini,coord,L,nsvs);

% Plots
energy_landscape_plot(4,Q1,Q2,DeltaE_21,Eb1,DeltaE_1t,Eb2);

tI = [1,4,201];
figure
plot(xr*1e6,clr(:,tI(1)),'k-o','LineWidth',1.5)
hold on
plot(xr*1e6,clr(:,tI(2)),'r-o','LineWidth',1.5)
plot(xr*1e6,clr(:,tI(3)),'b-o','LineWidth',1.5)
set(gca,'fontsize',16,'FontName','Times')
xlim([0,x(end)*1e6])
ylim([0,1.2*max(cs1,cs2)])
xlabel('$x$ \textrm{[$\mu$m]}','Interpreter','latex')
ylabel('$c(x,t)$ \textrm{[mol/m$^3$]}','Interpreter','latex')
legend({[num2str(ts(tI(1))*1000),' ms'],[num2str(ts(tI(2))*1000),' ms'],[num2str(ts(tI(3))*1000),' ms']},'FontSize',16,'Location','Northeast')
grid off
box on

figure          % ct(t)
plot(ts*1e3,cl(end/2+0.5,:),'LineWidth',1.5)
set(gca,'fontsize',16,'FontName','Times')
xlabel('Time [ms]')
ylabel('$c_t(t)$ \textrm{[mol/m$^3$]}','Interpreter','latex','fontsize',16)
grid off
box on

figure          % J
plot(ts,[Js1;Ji;Js2]/(L(1)*L(2)),'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Flux [mol/m^2s]')
legend({'Left boundary','Interface (left)','Interface (right)','Right boundary'},'FontSize',16,'Location','Northeast')
set(gca,'fontsize',16)
grid off
box on


%}

%%% 3. Diffuse trap
% Du 2011
% Fe - bcc-tetrahedral site - Sigma5 GB
%{

c1_ini = 0;                             % Initial concentration in the lattice [mol/m3]
c2_ini = 0;
cs1 = 1;                                % H concentration at free surface [mol/m3]
cs2 = 0;

aspacing = 2.85e-10;                    % Lattice parameter [m]
a = aspacing/4;                         % Layer spacing [m] (1/4*lattice parameter)
DeltaE_21 = 0;                          % Energy difference [J/mol]
D01 = 1.379e-8;                             % {He2017} Coefficient D0 [m2/s]
Q1 = 0.088*1.60218e-19*6.022e23;          % Enery barrier [J/mol]
S1 = 12/aspacing^3 / 6.022e23;          % Site density [mol/m^3]    
Eb1 = 0;                                    % ASSUMED - Interfacial energy barrier [J/mol]
D02 = D01;          
Q2 = Q1;
S2 = S1;
Eb2 = 0;
                        % Trap parameters
DeltaEt1 = [-0.13,-0.5,-0.13]*1.60218e-19*6.022e23;
D0t = [D01,D01,D01];
St = [1,1,1]*1/(4.67e-10*4.67e-10*a)/6.022e23;    % From the area of the supercall and thickness a
DeltaEbt1 = [0.07,-0.07,-0.07,0.07]*1.60218e-19*6.022e23;

xs = 5000*a;                            % Size of the simulated space [m]
xi = xs/2;                              % Interface position [m]
dx = xs/20;

T0 = 300 + 273;         % Temperature  [K]
tend = 2e-4;            % Time [s]
dt1 = 0.45*a^2/(D01*exp(-(0.06*1.60218e-19*6.022e23)/(8.31446261815324*T0)));
dt = dt1*(dx/a).^2;

coord = 1;              % Coordinate system (1 -> Cartesian, 2 -> Cylindrical, 3 -> Spherical)
L = [10e-6,10e-6];      % Sample dimensions (Cartesian [L,w], Cylindrical [L,~] or Spherical [~])
nsvs = 200;             % Number of saving points

% Simulation
[xr,x,t,ts,clr,cl,civ,Js1,Js2,Ji,H] = traps1D_perm_t_diff(xs,xi,T0,dx,dt,tend,S1,D01,Q1,S2,D02,Q2,DeltaE_21,St,DeltaEt1,D0t,DeltaEbt1,a,cs1,cs2,c1_ini,c2_ini,coord,L,nsvs);

% Plots
energy_landscape_plot_diff(4,1,1,Q1,Q2,DeltaE_21,ones(size(D0t)),DeltaEt1,DeltaEbt1);

tI = [1,4,201];
figure
plot(xr*1e6,clr(:,tI(1)),'k-o','LineWidth',1.5)
hold on
plot(xr*1e6,clr(:,tI(2)),'r-o','LineWidth',1.5)
plot(xr*1e6,clr(:,tI(3)),'b-o','LineWidth',1.5)
set(gca,'fontsize',16,'FontName','Times')
xlim([0,xr(end)*1e6])
ylim([0,1.2*max(cs1,cs2)])
yticks([0,0.5,1])
xlabel('$x$ \textrm{[$\mu$m]}','Interpreter','latex')
ylabel('$c(x,t)$ \textrm{[mol/m$^3$]}','Interpreter','latex')
legend({[num2str(ts(tI(1))*1000),' ms'],[num2str(ts(tI(2))*1000),' ms'],[num2str(ts(tI(3))*1000),' ms']},'FontSize',16,'Location','Northeast')
grid off
box on

figure          % ct(t)
plot(ts*1e3,civ(2,:),'LineWidth',1.5)
set(gca,'fontsize',16,'FontName','Times')
ylim([0,1.1*max(civ(2,:))])
xlabel('Time [ms]')
ylabel('$c_t(t)$ \textrm{[mol/m$^3$]}','Interpreter','latex')
grid off
box on

figure          % J
plot(ts,[Js1;Ji;Js2]/(L(1)*L(2)),'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Flux [mol/m^2s]')
legend({'Left boundary','Interface (left)','Interface (right)','Right boundary'},'FontSize',16,'Location','Northeast')
set(gca,'fontsize',16)
grid off
box on

%}

%%% 4. Inhomogeneous trap
% Di Stefano 2016
% Fe + TiC bicrystal interface 
%{

c1_ini = 0;                             % Initial concentration in the lattice [mol/m3]
c2_ini = 0;
cs1 = 1;                                % H concentration at free surface [mol/m3]
cs2 = 0;

aspacingFe = 2.85e-10;                  % Lattice parameter [m]
DeltaE_21 = 0.68*1.60218e-19*6.022e23;  % Energy difference [J/mol]
D01 = 1.379e-8;                         % {He2017} Coefficient D0 [m2/s]
Q1 = 0.09*1.60218e-19*6.022e23;         % Enery barrier [J/mol]
S1 = 12/aspacingFe^3 / 6.022e23;        % Site density [mol/m^3]
aspacingTiC = 4.34e-10;
D02 = D01;                                  % MISSING
Q2 = 0.47*1.60218e-19*6.022e23;
S2 = 32/aspacingTiC^3 / 6.022e23;
            % Trap parameters 
fd = 1/16;                              % Fraction of dislocation (corresponds to a separation of 4.6 nm, similar to the 4.2 experimental value of [Wei2006])
St = [S1*(1-fd),S1*fd]; 
DeltaE_1t = [0.32,0.49]*1.60218e-19*6.022e23;
D0t = zeros(size(St))+D01;                  % ASSUMED
Eb1 = zeros(size(St));
Eb2 = zeros(size(St));
Ebt1 = DeltaE_1t+Q1+Eb1;                    % ASSUMED
Ebt2 = DeltaE_1t+DeltaE_21+Q2+Eb2;          % ASSUMED

a = aspacingFe/4;                       % Layer spacing [m] (1/4*lattice parameter)
xs = 5000*a;                            % Size of the simulated space [m]
xi = xs/2;                              % Interface position [m]
dx = xs/20;

T0 = 300 + 273;                         % Temperature  [K]
tend = 5e-4;                            % Time [s]
dt1 = 0.45*a^2/(D01*exp(-Q1/(8.31446261815324*T0)));
dt = dt1*(dx/a).^2;

coord = 1;              % Coordinate system (1 -> Cartesian, 2 -> Cylindrical, 3 -> Spherical)
L = [10e-6,10e-6];      % Sample dimensions (Cartesian [L,w], Cylindrical [L,~] or Spherical [~])
nsvs = 200;             % Number of saving points

% Simulation
[xr,x,t,ts,clr,cl,ct,civ,Js1,Js2,Ji,H,k] = traps1D_perm_t_inh(xs,xi,T0,dx,dt,tend,S1,D01,Q1,Eb1,S2,D02,Q2,Eb2,St,DeltaE_21,DeltaE_1t,D0t,a,cs1,cs2,c1_ini,c2_ini,coord,L,nsvs);

% Plots
[f,ff] = deal(cell(1,length(DeltaE_1t))); % Energy landscape
for i = 1:length(DeltaE_1t)
    f{i} = energy_landscape_plot(4,Q1,Q2,DeltaE_21,Eb1(i),DeltaE_1t(i),Eb2(i));
    if i ~= 1
        ff{i} = findobj(f{i},'type','line');
        copyobj(ff{i},findobj(f{1},'type','axes'));
        close(f{i});
    end
end
ylim([-max(DeltaE_1t),max(max(DeltaE_21+Q2+Eb2),max(Q1+Eb1))])

tI = [1,3,201];
figure
plot(xr*1e6,clr(:,tI(1)),'k-o','LineWidth',1.5)
hold on
plot(xr*1e6,clr(:,tI(2)),'r-o','LineWidth',1.5)
plot(xr*1e6,clr(:,tI(3)),'b-o','LineWidth',1.5)
set(gca,'fontsize',16,'FontName','Times')
xlim([0,xr(end)*1e6])
ylim([0,1.2*max(cs1,cs2)])
xlabel('$x$ \textrm{[$\mu$m]}','Interpreter','latex')
ylabel('$c(x,t)$ \textrm{[mol/m$^3$]}','Interpreter','latex')
legend({[num2str(ts(tI(1))*1000),' ms'],[num2str(ts(tI(2))*1000),' ms'],[num2str(ts(tI(3))*1000),' ms']},'FontSize',16,'Location','Northeast')
grid off
box on

figure          % ct(t)
plot(ts*1e3,sum(ct),'LineWidth',1.5)
hold on
plot(ts*1e3,ct(1,:),'LineWidth',1.5)
plot(ts*1e3,ct(2,:),'LineWidth',1.5)
set(gca,'fontsize',16,'FontName','Times')
ylim([0,1.1*max(sum(ct))])
xlabel('Time [ms]')
ylabel('$c_t(t)$ \textrm{[mol/m$^3$]}','Interpreter','latex','fontsize',16)
legend({'Total','Trap site 1','Trap site 2'},'FontSize',16,'Location','Southeast')
grid off
box on

figure          % J
plot(ts,[Js1;Ji;Js2]/(L(1)*L(2)),'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Flux [mol/m^2s]')
legend({'Left boundary','Interface (left)','Interface (right)','Right boundary'},'FontSize',16,'Location','Northeast')
set(gca,'fontsize',16)
grid on
box on
    
%}

%%%%% Open system - BCs: Zero flux (left) and constant concentration (right)
%%% Monolayer trap (Made up values)
%{

c1_ini = 0;                             % Initial concentration in the lattice [mol/m3]
c2_ini = 0;
cs2 = 1;                                % H concentration at free surface [mol/m3]

aspacing = 2.85e-10;                    % Lattice parameter [m]
a = aspacing/4;                         % Layer spacing [m] (1/4*lattice parameter)
DeltaE_21 = 0.05*1.60218e-19*6.022e23;   % Energy difference [J/mol]
D01 = 1.379e-8;                         % Coefficient D0 [m2/s]
Q1 = 0.088*1.60218e-19*6.022e23;        % Enery barrier [J/mol]
S1 = 12/aspacing^3 / 6.022e23;          % Site density [mol/m^3]    
Eb1 = 0.2*1.60218e-19*6.022e23;         % Interfacial energy barrier [J/mol]
D02 = D01;          
Q2 = Q1;
S2 = S1;
Eb2 = 0.1*1.60218e-19*6.022e23;
DeltaE_1t = 0.43*1.60218e-19*6.022e23;  % Trap parameters
St = 1/(2.44e-10*4e-10*a) / 6.022e23;
D0t = D01;
Ebt1 = DeltaE_1t+Q1+Eb1;
Ebt2 = DeltaE_1t+DeltaE_21+Q2+Eb2;

xs = 5000*a;                            % Size of the simulated space [m]
xi = xs/2;                              % Interface position [m]
dx = xs/20;

T0 = 300 + 273;                         % Temperature  [K]
tend = 4e-4;                            % Time [s]
dt1 = 0.3*a^2/(D01*exp(-Q1/(8.31446261815324*T0)));
dt = dt1*(dx/a).^2;

coord = 3;              % Coordinate system (1 -> Cartesian, 2 -> Cylindrical, 3 -> Spherical)
L = [10e-6,10e-6];      % Sample dimensions (Cartesian [L,w], Cylindrical [L,~] or Spherical [~])
nsvs = 200;             % Number of saving points

% Simulation
[xr,x,t,ts,clr,cl,civ,Js1,Js2,Ji,H] = traps1D_open_t(xs,xi,T0,dx,dt,tend,S1,D01,Q1,Eb1,S2,D02,Q2,Eb2,St,DeltaE_21,DeltaE_1t,D0t,a,cs2,c1_ini,c2_ini,coord,L,nsvs);

% Plots
energy_landscape_plot(4,Q1,Q2,DeltaE_21,Eb1,DeltaE_1t,Eb2);

tI = [1,4,201];
figure
plot(xr*1e6,clr(:,tI(1)),'k-o','LineWidth',1.5)
hold on
plot(xr*1e6,clr(:,tI(2)),'r-o','LineWidth',1.5)
plot(xr*1e6,clr(:,tI(3)),'b-o','LineWidth',1.5)
set(gca,'fontsize',16,'FontName','Times')
xlim([0,x(end)*1e6])
ylim([0,1.2*max(max(clr([1,end],:)))])
xlabel('$x$ \textrm{[$\mu$m]}','Interpreter','latex')
ylabel('$c(x,t)$ \textrm{[mol/m$^3$]}','Interpreter','latex')
legend({[num2str(ts(tI(1))*1000),' ms'],[num2str(ts(tI(2))*1000),' ms'],[num2str(ts(tI(3))*1000),' ms']},'FontSize',16,'Location','Northeast')
grid off
box on

figure          % ct(t)
plot(ts*1e3,cl(end/2+0.5,:),'LineWidth',1.5)
set(gca,'fontsize',16,'FontName','Times')
xlabel('Time [ms]')
ylabel('$c_t(t)$ \textrm{[mol/m$^3$]}','Interpreter','latex','fontsize',16)
grid off
box on

figure          % J
plot(ts,[Ji;Js2]/(L(1)*L(2)),'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Flux [mol/m^2s]')
legend({'Interface (left)','Interface (right)','Right boundary'},'FontSize',16,'Location','Northeast')
set(gca,'fontsize',16)
grid off
box on

%}

%%%%% Closed system - BCs: Zero flux at both surfaces 
%%% Monolayer trap (Made up values)
%{

c1_ini = 1;                             % Initial concentration in the lattice [mol/m3]
c2_ini = 0;

aspacing = 2.85e-10;                    % Lattice parameter [m]
a = aspacing/4;                         % Layer spacing [m] (1/4*lattice parameter)
DeltaE_21 = -0.1*1.60218e-19*6.022e23;  % Energy difference [J/mol]
D01 = 1.379e-8;                         % Coefficient D0 [m2/s]
Q1 = 0.088*1.60218e-19*6.022e23;        % Enery barrier [J/mol]
S1 = 12/aspacing^3 / 6.022e23;          % Site density [mol/m^3]    
Eb1 = 0.4*1.60218e-19*6.022e23;         % Interfacial energy barrier [J/mol]
D02 = D01;          
Q2 = Q1;
S2 = S1;
Eb2 = 0.2*1.60218e-19*6.022e23;
DeltaE_1t = 0.43*1.60218e-19*6.022e23;  % Trap parameters
St = 1/(2.44e-10*4e-10*a) / 6.022e23;
D0t = D01;
Ebt1 = DeltaE_1t+Q1+Eb1;
Ebt2 = DeltaE_1t+DeltaE_21+Q2+Eb2;

xs = 5000*a;                            % Size of the simulated space [m]
xi = xs/2;                              % Interface position [m]
dx = xs/20;

T0 = 300 + 273;                         % Temperature  [K]
tend = 2e-4;                            % Time [s]
dt1 = 0.3*a^2/(D01*exp(-Q1/(8.31446261815324*T0)));
dt = dt1*(dx/a).^2;

coord = 1;              % Coordinate system (1 -> Cartesian, 2 -> Cylindrical, 3 -> Spherical)
L = [10e-6,10e-6];      % Sample dimensions (Cartesian [L,w], Cylindrical [L,~] or Spherical [~])
nsvs = 200;             % Number of saving points

% Simulation
[xr,x,t,ts,clr,cl,civ,Js1,Js2,Ji,H] = traps1D_closed_t(xs,xi,T0,dx,dt,tend,S1,D01,Q1,Eb1,S2,D02,Q2,Eb2,St,DeltaE_21,DeltaE_1t,D0t,a,c1_ini,c2_ini,coord,L,nsvs);

% Plots
energy_landscape_plot(4,Q1,Q2,DeltaE_21,Eb1,DeltaE_1t,Eb2);

tI = [1,4,201];
figure
plot(xr*1e6,clr(:,tI(1)),'k-o','LineWidth',1.5)
hold on
plot(xr*1e6,clr(:,tI(2)),'r-o','LineWidth',1.5)
plot(xr*1e6,clr(:,tI(3)),'b-o','LineWidth',1.5)
set(gca,'fontsize',16,'FontName','Times')
xlim([0,x(end)*1e6])
ylim([0,1.2*max(max(clr([1,end],:)))])
xlabel('$x$ \textrm{[$\mu$m]}','Interpreter','latex')
ylabel('$c(x,t)$ \textrm{[mol/m$^3$]}','Interpreter','latex')
legend({[num2str(ts(tI(1))*1000),' ms'],[num2str(ts(tI(2))*1000),' ms'],[num2str(ts(tI(3))*1000),' ms']},'FontSize',16,'Location','Northeast')
grid off
box on

figure          % ct(t)
plot(ts*1e3,cl(end/2+0.5,:),'LineWidth',1.5)
set(gca,'fontsize',16,'FontName','Times')
xlabel('Time [ms]')
ylabel('$c_t(t)$ \textrm{[mol/m$^3$]}','Interpreter','latex','fontsize',16)
grid off
box on

figure          % J
plot(ts,Ji/(L(1)*L(2)),'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Flux [mol/m^2s]')
legend({'Interface (left)','Interface (right)'},'FontSize',16,'Location','Northeast')
set(gca,'fontsize',16)
grid off
box on

%}



