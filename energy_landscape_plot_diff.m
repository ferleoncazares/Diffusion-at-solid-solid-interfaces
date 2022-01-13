%%%%%%% Drawing an energy landscape - diffuse interface
% Inputs:
% [1,2]     [Left, Right]
% n         Numbers of normal atomic sites to each side of the interface
% a         Separation between atomic sites in the bulk of each crystal
% Q         Energy barrier in the bulk of each crystal
% DeltaE_21 Difference in energy level between both crystals
% at        Separation between interfacial atomic sites
% DeltaEt1  Difference in energy level for each trap site
% Ebt       Height of the energy barrier between each interfacial site
%
% Output:
% f         Figure

function [f] = energy_landscape_plot_diff(n,a1,a2,Q1,Q2,DeltaE_21,at,DeltaEt1,Ebt)

dx = pi/1000*min(a1,a2);
xtu = cell(1,length(at));
ytu = xtu;
xtd = xtu;
ytd = xtu;

x1 = 0:dx:n*a1;
y1 = -Q1/2*cos(2*x1*pi/a1)+Q1/2;
xend = x1(end);
yend = y1(end);
xb1 = (0:dx:a1/2)+ xend;
yb1 = -Ebt(1)/2*cos(2*(xb1-xend)*pi/a1) + Ebt(1)/2+yend;
xend = xb1(end);
yend = yb1(end);
x = [x1,xb1];
y = [y1,yb1];

for i = 1:length(DeltaEt1)   
    xtd{i} = (0:dx:at(i)/2)+ xend;
    ytd{i} = (Ebt(i)-DeltaEt1(i))/2*cos(2*(xtd{i}-xend)*pi/at(i)) - (Ebt(i)-DeltaEt1(i))/2+yend;
    xend = xtd{i}(end);
    yend = ytd{i}(end);
    x = [x,xtd{i}];
    y = [y,ytd{i}];
    
    xtu{i} = (0:dx:at(i)/2)+ xend;
    ytu{i} = -(Ebt(i+1)-DeltaEt1(i))/2*cos(2*(xtu{i}-xend)*pi/at(i)) + (Ebt(i+1)-DeltaEt1(i))/2+yend;
    xend = xtu{i}(end);
    yend = ytu{i}(end);
    x = [x,xtu{i}];
    y = [y,ytu{i}];
end
xb2 = (0:dx:a2/2)+ xend;
yb2 = (Ebt(end)-DeltaE_21)/2*cos(2*(xb2-xend)*pi/a2) - (Ebt(end)-DeltaE_21)/2+yend;
xend = xb2(end);
yend = yb2(end);
x2 = (0:dx:n*a2)+ xend;
y2 = -Q2/2*cos(2*(x2-xend)*pi/a2)+Q2/2+yend;

x = [x,xb2,x2];
y = [y,yb2,y2];

f = figure;
plot(x,y,'Color','k')
set(gca,'visible','off')

end