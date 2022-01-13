%%%%%%% Drawing an energy landscape
% Inputs:
% [1,2]     [Left, Right]
% n         Numbers of normal atomic sites to each side of the interface
% Q         Energy barrier in the bulk of each crystal
% DeltaE_21 Difference in energy level between both crystals
% Eb        Energy barrier at each side of the interface
% DeltaE_1t Difference in energy level between crystal 1 and the interfacial trap site
%
% Output:
% f         Figure

function [f] = energy_landscape_plot(n,Q1,Q2,DeltaE_21,Eb1,DeltaE_1t,Eb2)

%%% No interface trap
if isempty(DeltaE_1t)
    dx = pi/(2*n)/1000;
    Eb2 = Q1+Eb1-Q2-DeltaE_21;
    x1 = 0:dx:(2*n-1)/(2*n)*pi;
    x2 = (2*n+1)/(2*n)*pi:dx:2*pi;
    xI1 = (2*n-1)/(2*n)*pi:dx:pi;
    xI2 = pi:dx:(2*n+1)/(2*n)*pi;

    y1 = Q1/2*cos(2*n*x1)+Q1/2;
    yI1 = (Q1+Eb1)/2*cos(2*n*xI1)+(Q1+Eb1)/2;
    yI2 = (Q2+Eb2)/2*cos(2*n*xI2)+(Q1+Eb1)-(Q2+Eb2)/2;
    y2 = Q2/2*cos(2*n*x2)+(Q1+Eb1)-(Q2+Eb2)+(Q2)/2;

    f = figure;
    plot(x1,y1,x2,y2,xI1,yI1,xI2,yI2,'Color','k')
    xlim([-1,2*pi+1])
    ylim([min([y1,y2,yI1,yI2])-1,max([y1,y2,yI1,yI2])+1])
    set(gca,'visible','off')

%%% With interface trap
else
    dx = pi/(2*n)/1000;
    x1 = 0:dx:(n-1)/n*pi;
    xI1 = (n-1)/n*pi:dx:(2*n-1)/(2*n)*pi;
    xT1 = (2*n-1)/(2*n)*pi:dx:pi;
    xT2 = pi:dx:(2*n+1)/(2*n)*pi;
    xI2 = (2*n+1)/(2*n)*pi:dx:(n+1)/n*pi;
    x2 = (n+1)/n*pi:dx:2*pi;
    y1 = -Q1/2*cos(2*n*x1)+Q1/2;
    yI1 = -(Q1+Eb1)/2*cos(2*n*xI1)+(Q1+Eb1)/2;
    yT1 = -(DeltaE_1t+Q1+Eb1)/2*cos(2*n*xT1)+(Q1+Eb1)-(DeltaE_1t+Q1+Eb1)/2;
    yT2 = -(DeltaE_1t+DeltaE_21+Q2+Eb2)/2*cos(2*n*xT2)+(DeltaE_1t+DeltaE_21+Q2+Eb2)/2-DeltaE_1t;
    yI2 = -(Q2+Eb2)/2*cos(2*n*xI2)+DeltaE_21+(Q2+Eb2)/2;
    y2 = -Q2/2*cos(2*n*x2)+DeltaE_21+Q2/2;

    f = figure;
    plot(x1,y1,x2,y2,xI1,yI1,xI2,yI2,xT1,yT1,xT2,yT2,'Color','k')
    xlim([-1,2*pi+1])
    ylim([min([y1,y2,yI1,yI2,yT1,yT2])-1,max([y1,y2,yI1,yI2,yT1,yT2])+1])
    set(gca,'visible','off')
end

end