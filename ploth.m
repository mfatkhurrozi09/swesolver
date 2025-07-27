%plot time series

clear all; close all; clc;

load uoutput.txt;
s1=uoutput;

nt=400;
nx=201;

for m=1:nx
    xx(m)=s1(m,2);
    eta(m)=s1(m,5);
    dep(m)=-s1(m,3);
end

%dep=dep-eta;
%dep2=-dep;

l=0;
for m=1:nt
    for n=1:nx
        elev1(n,m)=s1(l+n,4);
    end
    l=l+n;
end

for t=1:nt

    figure (1);

    %subplot(2,1,1)
    plot(xx,elev1(:,t))
    xlabel('position y (m)')
    ylabel('eta (m)')
    ylim([-6. 1.0]);
    xlim([5 200]);
    %grid on

    hold on
    plot(xx,dep,'k')

    title('Muka air ')

    set(gcf,'position',[50,100,1100,170])
    hold off
    %pause (0.02)

end

