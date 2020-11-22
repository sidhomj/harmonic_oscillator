%The following code was used to analyze the dimensionless system.
function dimensionless(tim,ch,range);

mu=-.2;delta=1;lambda=1;zeta=1;

options=odeset('RelTol',1e-12,'Stats','on');
%x1=.943877; x5=x1;
%x1=power(2,1/3)/2+1; x5=1-power(2,1/3)/2;
Xo=[1;0;pi;.2;2;0;0;-.2];
tspan=[0,tim];

tic
[t,X]=ode113(@dimensionlessode,tspan,Xo,options);
toc

[X3,X1]=pol2cart(X(:,3),X(:,1));
[X4,X2]=pol2cart(X(:,4),X(:,2));
[X7,X5]=pol2cart(X(:,7),X(:,5));
[X8,X6]=pol2cart(X(:,8),X(:,6));

if ch==1
p1=plot(X3,X1);
set(p1,'Color','blue');
hold on;
p2=plot(X7,X5);
set(p2,'Color','red')
title('Trajectories of Masses');
xlabel('X-Position');
ylabel('Y-Position');

elseif ch==2;
     for j=1:5:size(X(:,1),1);
        p1=plot(X3(1:j),X1(1:j));
        hold on;
        set(p1,'Color','blue');
        hold on;
        p2=plot(X7(1:j),X5(1:j));
        set(p2,'Color','red');
        M(j)=getframe;
    end

numtimes=3;
fps=10;
movie(M,numtimes,fps)

elseif ch==3;
    
    dis=(sqrt(power(X(:,1),2)-2.*X(:,1).*X(:,5).*cos(X(:,3)-X(:,7))+power(X(:,5),2)));
    
    %Plot Energy Profile for Ball 1
     figure
    eb1rot=power(X(:,1),2).*power(X(:,4),2);
    eb1trans=power(X(:,2),2);
    eb1spring=power(X(:,1)-1,2);
    eb1elec=(2*mu)./dis;
    p5=plot(t,eb1rot);
    hold on;
    set(p5,'Color','red');
    p6=plot(t,eb1trans);
    hold on;
    set(p6,'Color','blue');
    p7=plot(t,eb1spring);
    hold on;
    set(p7,'Color','c');
    hold on;
    p9=plot(t,eb1elec);
    set(p9,'Color','m');
    hold on;
    p8=plot(t,eb1rot+eb1spring+eb1trans+eb1elec);
    hold on;
    set(p8,'Color','black');
    legend('Rotational Energy','Translational','Spring Energy','Electrical Energy','Total Energy');
    title('Ball 1 Energy');
    axis([0 tim 0 (1.1)*max(eb1rot+eb1spring+eb1trans+eb1elec)])
    
    %Plot Energy for Ball 2
     figure
    eb2rot=(power(X(:,5),2).*power(X(:,8),2))/lambda;
    eb2trans=(power(X(:,6),2))/lambda;
    eb2spring=(power(X(:,5)-(1/delta),2))/zeta;
    eb2elec=(2*mu)./dis;
    p5=plot(t,eb2rot);
    hold on;
    set(p5,'Color','red');
    p6=plot(t,eb2trans);
    hold on;
    set(p6,'Color','blue');
    p7=plot(t,eb2spring);
    hold on;
    set(p7,'Color','c');
    hold on;
    p9=plot(t,eb2elec);
    set(p9,'Color','m');
    hold on;
    p8=plot(t,eb2rot+eb2spring+eb2trans+eb2elec);
    hold on;
    set(p8,'Color','black');
    legend('Rotational Energy','Translational','Spring Energy','Electrical Energy','Total Energy');
    title('Ball 2 Energy');
    axis([0 tim 0 (1.1)*max(eb2rot+eb2spring+eb2trans+eb2elec)])
    
    figure
    eball1=eb1rot+eb1spring+eb1trans+eb1elec;
    eball2=eb2rot+eb2spring+eb2trans+eb2elec;
    p1=plot(t,eball1);
    set(p1,'Color','blue');
    hold on;
    p2=plot(t,eball2);
    set(p2,'Color','red');
    hold on;
    p3=plot(t,eball1+eball2);
    set(p3,'Color','black');
    legend('Ball 1','Ball 2','Total');
    title('Total Energy in System');

elseif ch==4;
    
     angmom1=power(X(:,1),2).*X(:,4);
     angmom2=(power(X(:,5),2).*X(:,8))/lambda;
     
     figure
    p1=plot(t,angmom1);
    set(p1,'Color','blue');
    hold on;
    p2=plot(t,angmom2);
    set(p2,'Color','red');
    hold on;
    p3=plot(t,angmom1+angmom2);
    set(p3,'Color','c');
    title('Angular Momentum');
    legend('Ball 1', 'Ball 2','Total');
    
elseif ch==5;
    figure
    p1=plot(t,X(:,1));
    title('Ball 1 Position');
    xlabel('Time');
    ylabel('Position');
    
    figure
    p2=plot(t,X(:,2));
    title('Ball 1 Radial Velocity');
    xlabel('Time');
    ylabel('Radial Velocity');
    
    figure
    p3=plot(t,X(:,3));
    title('Ball 1 Theta');
    
    figure
    p4=plot(t,X(:,4));
    title('Ball 1 Omega');
    xlabel('Time');
    ylabel('Angular Velocity');
    
else

end

return