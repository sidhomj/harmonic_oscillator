%The following code was used to analyze the dimensional system.
function main(tim,ch,range)

K1=10; K2=10; q1=1e-6; q2=1e-6; m1=1; m2=1;R1=1;R2=1;Ke= 8.987e9;
options=odeset('RelTol',1e-12,'Stats','on');
Xo=[1;0;pi;-1;1;0;0;-.5];
tspan=[0,tim];

tic
[t,X]=ode113(@Testfunction,tspan,Xo,options);
toc

[X3,X1]=pol2cart(X(:,3),X(:,1));
[X7,X5]=pol2cart(X(:,7),X(:,5));

if ch==1
p1=plot(X3,X1);
set(p1,'Color','blue');
hold on;
p2=plot(X7,X5);
set(p2,'Color','red');

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
% 
numtimes=3;
fps=10;
movie(M,numtimes,fps)

    
elseif ch==3;
dis=(sqrt(power(X(:,1),2)-2.*X(:,1).*X(:,5).*cos(X(:,3)-X(:,7))+power(X(:,5),2)));
    eb1=(1/2)*m1.*power(X(:,1),2).*power(X(:,4),2)+(1/2)*m1.*power(X(:,2),2)+(1/2)*K1.*power(X(:,1)-R1,2)+(Ke*q1*q2)./(sqrt(power(X(:,1),2)-2.*X(:,1).*X(:,5).*cos(X(:,3)-X(:,7))+power(X(:,5),2)));
    p2=plot(t,eb1,'*');
    set(p2,'Color','c');
    hold on;
    eb2=(1/2)*m2.*power(X(:,5),2).*power(X(:,8),2)+(1/2)*m2.*power(X(:,6),2)+(1/2)*K2.*power(X(:,5)-R2,2)+(Ke*q1*q2)./(sqrt(power(X(:,1),2)-2.*X(:,1).*X(:,5).*cos(X(:,3)-X(:,7))+power(X(:,5),2)));
    p3=plot(t,eb2);
    set(p3,'Color','b');
    hold on;
    total=eb1+eb2;
    p1=plot(t,total);
    maxi=max(total);
    set(p1,'Color','red');
    axis([0 tim 0 (1.1)*maxi]);
    xlabel('Time');
    ylabel('Energy');
    legend('Energy Ball 1','Energy Ball 2','Energy Total');
    
    figure
    eb1rot=(1/2)*m1.*power(X(:,1),2).*power(X(:,4),2);
    eb1trans=(1/2)*m1*power(X(:,2),2);
    eb1spring=(1/2)*K1.*power(X(:,1)-R1,2);
    eb1elec=(Ke*q1*q2)./dis;
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
     xlabel('Time (s)');
    ylabel('Energy (J)');
    
    figure
    eb2rot=(1/2)*m2.*power(X(:,5),2).*power(X(:,8),2);
    eb2trans=(1/2)*m2*power(X(:,6),2);
    eb2spring=(1/2)*K2.*power(X(:,5)-R2,2);
    eb2elec=(Ke*q1*q2)./dis;
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
    xlabel('Time (s)');
    ylabel('Energy (J)');
    
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
    title('Total Energy within the System');
     xlabel('Time (s)');
    ylabel('Energy (J)');
   

    elseif ch==5;
        
        p1=plot(t,X(:,4));
        title('ang speed ball 1')
        set(p1,'Color','red');
        hold on;
        p4=plot(t,X(:,8));
        set(p4,'Color','blue');
        title('ang speed ball 2');
        figure
        p2=plot(t,X(:,2));
        title('trans speed ball 1');
        hold on;
        p3=plot(t,X(:,6));
        title('trans speed ball 2');
        set(p2,'Color','blue');
        set(p3,'Color','red');
        


elseif ch==6;
    
    %Calculate angular Momentum of Ball 1
    angmom1=m1*power(X(:,1),2).*X(:,4);
    linmom1=m1*X(:,2);
    
    plot(t,angmom1);
    title('Angular Momentum of Ball 1');
    figure
    plot(t,linmom1);
    title('Linear Momentum of Ball 1');
    
    Calculate the angular Momentum of Ball 2
    angmom2=m2*power(X(:,5),2).*X(:,8);
    linmom2=m2*X(:,6);
    
    figure
    plot(t,angmom2);
    title('Angular Momentum of Ball 2');
    figure
    plot(t,linmom2);
    title('Linear Momentum of Ball 2');
    
    %Calculate the Total Linear Momentum
%     
    figure
    p1=plot(t,linmom1);
    set(p1,'Color','blue');
    hold on;
    p2=plot(t,linmom2);
    set(p2,'Color','red');
    hold on;
    p3=plot(t,linmom1+linmom2);
    set(p3,'Color','c');
    title('Linear Momentum');
    legend('Ball 1', 'Ball 2','Total');

    %Calculate the Total Angular Momentum
    
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
    xlabel('Time (s');
    ylabel('Momentum (N*m*s)');
    
     
    
else

end

return

