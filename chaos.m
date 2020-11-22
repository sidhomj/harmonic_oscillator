%The following code was used to analyze the presence of sensitive
%dependence on initial conditions
function chaos(tim1,tim2,k,p,tstep);

% mu=.1;delta=1;lambda=1;zeta=1;

options=odeset('RelTol',1e-13,'Stats','on');
%x1=.943877; x5=x1;
%x1=power(2,1/3)/2+1; x5=1-power(2,1/3)/2;
x1=2; x5=2;
X01=[x1;1;pi;1;x5;1;0;-1];
tspan=[0:tstep:tim2];

tic
[t1,X1]=ode113(@dimensionlessode,tspan,X01,options);
toc

X02= X01+p;
tic
[t2,X2]=ode113(@dimensionlessode,tspan,X02,options);
toc


%Calculate  liaponov exponent for x
delta=log(abs(X1(:,k)-X2(:,k)));
init=(tim1/tstep);
fin=(tim2/tstep);
plot(t1(init:fin,1),delta(init:fin,1));
xlabel('t');
ylabel('ln(delta)');
title('Exponential Divergence of Nearby Trajectories (x1)');

figure
p1=plot(t1,X1(:,k));
set(p1,'Color','blue');
hold on;
p2=plot(t1,X2(:,k));
set(p2,'Color','red');
xlabel('t');
ylabel('Values of X01,X02');
title('Exponential Divergence of Nearby Trajectories');
legend('X01','X02');


end
