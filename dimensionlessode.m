%Dimensionless ODE's
function [dx_dt]=dimensionlessode(t,x);

mu=-.2;delta=1;lambda=1;zeta=1;

dis=(x(1)^2-2*x(1)*x(5)*cos(x(3)-x(7))+x(5)^2)^(3/2);
A=x(1)-x(5)*cos(x(3)-x(7));
B=x(1)*x(5)*sin(x(3)-x(7));
C=x(5)-x(1)*cos(x(3)-x(7));
D=x(1)*x(5)*sin(x(7)-x(3));


dx_dt(1)=x(2);
dx_dt(2)=mu*(A/dis)+x(1)*(x(4)^2)-(x(1)-1);
dx_dt(3)=x(4);
dx_dt(4)=(mu/((x(1))^2))*(B/dis)-(2*x(2)*x(4))/x(1);
dx_dt(5)=x(6);
dx_dt(6)=mu*lambda*(C/dis)+x(5)*((x(8))^2)-(lambda/zeta)*(x(5)-(1/delta));
dx_dt(7)=x(8);
dx_dt(8)=((mu*lambda)/((x(5))^2))*(D/dis)-(2*x(6)*x(8))/x(5);

dx_dt=transpose(dx_dt);

return