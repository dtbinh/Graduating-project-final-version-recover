function DXDT = dxdt(t,x) 

dx1dt = -2*x(1)+4*x(2);
dx2dt = (-1-(20*(1-sin(x(1)))/2))*x(1)-2*x(2);

DXDT = [dx1dt dx2dt]';


