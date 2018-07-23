function X = dxdt_fsolve(X0)

if nargin == 0
    X0 = [0 0 0 0];
end

X = fsolve(@dxdt,X0);

function DXDT = dxdt(x) 

global A_ k g gamma a v;

v(1) = 1;
v(2) = 1;

dx1dt = (-a(1)*sqrt(2*g*x(1)) + a(3)*sqrt(2*g*x(3)) + gamma(1)...
            *k(1)*v(1))/A_;
dx2dt = (-a(2)*sqrt(2*g*x(2)) + a(4)*sqrt(2*g*x(4)) + gamma(2)...
            *k(2)*v(2))/A_;
dx3dt = (-a(3)*sqrt(2*g*x(3)) + ((1 - gamma(2))*k(2)*v(2)))/A_;
dx4dt = (-a(4)*sqrt(2*g*x(4)) + ((1 - gamma(1))*k(1)*v(1)))/A_;

DXDT = [dx1dt dx2dt dx3dt dx4dt]';