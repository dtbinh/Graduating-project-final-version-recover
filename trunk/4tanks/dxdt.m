function DXDT = dxdt(x)

% states
%x1 = h1 \in [0.1; 20] [cm]
%x2 = h2 \in [0.1; 20] [cm]
%x3 = h3 \in [0.1; 20] [cm]
%x4 = h4 \in [0.1; 20] [cm]

%parametros da planta (constantes do modelo)
global A_ k g gamma a v;

var_degrau = 1;

if var_degrau == 1
    v1 = f_input(t,v(1));
    v2 = f_input(t,v(2));
else
    v1 = v(1);
    v2 = v(2);
end

%function
dx1dt = (-a(1)*sqrt(2*g*x(1)) + a(3)*sqrt(2*g*x(3)) + gamma(1)...
            *k(1)*v1)/A_;
dx2dt = (-a(2)*sqrt(2*g*x(2)) + a(4)*sqrt(2*g*x(4)) + gamma(2)...
            *k(2)*v2)/A_;
dx3dt = (-a(3)*sqrt(2*g*x(3)) + ((1 - gamma(2))*k(2)*v2))/A_;
dx4dt = (-a(4)*sqrt(2*g*x(4)) + ((1 - gamma(1))*k(1)*v1))/A_;

DXDT = [dx1dt dx2dt dx3dt dx4dt]';
%%
function fe = f_input(t,x)
% disturbio de entrada
% x = [x1 x2 t1 t2]
if t < 100 || t > 300
    fe = x;
else
    fe = x*1.2;
end





