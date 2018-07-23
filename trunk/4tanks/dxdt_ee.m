function DXDT = dxdt_ee(t, x)
% states
%x1 = h1 \in [0.1; 20] [cm]
%x2 = h2 \in [0.1; 20] [cm]
%x3 = h3 \in [0.1; 20] [cm]
%x4 = h4 \in [0.1; 20] [cm]

%parametros da planta (constantes do modelo)
global A_ k g gamma a v h0;

delta_x(1) = x(1) - h0(1);
delta_x(2) = x(2) - h0(2);
delta_x(3) = x(3) - h0(3);
delta_x(4) = x(4) - h0(4);


var_degrau = 1;

if var_degrau == 1
    delta_v1 = f_input(t,v(1));
    delta_v2 = f_input(t,v(2));
else
    delta_v1 = v(1);
    delta_v2 = v(2);
end



%function
dx1dt = (-a(1)*sqrt(2*g)*delta_x(1))/(A_*2*sqrt(h0(1))) + (a(3)*sqrt(2*g)*delta_x(3))/(A_*2*sqrt(h0(3)))...
        + (gamma(1)*k(1)*delta_v1)/A_;
dx2dt = (-a(2)*sqrt(2*g)*delta_x(2))/(A_*2*sqrt(h0(2))) + (a(4)*sqrt(2*g)*delta_x(4))/(A_*2*sqrt(h0(4)))...
        + gamma(2)*k(2)*delta_v2/A_;
dx3dt = (-a(3)*sqrt(2*g)*delta_x(3))/(A_*2*sqrt(h0(3))) + ((1 - gamma(2))*k(2)*delta_v2)/A_;
dx4dt = (-a(4)*sqrt(2*g)*delta_x(4))/(2*A_*sqrt(h0(4))) + ((1 - gamma(1))*k(1)*delta_v1)/A_;

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