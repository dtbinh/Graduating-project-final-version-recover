function X = dxdt_fsolve(X0)

if nargin == 0
    N = 0; %0,3,5
    X0 = [0 0];
end

X = fsolve(@dxdt,X0);

function DXDT = dxdt(x) 

% states
% x1    \in [-2; 2]
% x2    \in [-2; 2]

dx1dt = -2*x(1)+4*x(2); 
dx2dt = -11*x(1)+10*x(1)*sin(x(1))-2*x(2);

DXDT = [dx1dt dx2dt]';