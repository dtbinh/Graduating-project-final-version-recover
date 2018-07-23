function DXDT = dxdt_fuzzy(~, x)

A1 = [-2 4;-1 -2];
A2 = [-2 4; -21 -2];
h1 = (1 + sin(x(1)))/2;
if h1 > 1
    h1 = 1;
elseif h1 < 0
    h1 = 0;
end
h2 = 1 - h1;
DXDT = (h1*A1+h2*A2)*x;