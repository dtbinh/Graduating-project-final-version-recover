function [xi, A, x1, x2, h, n] = problem_definition()

xi = [-2; 2];

A(:, :, 1) = [-3 2;0 -0.9];
A(:, :, 2) = [-0.8 3;0 -0.9];
A(:, :, 3) = [-1.9 2;-0.5 0.1];
A(:, :, 4) = [0.1 3;-0.5 -2];

syms x1 x2;
h(1) = ((4 - power(x1,2))*(4 - power(x2,2)))/16;
h(2) = ((4 - power(x1,2))*power(x2,2))/16;
h(3) = (power(x1,2)*(4 - power(x2,2)))/16;
h(4) = (power(x1,2) * power(x2,2))/16;

% number of state variables
n = 2;

end

