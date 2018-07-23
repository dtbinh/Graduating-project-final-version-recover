function [h, A_fuzzy] = TS_fuzzy_model(x1, xi, lambda)

% nonlinearities number
N = 1;
z = sin(x1);

x1_range = xi(1):2*(xi(2)-xi(1))/100:xi(2);
syms y;
for i = 1:length(x1_range)
    s = vpasolve([y == z, x1 == x1_range(i)]);
    z_lim(i) = eval(s.y);
end
z_max = max(z_lim);
z_min = min(z_lim);

h(1) = (z - z_min)/(z_max - z_min);
h(2) = (z_max - z)/(z_max - z_min);

A_fuzzy(:, :, 1) = [-2 4;  -(1 + ((lambda*(1-z_max))/(2))) -2];
A_fuzzy(:, :, 2) = [-2 4;  -(1 + ((lambda*(1-z_min))/(2))) -2];

end

