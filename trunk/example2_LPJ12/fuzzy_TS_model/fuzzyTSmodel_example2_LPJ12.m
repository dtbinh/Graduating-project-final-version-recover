% Dinamic simulation of fuzzy Takagi-Sugeno modeling
clear; clear all; clc;

ts = [0 4];
options = odeset('RelTol',1e-4,'AbsTol',1e-9);

x_init = [-pi/3,-pi/3];
[t,x] = ode45('dxdt', ts, x_init, options);
[t_fuzzy,x_fuzzy] = ode23t('dxdt_fuzzy',ts, x_init, options);

hold off;
figure(1);
plot(t(:,1), x(:,1),'-', t(:,1), x(:,2), '-');
hold on;
plot(t_fuzzy(:,1), x_fuzzy(:,1),'.', t_fuzzy(:,1), x_fuzzy(:,2), '.');
xlabel('t');
ylabel('x');