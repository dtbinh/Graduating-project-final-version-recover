%function out = main
clear all; close all; clc;

global X z_lim

% Regime permanente
% x1 = Pf    \in [0; 25000]
% x2 = Qf    \in [-70000; 5000]
% x3 = delta \in [-0.02; 0.1]

[wf, V, Eref, m, n, Ro, Pref, Qref] = parameters;

x1i = [0; 25000]';
x2i = [-70000; 5000]';
x3i = [-0.02; 0.1]';

N=10;
x1r=linspace(x1i(1),x1i(2),N);
x2r=linspace(x2i(1),x2i(2),N);
x3r=linspace(x3i(1),x3i(2),N);
            
% X0 = [10000 1000 0.01];  
X0 = [1000000 -100000 0.3];
%X0 = [x1r(i) x2r(j) x3r(k)];

Pf0 = X0(1);
Qf0 = X0(2);
delta0 = X0(3);

% Tempo de simulação
ts = [0 0.04];

%condicoes iniciais ponto de equilibrio
X = dxdt_fsolve;

z_lim = bounds_membership(X);

% Regime permanente
% x1 = Pf    \in [0; 25000]
% x2 = Qf    \in [-70000; 5000]
% x3 = delta \in [-0.02; 0.1]

%Pf0 = X(1);
%Qf0 = X(2);
%delta0 = X(3);

%Pf0 = 0.1;
%Qf0 = -70000;
%delta0 = -0.2;

options = odeset('RelTol',1e-4,'AbsTol',1e-9);
[t,x] = ode45('dxdt',ts,[Pf0 Qf0 delta0],options);

% ploting the phase portrait
xend = [];
xend = [xend; x(end,:)];
% deviation variable
% x(:,1) = x(:,1) - X(1);
% phase_portrait(x, xend-X);
% disp(xend - X);
phase_portrait(x, xend);
grid on;
disp(xend);

out.t = t;
out.naolinear.Pf = x(:,1);
out.naolinear.Qf = x(:,2);
out.naolinear.delta = x(:,3);

[t_fuzzy,x_fuzzy] = ode45('dxdt_fuzzy',ts, [Pf0-X(1) Qf0 delta0], options);
x_fuzzy(:,1) = x_fuzzy(:,1)+X(1);

out.t_fuzzy = t_fuzzy;
out.fuzzy.Pf = x_fuzzy(:,1)+X(1);
out.fuzzy.Qf = x_fuzzy(:,2);
out.fuzzy.delta = x_fuzzy(:,3);

% Visualização da simulação
figure;
plot(t,x(:,1),'-',t,x(:,2),'-',t,x(:,3),'-',...
     t_fuzzy,x_fuzzy(:,1),'.',t_fuzzy,x_fuzzy(:,2),'.',t_fuzzy,x_fuzzy(:,3),'.');
xlabel('t'); ylabel ('P_f, Q_f, \delta'); grid;

figure;
plot(t,x(:,1),'-',t_fuzzy,x_fuzzy(:,1),'.'); title ('P_f' );
xlabel('t'); ylabel ('P_f');
figure;
plot(t,x(:,2),'-',t_fuzzy,x_fuzzy(:,2),'.'); title ('Q_f' );
xlabel('t'); ylabel ('Q_f');
figure;
plot(t,x(:,3),'-',t_fuzzy,x_fuzzy(:,3),'.'); title ('\delta' );
xlabel('t'); ylabel ('\delta');