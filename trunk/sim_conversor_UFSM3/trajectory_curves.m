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

% Tempo de simulação
ts = [0 1e-1];

%condicoes iniciais ponto de equilibrio
X = dxdt_fsolve;

for i = 1:N
    for j = 1:N
        for k = 1:N
            X0 = [x1r(i) x2r(j) x3r(k)];
%             X0 = [10000 1000 0.01];
            Pf0 = X0(1);
            Qf0 = X0(2);
            delta0 = X0(3);

            options = odeset('RelTol',1e-4,'AbsTol',1e-9);
            [t,x] = ode45('dxdt',ts,[Pf0 Qf0 delta0],options);

            % ploting the phase portrait
            xend = [];
            xend = [xend; x(end,:)];
            phase_portrait(x);
            hold on;
            disp(xend);
        end
    end
end
hold on;
grid;
line(X(1),X(2),X(3),'marker','o','linestyle','none','markerfacecolor','r')