function output = processo4tanques

close all; clc;

global A_ k g gamma a h0 v h_init h_min h_max nro_tanques;

fase_minima = 1;

[A_, k, g, gamma, a, h0, v, h_init, h_min, h_max, nro_tanques]...
                            = planta4tanquesParametros(fase_minima);
% Tempo de simulação
%ts = 0.0:0.01:60.0; % h
ts = [0 1e3];

%%  Simulação dinamica não-linear

% AbsTol and RelTol options to specify absolute and relative error tolerances
options = odeset('RelTol',1e-4,'AbsTol',1e-9); 
[t,x] = ode23t('dxdt',ts, h_init, options);
out.naolinear.t = t;
out.naolinear.h1 = x(:,1);
out.naolinear.h2 = x(:,2);
out.naolinear.h3 = x(:,3);
out.naolinear.h4 = x(:,4);

figure(1);
plot(t(:,1), x(:,1), t(:,1), x(:,2), t(:,1), x(:,3), t(:,1), x(:,4));
xlabel('t (s)'); ylabel ('Altura (cm)'); grid;
grid;
% 
% %% Simulação dinâmica em espaço de estados
% [t_ee,x_ee] = ode23t('dxdt_ee',ts, h_init, options);
% out.ee.t = t_ee;
% out.ee.h1 = x_ee(:,1); + h0(1);
% out.ee.h2 = x_ee(:,2); + h0(2);
% out.ee.h3 = x_ee(:,3); + h0(3);
% out.ee.h4 = x_ee(:,4); + h0(4);
% hold on;
% plot(t_ee(:,1), x_ee(:,1),'--', t_ee(:,1), x_ee(:,2), '--',...
%     t_ee(:,1), x_ee(:,3), '--', t_ee(:,1), x_ee(:,4), '--');

%% Simulação dinâmica modelagem fuzzy Takagi-Sugeno
[t_fuzzy,x_fuzzy] = ode23t('dxdt_fuzzy',ts, h_init, options);
out.fuzzy.t = t_fuzzy;
out.fuzzy.h1 = x_fuzzy(:,1);
out.fuzzy.h2 = x_fuzzy(:,2);
out.fuzzy.h3 = x_fuzzy(:,3);
out.fuzzy.h4 = x_fuzzy(:,4);
hold on;
plot(t_fuzzy(:,1), x_fuzzy(:,1),'.', t_fuzzy(:,1), x_fuzzy(:,2), '.',...
    t_fuzzy(:,1), x_fuzzy(:,3), '.', t_fuzzy(:,1), x_fuzzy(:,4), '.');
