% function sol = stability_conditions()

clear all; close all; clc;

global z_lim X;

% at this function the following LMIs will be solved, if there is
% an existing solution, to verify if the nonlinear system is stable.
% If there is not a solution (an exiting matrix P) so the system is
% not stable.
%
% A(alpha)' P(alpha) + P(alpha) A(alpha) + Q(gamma) J(theta) A(alpha) + ...
% ... X(alpha) 1' J(theta) A(alpha) + A(alpha)' J(theta)' 1 X(alpha)' < 0
%
% Where Q = sum(gamma^k *[P_1*x^k ... P_N*x^k])
% the vertices of X(alpha) and P(alpha) are the variables

% Breakeven point - initial conditions
X = dxdt_fsolve;
z_lim = bounds_membership(X);

n = 3; % number of state variables

A_alpha = vertices(z_lim);
% n_alpha is the number of vertices of A(alpha)
[~, ~, n_alpha] = size(A_alpha);

J_theta = jacobians_vertices();
% n_theta is the number of vertices of J(theta)
[~, ~, n_theta] = size(J_theta);

x_k = polyhedral_set();
% n_gamma is the number of vertices of Q(gamma)
[~, n_gamma] = size(x_k);

A = rolmipvar_Matrix(A_alpha, 'A', n_alpha, n_theta, n_gamma, [1 0 0]);

J = rolmipvar_Matrix(J_theta, 'J', n_alpha, n_theta, n_gamma, [0 1 0]);

theta = zeros(1, n_theta);
gamma = zeros(1, n_gamma);
for i = 1:n_alpha
    alpha = zeros(1, n_alpha);
    alpha(i) = 1;
    Pi_{i} = sdpvar(n,n,'symmetric');
    Pi{i} = {alpha, theta, gamma, Pi_{i}};
end
P = rolmipvar(Pi,'P', [n_alpha n_theta n_gamma], [1 0 0]);

alpha = zeros(1, n_alpha);
theta = zeros(1, n_theta);
for k = 1:n_gamma
    for i = 1:n_alpha
        if i == 1
            Qi_{k} = Pi_{i}*x_k(:, n_gamma);
        else
            Qi_{k} = horzcat(Qi_{k}, Pi_{i}*x_k(:, n_gamma));
        end
    end
    gamma = zeros(1, n_gamma);
    gamma(k) = 1;
    Qi{k} = {alpha, theta, gamma, Qi_{k}};
end
Q = rolmipvar(Qi,'Q', [n_alpha n_theta n_gamma], [0 0 1]);

theta = zeros(1, n_theta);
gamma = zeros(1, n_gamma);
for i = 1:n_alpha
    alpha = zeros(1, n_alpha);
    alpha(i) = 1;
    Xi{i} = {alpha, theta, gamma, sdpvar(n,n,'full')};
end
X_ = rolmipvar(Xi, 'X', [n_alpha n_theta n_gamma], [1 0 0]);

LMIs = [P > 0];
LMIs = [LMIs, [A'*P+P*A+Q*J*A+X_*ones(16,n)'*J*A+A'*J'*ones(16,n)*X_']<0];

apply_theorem = 1;

if apply_theorem == 1
    % LMI (9)
%     poly = polytope(x_k');
%     [H,K] = double(poly);
%     [q, ~] = size(H);
%     b_k = [];
%     for k = 1:q
%         b_k = [b_k (H(k, :)/K(k))'];
%     end
%     [r, c] = size(b_k');
% 
%     for k = 1:length(b_k)
%         LMIs = [LMIs, [1 b_k(:,k)' ;b_k(:,k) P] >=0];
%     end
%     [LMIs, crit] = EnlargementOfLargestInvariantSet(LMIs, P);
%     solvesdp(LMIs,crit,sdpsettings('solver','sedumi','verbose',0));
    solvesdp(LMIs,[],sdpsettings('solver','sedumi','verbose',0));
    [p,d]=checkset(LMIs);
    pmin = min(checkset(LMIs));
    display(pmin)
    maxViolation = 1e-7;
    if pmin > -maxViolation
        msgbox 'Stable  (method 4 + set enlargement')'
        output.P = double(P);
        P_n = verticesP(output.P, n_alpha, n_theta, n_gamma);
        %level_curve(P_n, 1, 'g');
%         level_curve_v(P_n, 1, xi, h, 'g');
    else
        msgbox 'Not stable (method 4 + set enlargement)'
    end
    else
    % theorem 2
    [LMIs, T, gamma] = theorem2(xi, P, LMIs, n_alpha, n_theta, n_gamma, n);
    omega1 = 1;
    omega2 = 1;
    %gamma = sdpvar(1,1,'symmetric');
    crit = omega1 * trace(T) + omega2 * gamma;
    crit = crit([0]);
    crit = crit{1};

    solvesdp(LMIs,crit,sdpsettings('solver','sedumi','verbose',0));
    [p,d]=checkset(LMIs);
    pmin = min(checkset(LMIs));
    display(pmin)
    maxViolation = 1e-7;
    if pmin > -maxViolation
        msgbox 'Stable  (method 4 (with theorem 2))'
        output.P = double(P);
        P_n = verticesP(output.P, n_alpha, n_theta, n_gamma);
%         level_curve(P_n, 1/double(gamma), 'r');
        level_curve_v(P_n, 1/double(gamma), xi, h, 'y');
    else
        msgbox 'Not stable (method 4 (with theorem 2))'
    end
end