clear all; clear; clc;

lambda = 0;
pmin = 1e-20;
while pmin > 0
    lambda = lambda - 0.1;
    [xi, A, x1, x2, h, n] = problem_definition(lambda);

    x_k = StateVariablesVertices(xi);

    % method 4) Fuzzy dynamics TS + Lyapunov method with P(alpha): LMIs (8)-(9)

    % LMI (8)
    [LMIs, P, n_alpha, n_theta, n_gamma] = ...
                                    fuzzy_TS_dinamic_Lyapunov_with_P_alpha(...
                                                    A, h(1), h(2), x1, x2, x_k);
    solvesdp(LMIs,[],sdpsettings('solver','sedumi','verbose',0));
    [p,d]=checkset(LMIs);
    pmin = min(checkset(LMIs));
    display(pmin)
    display(lambda)
    if pmin > 0
        %         msgbox 'Stable  (method 4 - stability)'
        lambda_min = lambda;
    end
end