clear all; clear; clc;

lambda = 4e+12;
pmin = -1e-20;
while pmin < 0
    [xi, A_fuzzy, x1, x2, h, n] = problem_definition(lambda);

    A_ = h(1)*A_fuzzy(:, :, 1) + h(2)*A_fuzzy(:,:,2);
    A_lin = taylor(A_, x1);
    syms y;
    s = vpasolve([y == A_lin(2, 1), x1 == 0]);
    A = [A_lin(1, 1) A_lin(1, 2); double(s.y) A_lin(2,2)];
    A = eval(A);

    P = sdpvar(n,n,'symmetric');

    LMIs = [];
    LMIs = [LMIs, P >= 0];
    LMIs = [LMIs, A'*P + P*A <= 0];

    solvesdp(LMIs, [], sdpsettings('solver', 'sedumi', 'verbose', 0));
    [p,d]=checkset(LMIs);
    pmin = min(checkset(LMIs));
    display(pmin)
    if pmin > 0
        lambda_max = lambda;
        lambda = lambda - 1;
    end
end