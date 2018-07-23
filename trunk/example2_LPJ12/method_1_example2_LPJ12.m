clear all; clear; clc;

lambda = 20;
[xi, A_fuzzy, x1, x2, h, n] = problem_definition(lambda);

x_k = StateVariablesVertices(xi);
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

LMIs = LargestInvariantSetContainedInPolytope(LMIs, x_k, P);
[LMIs, crit] = EnlargementOfLargestInvariantSet(LMIs, P);

solvesdp(LMIs, crit, sdpsettings('solver', 'sedumi', 'verbose', 0));
[p,d]=checkset(LMIs);
pmin = min(checkset(LMIs));
display(pmin)
maxViolation = 1e-7; %minimization problem
if pmin  > -maxViolation %adopted precision for the minimum primal residual
    msgbox 'Stable (method 1)';
    output.P = double(P);
    P_n = {};
    P_n{1} = output.P;
    attraction_region_nonlinear_system(P_n{1}, A_*[x1;x2], xi, x1, x2, 'm');
    level_curve(P_n, 1/152, 'y');
else
    msgbox 'Not stable (method 1)';
end