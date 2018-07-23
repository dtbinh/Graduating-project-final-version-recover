clear all; clear; clc;

lambda = 0;
pmin = 1e-20;
while pmin > 0
    [xi, ~, x1, x2, ~, n] = problem_definition(lambda);
    [h, A] = TS_fuzzy_model(x1, xi, lambda);
    
    poly_A = rolmipvar(A,'A',2,1);
    P = sdpvar(n,n,'symmetric');
    poly_P = rolmipvar(P,'P',2,0);

    LMIs = [];
    LMIs = [LMIs, poly_P >= 0];
    LMIs = [LMIs, poly_A'*poly_P + poly_P*poly_A <= 0];

    solvesdp(LMIs, [], sdpsettings('solver', 'sedumi', 'verbose', 0));
    [p,d]=checkset(LMIs);
    pmin = min(checkset(LMIs));
    display(pmin)
    if pmin > 0
        lambda_min = lambda;
        lambda = lambda - 0.1;
    end
end