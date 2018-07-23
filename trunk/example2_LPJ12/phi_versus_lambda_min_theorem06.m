clear all; clear; clc;

phi = 0.0;
points = 0;
while phi < 10
    phi = phi + 0.1;
    points = points + 1;
    phi_max = [phi phi];
    lambda = 0;
    pmin = 1e-20;
    
    phi_range(points) = phi;
     
    while pmin > 0
        lambda = lambda - 0.1;
        [xi, A, x1, x2, h, n] = problem_definition(lambda);

        poly_A = rolmipvar(A,'A',2,1);
        Pi = {};
        P_ = {};
        for i = 1:n
            alpha = zeros(1, n);
            alpha(i) = 1;
            Pi{i} = sdpvar(n, n, 'symmetric');
            P_{i} = {alpha, Pi{i}};
        end
        poly_P = rolmipvar(P_,'P', n, 1);

        LMIs = [];
        LMIs = Theorem06(LMIs, A, Pi, n, phi_max);

        solvesdp(LMIs,[],sdpsettings('solver','sedumi','verbose',0));
        [p,d]=checkset(LMIs);
        pmin = min(checkset(LMIs));
        display(pmin)
        display(lambda)
        display(phi_max)
        if pmin > 0
        %         msgbox 'Stable  (method 3 - stability)'
            lambda_min(points) = lambda;
        end
    end
end
plot(phi_range, lambda_min);