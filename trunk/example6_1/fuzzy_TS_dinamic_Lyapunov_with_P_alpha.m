function [LMIs, P, n_alpha, n_theta, n_gamma] = ...
                            fuzzy_TS_dinamic_Lyapunov_with_P_alpha(A_,...
                                                    h, x1, x2, x_k, xi)

    % LMI (8): A(alpha)' P(alpha) + P(alpha) A(alpha)...
    % + Q(gamma) J(theta) A(alpha) + X(alpha) 1' J(theta) A(alpha)...
    % + A(alpha)' J(theta)' 1 X(alpha)' < 0
    %
    % Where Q = sum(gamma^k *[P_1*x^k ... P_N*x^k])
    % the vertices of X(alpha) and P(alpha) are the variables

    [~, ~, n_alpha] = size(A_);

    dh_dx = jacobian([h(1); h(2); h(3); h(4)],[x1 x2]);

    points = 40;
    x1_range = linspace(xi(1), xi(2), points);
    x2_range = linspace(xi(1), xi(2), points);

    syms y;
    [rows, columns] = size(dh_dx);
    J_ = zeros(rows, columns, 2);
    for r = 1:rows
        for c = 1:columns
            for p1 = 1:points
                for p2 = 1:points
                    s = vpasolve([y == dh_dx(r, c), x1 == x1_range(p1),...
                                        x2 == x2_range(p2)],[x1, x2, y]);
                    J_values(p1, p2) = double(s.y);
                end
                J_range(r, c, 1, p1) = min(J_values(p1, :));
                J_range(r, c, 2, p1) = max(J_values(p1, :));
            end
            J_(r, c, 1) = min(J_range(r,c,1,:));
            J_(r, c, 2) = max(J_range(r,c,2,:));
        end
    end
    [~, ~, n_theta] = size(J_);
    
    [~, n_gamma] = size(x_k);

    A = rolmipvar_Matrix(A_, 'A', n_alpha, n_theta, n_gamma, [1 0 0]);

    J = rolmipvar_Matrix(J_, 'J', n_alpha, n_theta, n_gamma, [0 1 0]);

    [P, Q, X_] = vertices_to_be_determined(n_alpha, n_theta, n_gamma, x_k, 2);
    
    [one_rows, ~] = size(J);
    [~, one_columns] = size(X_);
    LMIs = [];
    LMIs = [LMIs, P>0];
    LMIs = [LMIs, A'*P + P*A + Q*J*A + X_* ones(one_rows,one_columns)'*J*A ...
                            + A'*J'*ones(one_rows,one_columns)*X_'< 0];
end

