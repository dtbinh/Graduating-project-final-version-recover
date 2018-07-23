function [DXDT_FUZZY ] = dxdt_fuzzy(t, x)


% states
% x1 = Pf    \in [0; 25000]
% x2 = Qf    \in [-70000; 5000]
% x3 = delta \in [-0.02; 0.1]

% nao linearidades
%
% z1 = (-wf*(1 + V*n*cox(x3)/Ro))
%       + (wf*V*(Eref + n*(-Pf0 + Pref))*(cos(x3)/x1)/Ro)*x1
%       - wf*(X(1)+V^2/Ro)
% z2 = sen(x3)
% z3 = sen(x3)/x3

% ponto de equilíbrio
global X z_lim

[wf, V, Eref, ~, n, Ro, Pref, ~] = parameters();

% graus de petinencia
N = 4;
Mij=zeros(2,N);
z  = zeros(1, N);
for j = 1:N
    if j == 1
        z(j) = cos(x(3));
    elseif j == 2
        z(j) = ((wf*(V*(Eref-n*(X(1)-Pref))*cos(x(3))-V^2))/Ro -wf*X(1))/x(3);
    elseif j == 3
        z(j) = sin(x(3));
    elseif j == 4
        z(j) = sin(x(3))/x(3);
    end
    Mij(1,j) = (z_lim{j,1} - z(j))/(z_lim{j,1} - z_lim{j,2});
    if (Mij(1,j) < 0)
        Mij(1,j) = 0;
        Mij(2,j) = 1;
    else
        if (Mij(1,j) > 1)
            Mij(1,j) = 1;
            Mij(2,j) = 0;
        else
            Mij(2,j) = (z(j) - z_lim{j,2})/(z_lim{j,1} - z_lim{j,2});
            if (Mij(2,j) < 0)
                Mij(2,j) = 0;
                Mij(1,j) = 1;
            else
                if (Mij(2,j) > 1)
                    Mij(2,j) = 1;
                    Mij(1,j) = 0;
                end
            end
        end
    end
end

% funcoes de pertinencia
mi = ones(1,2^N);

i = 0;
for j = 1:2
    for k = 1:2
        for l = 1:2
            for p = 1:2
                i = i + 1;
                %funcoes de pertinencia
                mi(i)  = Mij(j,1)*Mij(k,2)*Mij(l,3)*Mij(p,4);
            end
        end
    end
end

%vertices
Ai = vertices(z_lim);

% modelo fuzzy
dxdt = zeros(3, 1);

for i = 1:2^N
    dxdt = dxdt + mi(i)*(Ai(:, :, i)*x);
end

DXDT_FUZZY = [dxdt(1,:) dxdt(2,:) dxdt(3,:)]';


