function DXDT = dxdt_fuzzy(t, x)

global A_ k g gamma a v h_min h_max nro_tanques;

%% variacoes das nao linearidades z1 = z2 =z3 = z4 = z = sqrt(hi)/hi; i = 1, 2, 3, 4. 
h_range = h_min:0.1:h_max;
z_range = zeros(size(h_range));
for i = 1:length(h_range)
    z_range(i) = sqrt(h_range(i))/h_range(i);
end
z_max = max(z_range);
z_min = min(z_range);

%% graus de petinencia
Mij=zeros(2,nro_tanques);
z = zeros(nro_tanques);
for j = 1: nro_tanques
    z(j) = sqrt(x(j))/x(j);
    Mij(1,j) = (z_max - z(j))/(z_max - z_min);
    if (Mij(1,j) < 0)
        Mij(1,j) = 0;
        Mij(2,j) = 1;
    else
        if (Mij(1,j) > 1)
            Mij(1,j) = 1;
            Mij(2,j) = 0;
        else
            Mij(2,j) = (z(j) - z_min)/(z_max - z_min);
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
    
%% funcoes de pertinencia
mi = ones(1,2^nro_tanques);

mi(1)  = Mij(1,1)*Mij(1,2)*Mij(1,3)*Mij(1,4);
mi(2)  = Mij(1,1)*Mij(1,2)*Mij(1,3)*Mij(2,4);
mi(3)  = Mij(1,1)*Mij(1,2)*Mij(2,3)*Mij(1,4);
mi(4)  = Mij(1,1)*Mij(1,2)*Mij(2,3)*Mij(2,4);
mi(5)  = Mij(1,1)*Mij(2,2)*Mij(1,3)*Mij(1,4);
mi(6)  = Mij(1,1)*Mij(2,2)*Mij(1,3)*Mij(2,4);
mi(7)  = Mij(1,1)*Mij(2,2)*Mij(2,3)*Mij(1,4);
mi(8)  = Mij(1,1)*Mij(2,2)*Mij(2,3)*Mij(2,4);
mi(9)  = Mij(2,1)*Mij(1,2)*Mij(1,3)*Mij(1,4);
mi(10) = Mij(2,1)*Mij(1,2)*Mij(1,3)*Mij(2,4);
mi(11) = Mij(2,1)*Mij(1,2)*Mij(2,3)*Mij(1,4);
mi(12) = Mij(2,1)*Mij(1,2)*Mij(2,3)*Mij(2,4);
mi(13) = Mij(2,1)*Mij(2,2)*Mij(1,3)*Mij(1,4);
mi(14) = Mij(2,1)*Mij(2,2)*Mij(1,3)*Mij(2,4);
mi(15) = Mij(2,1)*Mij(2,2)*Mij(2,3)*Mij(1,4);
mi(16) = Mij(2,1)*Mij(2,2)*Mij(2,3)*Mij(2,4);

%% vertices
[Ai, Bi] = vertices_fuzzy(z_max, z_min, A_, k, g, gamma, a);

%% tempo de simulacao

var_degrau = 1;

if var_degrau == 1
    v1 = f_input(t,v(1));
    v2 = f_input(t,v(2));
else
    v1 = v(1);
    v2 = v(2);
end

%% modelo fuzzy

DXDT = zeros(1, 4)';

for i = 1:2^nro_tanques
    DXDT = DXDT + mi(i)*(Ai(:, :, i)*x + Bi*[v1 v2]');
end

%%
function fe = f_input(t,x)
% disturbio de entrada
% x = [x1 x2 t1 t2]
if t < 100 || t > 300
    fe = x;
else
    fe = x*1.2;
end
