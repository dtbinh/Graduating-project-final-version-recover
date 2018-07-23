function [A_, k, g, gamma, a, h0, v, h_init, h_min, h_max, nro_tanques]...
                                = planta4tanquesParametros(fase_minima)
                            
% numero de tanques
nro_tanques = 4;
                            
% entradas
v(1) = 1;
v(2) = 1;

% valores iniciais
h_init = ones(1,4);

% areas das secao transversal dos tanques em [cm^2]
% A1 =~ A2 =~ A3 =~ A4 = A_
A_ = 47.6;

% ganhos das bombas [cm^2/V*s]
k(1) = 8.63;
k(2) = 12.72;

% aceleracao da gravidade em Brasília [cm/s^2]
g = 978;

% para fase minima
    % taxa de liquido desviado para o tanque 1
    gamma(1) = 0.7363;
    % taxa de liquido desviado para o tanque 2
    gamma(2) = 0.7577;


    % altura estacionaria do nivel de liquido nos tanques [cm]
    h0(1) = 8.9;
    h0(2) = 9.97;
    h0(3) = 8.65;
    h0(4) = 9.67;
    
    %h0(1) = 0;
    %h0(2) = 0;
    %h0(3) = 0;
    %h0(4) = 0;
    
% para fase nao-minima
if ~fase_minima
    % taxa de liquido desviado para o tanque 1
    gamma(1) = 0.1446;
    % taxa de liquido desviado para o tanque 2
    gamma(2) = 0.1330;

    % altura estacionaria do nivel de liquido nos tanques [cm]
    h0(1) = 8.5761;
    h0(2) = 9.0253;
    h0(3) = 8.1070;
    h0(4) = 8.8679;
end

% areas dos tubos que flui pra fora do tanque i [cm^2]
a(1) = 0.071;
a(2) = 0.057;
a(3) = 0.071;
a(4) = 0.057;

% valores minimo e maximo que as alturas dos tanques podem atingir [cm]
h_min = 0.1; % nao se escolhe h_min = 0, porque este é um ponto de descontinuidade
h_max = 23;

end