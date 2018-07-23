function [Ai, Bi] = vertices_fuzzy(z_max, z_min, A_, k, g, gamma, a)

%       [z(1)*((-a(1)*sqrt(2*g))/A_)    0       z(3)((a(3)*sqrt(2*g))/A_)        0
%   A =         0           z(2)((-a(2)*sqrt(2*g))/A_)      0       z(4)((a(4)*sqrt(2*g))/A_)
%               0                       0       z(3)((-a(3)*sqrt(2*g))/A_)      0
%               0                       0                   0       z(4)((-a(4)*sqrt(2*g))/A_)]

%       [gamma(1)*k(1)/A_       0
%   B =         0        gamma(2)*k(2)/A_
%               0     (1 - gamma(2))*k(2)/A_
%       (1 - gamma(1))*k(1)/A_  0           ]

Ai = zeros(4, 4, 4);

% z(1) = z_min; z(2) = z_min; z(3) = z_min; z(4) = z_min
Ai(:, :, 1) = [z_min*((-a(1)*sqrt(2*g))/A_) 0   z_min*((a(3)*sqrt(2*g))/A_)   0
        0   z_min*((-a(2)*sqrt(2*g))/A_)  0   z_min*((a(4)*sqrt(2*g))/A_)
        0   0   z_min*((-a(3)*sqrt(2*g))/A_)  0
        0   0   0   z_min*((-a(4)*sqrt(2*g))/A_)];
    
% z(1) = z_min; z(2) = z_min; z(3) = z_min; z(4) = z_max
Ai(:, :, 2) = [z_min*((-a(1)*sqrt(2*g))/A_) 0   z_min*((a(3)*sqrt(2*g))/A_)   0
        0   z_min*((-a(2)*sqrt(2*g))/A_)  0   z_max*((a(4)*sqrt(2*g))/A_)
        0   0   z_min*((-a(3)*sqrt(2*g))/A_)  0
        0   0   0   z_max*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_min; z(2) = z_min; z(3) = z_max; z(4) = z_min
Ai(:, :, 3) = [z_min*((-a(1)*sqrt(2*g))/A_) 0   z_max*((a(3)*sqrt(2*g))/A_)   0
        0   z_min*((-a(2)*sqrt(2*g))/A_)  0   z_min*((a(4)*sqrt(2*g))/A_)
        0   0   z_max*((-a(3)*sqrt(2*g))/A_)  0
        0   0   0   z_min*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_min; z(2) = z_min; z(3) = z_max; z(4) = z_max
Ai(:, :, 4) = [z_min*((-a(1)*sqrt(2*g))/A_) 0   z_max*((a(3)*sqrt(2*g))/A_)   0
        0   z_min*((-a(2)*sqrt(2*g))/A_)  0   z_max*((a(4)*sqrt(2*g))/A_)
        0   0   z_max*((-a(3)*sqrt(2*g))/A_)  0
        0   0   0   z_max*((-a(4)*sqrt(2*g))/A_)];
    
% z(1) = z_min; z(2) = z_max; z(3) = z_min; z(4) = z_min
Ai(:, :, 5) = [z_min*((-a(1)*sqrt(2*g))/A_) 0   z_min*((a(3)*sqrt(2*g))/A_)   0
        0   z_max*((-a(2)*sqrt(2*g))/A_)  0   z_min*((a(4)*sqrt(2*g))/A_)
        0   0   z_min*((-a(3)*sqrt(2*g))/A_)  0
        0   0   0   z_min*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_min; z(2) = z_max; z(3) = z_min; z(4) = z_max
Ai(:, :, 6) = [z_min*((-a(1)*sqrt(2*g))/A_) 0   z_min*((a(3)*sqrt(2*g))/A_)   0
        0   z_max*((-a(2)*sqrt(2*g))/A_)  0   z_max*((a(4)*sqrt(2*g))/A_)
        0   0   z_min*((-a(3)*sqrt(2*g))/A_)  0
        0   0   0   z_max*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_min; z(2) = z_max; z(3) = z_max; z(4) = z_min
Ai(:, :, 7) = [z_min*((-a(1)*sqrt(2*g))/A_) 0   z_max*((a(3)*sqrt(2*g))/A_)   0
        0   z_max*((-a(2)*sqrt(2*g))/A_)  0   z_min*((a(4)*sqrt(2*g))/A_)
        0   0   z_max*((-a(3)*sqrt(2*g))/A_)  0
        0   0   0   z_min*((-a(4)*sqrt(2*g))/A_)];
    
% z(1) = z_min; z(2) = z_max; z(3) = z_max; z(4) = z_max
Ai(:, :, 8) = [z_min*((-a(1)*sqrt(2*g))/A_) 0   z_max*((a(3)*sqrt(2*g))/A_)   0
        0   z_max*((-a(2)*sqrt(2*g))/A_)  0   z_max*((a(4)*sqrt(2*g))/A_)
        0   0   z_max*((-a(3)*sqrt(2*g))/A_)  0
        0   0   0   z_max*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_max; z(2) = z_min; z(3) = z_min; z(4) = z_min
Ai(:, :, 9) = [z_max*((-a(1)*sqrt(2*g))/A_) 0   z_min*((a(3)*sqrt(2*g))/A_)   0
        0   z_min*((-a(2)*sqrt(2*g))/A_)  0   z_min*((a(4)*sqrt(2*g))/A_)
        0   0   z_min*((-a(3)*sqrt(2*g))/A_)  0
        0   0   0   z_min*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_max; z(2) = z_min; z(3) = z_min; z(4) = z_max
Ai(:, :, 10) = [z_max*((-a(1)*sqrt(2*g))/A_) 0   z_min*((a(3)*sqrt(2*g))/A_)   0
         0   z_min*((-a(2)*sqrt(2*g))/A_)  0   z_max*((a(4)*sqrt(2*g))/A_)
         0   0   z_min*((-a(3)*sqrt(2*g))/A_)  0
         0   0   0   z_max*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_max; z(2) = z_min; z(3) = z_max; z(4) = z_min
Ai(:, :, 11) = [z_max*((-a(1)*sqrt(2*g))/A_) 0   z_max*((a(3)*sqrt(2*g))/A_)   0
         0   z_min*((-a(2)*sqrt(2*g))/A_)  0   z_min*((a(4)*sqrt(2*g))/A_)
         0   0   z_max*((-a(3)*sqrt(2*g))/A_)  0
         0   0   0   z_min*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_max; z(2) = z_min; z(3) = z_max; z(4) = z_max
Ai(:, :, 12) = [z_max*((-a(1)*sqrt(2*g))/A_) 0   z_max*((a(3)*sqrt(2*g))/A_)   0
         0   z_min*((-a(2)*sqrt(2*g))/A_)  0   z_max*((a(4)*sqrt(2*g))/A_)
         0   0   z_max*((-a(3)*sqrt(2*g))/A_)  0
         0   0   0   z_max*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_max; z(2) = z_max; z(3) = z_min; z(4) = z_min
Ai(:, :, 13) = [z_max*((-a(1)*sqrt(2*g))/A_) 0   z_min*((a(3)*sqrt(2*g))/A_)   0
         0   z_max*((-a(2)*sqrt(2*g))/A_)  0   z_min*((a(4)*sqrt(2*g))/A_)
         0   0   z_min*((-a(3)*sqrt(2*g))/A_)  0
         0   0   0   z_min*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_max; z(2) = z_max; z(3) = z_min; z(4) = z_max
Ai(:, :, 14) = [z_max*((-a(1)*sqrt(2*g))/A_) 0   z_min*((a(3)*sqrt(2*g))/A_)   0
         0   z_max*((-a(2)*sqrt(2*g))/A_)  0   z_max*((a(4)*sqrt(2*g))/A_)
         0   0   z_min*((-a(3)*sqrt(2*g))/A_)  0
         0   0   0   z_max*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_max; z(2) = z_max; z(3) = z_max; z(4) = z_min
Ai(:, :, 15) = [z_max*((-a(1)*sqrt(2*g))/A_) 0   z_max*((a(3)*sqrt(2*g))/A_)   0
         0   z_max*((-a(2)*sqrt(2*g))/A_)  0   z_min*((a(4)*sqrt(2*g))/A_)
         0   0   z_max*((-a(3)*sqrt(2*g))/A_)  0
         0   0   0   z_min*((-a(4)*sqrt(2*g))/A_)];

% z(1) = z_max; z(2) = z_max; z(3) = z_max; z(4) = z_max
Ai(:, :, 16) = [z_max*((-a(1)*sqrt(2*g))/A_) 0   z_max*((a(3)*sqrt(2*g))/A_)   0
         0   z_max*((-a(2)*sqrt(2*g))/A_)  0   z_max*((a(4)*sqrt(2*g))/A_)
         0   0   z_max*((-a(3)*sqrt(2*g))/A_)  0
         0   0   0   z_max*((-a(4)*sqrt(2*g))/A_)];

Bi = [gamma(1)*k(1)/A_	0
        0	gamma(2)*k(2)/A_
        0	(1 - gamma(2))*k(2)/A_
        (1 - gamma(1))*k(1)/A_	0];
     
end

