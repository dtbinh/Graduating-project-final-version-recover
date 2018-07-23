function Ji = jacobians_vertices()

global z_lim X;
syms x1 x2 x3 y;

% non-linearities number
N = 4;

% membership degrees
Mij = membership_degrees([x1, x2, x3], N);

% nonlinearity \gamma2 = z2/x3 - k2*z4  (1st jacobian nonlinearity)
% \gamma4 = cos(x3)/x3 - sen(x3)/(x3)^2 (2nd jacobian nonlinearity)

% number of nonlinearities of jacobians
NJ = 6;

range_points = 100;
x3_range = linspace(-0.02, 0.1, range_points);
gamma2 = zeros(1,range_points);
gamma4 = zeros(1,range_points);
Mij_range = zeros(4,range_points);
[wf, V, Eref, ~, n, Ro, Pref, ~] = parameters;
k2 = wf*(V*(Eref-n*(X(1)-Pref)))/Ro;
for i = 1:range_points
    z2 = ((X(1)*wf+((wf*V^2)/Ro))/x3_range(i))+(cos(x3_range(i))/x3_range(i))*((-wf*V*(Eref+n*Pref-n*X(1)))/Ro);
    z4 = sin(x3_range(i))/x3_range(i);
    gamma2(i) = z2/x3_range(i) - k2*z4;
    gamma4(i) = cos(x3_range(i))/x3_range(i) - sin(x3_range(i))/(x3_range(i))^2;
    for c = 1:4
        s = vpasolve([y == Mij(1, c), x3 == x3_range(i)], [x3, y]);
        Mij_range(c,i) = double(s.y);
    end
end

% membership degrees limits (4 jacobians nonlinearities: M11 M21 M31 e M41)
% Mi2 = 1 - Mi1
Mij_lim = zeros(4,2);
for i = 1:4
    Mij_lim(i,1) = min(Mij_range(i,:));
    if(Mij_lim(i,1) <0)
        Mij_lim(i,1) = 0;
    end
    Mij_lim(i,2) = max(Mij_range(i,:));
    if(Mij_lim(i,2) > 1)
    	Mij_lim(i,2) = 1;
    end
end
gamma2_lim(1) = min(gamma2);
gamma2_lim(2) = max(gamma2);

gamma4_lim(1) = min(gamma4);
gamma4_lim(2) = max(gamma4);

delta_z1 = z_lim{1,2} - z_lim{1,1};
delta_z2 = z_lim{2,2} - z_lim{2,1};
delta_z3 = z_lim{3,2} - z_lim{3,1};
delta_z4 = z_lim{4,2} - z_lim{4,1};
% jacobians vertices of membership functions
VN = 2^NJ;
Ji = zeros(16,3,VN);

i = 0;
for m1 = 1:2
    for m2 = 1:2
        for m3 = 1:2
            for m4 = 1:2
                for g2 = 1:2
                    for g4 = 1:2
                        i = i+1;
    % alpha1 = M11M21M31M41
    % dalpha1_dx3 = M11'M21M31M41+M11M21'M31M41+M11M21M31'M41+M11M21M31M41'
	dalpha1_dx3 = -((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*Mij_lim(2,m2)*Mij_lim(3,m3)*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*(gamma2_lim(g2)/delta_z2)*Mij_lim(3,m3)*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*Mij_lim(2,m2)*(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3)*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*Mij_lim(2,m2)*Mij_lim(3,m3)*gamma4_lim(g4)/delta_z4;
	Ji(1, :, i) = [dalpha1_dx3 0 0];
    
    % alpha2 = M11M21M31M42
    % dalpha2_dx3 = M11'M21M31M42+M11M21'M31M42+M11M21M31'M42+M11M21M31M42'
    dalpha2_dx3 = -((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*Mij_lim(2,m2)*Mij_lim(3,m3)*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*(gamma2_lim(g2)/delta_z2)*Mij_lim(3,m3)*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*Mij_lim(2,m2)*(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3)*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*Mij_lim(2,m2)*Mij_lim(3,m3)*(-gamma4_lim(g4)/delta_z4);
	Ji(2, :, i) = [dalpha2_dx3 0 0];
    
    % alpha3 = M11M21M32M41
    % dalpha3_dx3 = M11'M21M32M41+M11M21'M32M41+M11M21M32'M41+M11M21M32M41'
    dalpha3_dx3 = -((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*Mij_lim(2,m2)*(1-Mij_lim(3,m3))*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*(gamma2_lim(g2)/delta_z2)*(1-Mij_lim(3,m3))*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*Mij_lim(2,m2)*(-(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3))*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*Mij_lim(2,m2)*(1-Mij_lim(3,m3))*(gamma4_lim(g4)/delta_z4);
	Ji(3, :, i) = [dalpha3_dx3 0 0];
    
    % alpha4 = M11M21M32M42
    % dalpha4_dx3 = M11'M21M32M42+M11M21'M32M42+M11M21M32'M42+M11M21M32M42'
    dalpha4_dx3 = -((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*Mij_lim(2,m2)*(1-Mij_lim(3,m3))*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*(gamma2_lim(g2)/delta_z2)*(1-Mij_lim(3,m3))*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*Mij_lim(2,m2)*(-(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3))*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*Mij_lim(2,m2)*(1-Mij_lim(3,m3))*(-gamma4_lim(g4)/delta_z4);
	Ji(4, :, i) = [dalpha4_dx3 0 0];
    
    % alpha5 = M11M22M31M41
    % dalpha5_dx3 = M11'M22M31M41+M11M22'M31M41+M11M22M31'M41+M11M22M31M41'
    dalpha5_dx3 = -((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*(1-Mij_lim(2,m2))*Mij_lim(3,m3)*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*(-gamma2_lim(g2)/delta_z2)*Mij_lim(3,m3)*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*(1-Mij_lim(2,m2))*(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3)*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*(1-Mij_lim(2,m2))*Mij_lim(3,m3)*gamma4_lim(g4)/delta_z4;
	Ji(5, :, i) = [dalpha5_dx3 0 0];
    
    % alpha6 = M11M22M31M42
    % dalpha6_dx3 = M11'M22M31M42+M11M22'M31M42+M11M22M31'M42+M11M22M31M42'
    dalpha6_dx3 = -((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*(1-Mij_lim(2,m2))*Mij_lim(3,m3)*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*(-gamma2_lim(g2)/delta_z2)*Mij_lim(3,m3)*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*(1-Mij_lim(2,m2))*(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3)*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*(1-Mij_lim(2,m2))*Mij_lim(3,m3)*(-gamma4_lim(g4)/delta_z4);
	Ji(6, :, i) = [dalpha6_dx3 0 0];
    
    % alpha7 = M11M22M32M41
    % dalpha7_dx3 = M11'M22M32M41+M11M22'M32M41+M11M22M32'M41+M11M22M32M41'
    dalpha7_dx3 = -((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*(1-Mij_lim(2,m2))*(1-Mij_lim(3,m3))*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*(-gamma2_lim(g2)/delta_z2)*(1-Mij_lim(3,m3))*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*(1-Mij_lim(2,m2))*(-(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3))*Mij_lim(4,m4)...
                  + Mij_lim(1,m1)*(1-Mij_lim(2,m2))*(1-Mij_lim(3,m3))*(gamma4_lim(g4)/delta_z4);
	Ji(7, :, i) = [dalpha7_dx3 0 0];
    
    % alpha8 = M11M22M32M42
    % dalpha8_dx3 = M11'M22M32M42+M11M22'M32M42+M11M22M32'M42+M11M22M32M42'
    dalpha8_dx3 = -((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*(1-Mij_lim(2,m2))*(1-Mij_lim(3,m3))*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*(-gamma2_lim(g2)/delta_z2)*(1-Mij_lim(3,m3))*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*(1-Mij_lim(2,m2))*(-(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3))*(1-Mij_lim(4,m4))...
                  + Mij_lim(1,m1)*(1-Mij_lim(2,m2))*(1-Mij_lim(3,m3))*(-gamma4_lim(g4)/delta_z4);
	Ji(8, :, i) = [dalpha8_dx3 0 0];
    
    % alpha9 = M12M21M31M41
    % dalpha9_dx3 = M12'M21M31M41+M12M21'M31M41+M12M21M31'M41+M12M21M31M41'
    dalpha9_dx3 = ((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*Mij_lim(2,m2)*Mij_lim(3,m3)*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*(gamma2_lim(g2)/delta_z2)*Mij_lim(3,m3)*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*Mij_lim(2,m2)*(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3)*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*Mij_lim(2,m2)*Mij_lim(3,m3)*gamma4_lim(g4)/delta_z4;
	Ji(9, :, i) = [dalpha9_dx3 0 0];
    
    % alpha10 = M12M21M31M42
    % dalpha10_dx3 = M12'M21M31M42+M12M21'M31M42+M12M21M31'M42+M12M21M31M42'
    dalpha10_dx3 = ((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*Mij_lim(2,m2)*Mij_lim(3,m3)*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*(gamma2_lim(g2)/delta_z2)*Mij_lim(3,m3)*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*Mij_lim(2,m2)*(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3)*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*Mij_lim(2,m2)*Mij_lim(3,m3)*(-gamma4_lim(g4)/delta_z4);
	Ji(10, :, i) = [dalpha10_dx3 0 0];
    
    % alpha11 = M12M21M32M41
    % dalpha11_dx3 = M12'M21M32M41+M12M21'M32M41+M12M21M32'M41+M12M21M32M41'
    dalpha11_dx3 = ((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*Mij_lim(2,m2)*(1-Mij_lim(3,m3))*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*(gamma2_lim(g2)/delta_z2)*(1-Mij_lim(3,m3))*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*Mij_lim(2,m2)*(-(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3))*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*Mij_lim(2,m2)*(1-Mij_lim(3,m3))*gamma4_lim(g4)/delta_z4;
	Ji(11, :, i) = [dalpha11_dx3 0 0];
    
    % alpha12 = M12M21M32M42
    % dalpha12_dx3 = M12'M21M32M42+M12M21'M32M42+M12M21M32'M42+M12M21M32M42'
    dalpha12_dx3 = ((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*Mij_lim(2,m2)*(1-Mij_lim(3,m3))*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*(gamma2_lim(g2)/delta_z2)*(1-Mij_lim(3,m3))*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*Mij_lim(2,m2)*(-(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3))*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*Mij_lim(2,m2)*(1-Mij_lim(3,m3))*(-gamma4_lim(g4)/delta_z4);
	Ji(12, :, i) = [dalpha12_dx3 0 0];
    
    % alpha13 = M12M22M31M41
    % dalpha13_dx3 = M12'M22M31M41+M12M22'M31M41+M12M22M31'M41+M12M22M31M41'
    dalpha13_dx3 = ((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*(1-Mij_lim(2,m2))*Mij_lim(3,m3)*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*(-gamma2_lim(g2)/delta_z2)*Mij_lim(3,m3)*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*(1-Mij_lim(2,m2))*(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3)*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*(1-Mij_lim(2,m2))*Mij_lim(3,m3)*gamma4_lim(g4)/delta_z4;
	Ji(13, :, i) = [dalpha13_dx3 0 0];
    
    % alpha14 = M12M22M31M42
    % dalpha14_dx3 = M12'M22M31M42+M12M22'M31M42+M12M22M31'M42+M12M22M31M42'
    dalpha14_dx3 = ((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*(1-Mij_lim(2,m2))*Mij_lim(3,m3)*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*(-gamma2_lim(g2)/delta_z2)*Mij_lim(3,m3)*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*(1-Mij_lim(2,m2))*(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3)*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*(1-Mij_lim(2,m2))*Mij_lim(3,m3)*(-gamma4_lim(g4)/delta_z4);
	Ji(14, :, i) = [dalpha14_dx3 0 0];
    
    % alpha15 = M12M22M32M41
    % dalpha15_dx3 = M12'M22M32M41+M12M22'M32M41+M12M22M32'M41+M12M22M32M41'
    dalpha15_dx3 = ((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*(1-Mij_lim(2,m2))*(1-Mij_lim(3,m3))*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*(-gamma2_lim(g2)/delta_z2)*(1-Mij_lim(3,m3))*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*(1-Mij_lim(2,m2))*(-(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3))*Mij_lim(4,m4)...
                  + (1-Mij_lim(1,m1))*(1-Mij_lim(2,m2))*(1-Mij_lim(3,m3))*gamma4_lim(g4)/delta_z4;
	Ji(15, :, i) = [dalpha15_dx3 0 0];
    
    % alpha16 = M12M22M32M42
    % dalpha16_dx3 = M12'M22M32M42+M12M22'M32M42+M12M22M32'M42+M12M22M32M42'
    dalpha16_dx3 = ((Mij_lim(3,m3)*delta_z3/+z_lim{3,1})/delta_z1)*(1-Mij_lim(2,m2))*(1-Mij_lim(3,m3))*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*(-gamma2_lim(g2)/delta_z2)*(1-Mij_lim(3,m3))*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*(1-Mij_lim(2,m2))*(-(Mij_lim(1,m1)*delta_z1/delta_z3+z_lim{1,1}/delta_z3))*(1-Mij_lim(4,m4))...
                  + (1-Mij_lim(1,m1))*(1-Mij_lim(2,m2))*(1-Mij_lim(3,m3))*(-gamma4_lim(g4)/delta_z4);
	Ji(16, :, i) = [dalpha16_dx3 0 0];
                    end
                end
            end
        end
    end
end

end