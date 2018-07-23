% [\dot{x1}  = [-2               4  * [x1
%  \dot{x2}]    -11+10sin(x1)   -2]    x2]

% ponto de equilibrio
X = dxdt_fsolve;

x1_range = -pi/2:pi/20:pi/2;
x2_range = -pi/2:pi/20:pi/2;
[X1,X2]=meshgrid(x1_range,x2_range);

F1=-2*X1+4*X2;
F2=-11*X1+10.*X1.*sin(X1)-2.*X2;

hold on
figure(1);clf
streamslice(X1,X2,F1,F2,1);
axis([min(x1_range),max(x1_range),min(x2_range),max(x2_range)])
xlabel('x_1');ylabel('x_2');
line(X(1), X(2),'marker','o','linestyle','none','markerfacecolor','r');