
X = dxdt_fsolve;

x1_range = 0:100:25000;
x2_range = 70000:100:5000;
x3_range = -0.02:10:0.1;
[X1,X2,X3] = meshgrid(x1_range,x2_range, x3_range);

wf = 2*pi*60;
V=311;
Eref=311;
m=3.77/22000;
n=20/22000;
Ro=0.1;
Pref=22000;
Qref=0;

F1 = -wf*X1 + (wf*(V*(Eref-n*(X1-Pref))*cos(X3)-V^2))/Ro; 
F2 = -wf*X2 - (wf*V*(Eref-n*(X1-Pref))*sin(X3))/Ro; 
F3 = m*(X2-Qref); 

hold on
figure(1);clf
streamslice(X1,X2,X3,F1,F2,F3,[],[],[5]);
axis([min(x1_range),max(x1_range),min(x2_range),max(x2_range),min(x3_range),max(x3_range)])
xlabel('x_1');ylabel('x_2');zlabel('x_3');
line(X(1), X(2),'marker','o','linestyle','none','markerfacecolor','r');