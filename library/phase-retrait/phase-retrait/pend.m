g=1;
l=1;
k=0;
m=1;

x1=-pi/2:.1:3*pi/2;
x2=-3:.01:3;
[X1,X2]=meshgrid(x1,x2);

% phase portrait nonlinear e.d.o
% dX = F(X)
F1=X2;
F2=-(g/l)*sin(X1)-(k/m)*X2;

figure(1);clf
streamslice(X1,X2,F1,F2,2);
xlabel('x_1');ylabel('x_2');
line([0,pi],[0,0],'marker','o','linestyle','none','markerfacecolor','r')
set(gca,'ytick',[0]);
set(gca,'xtick',[0, pi]);
axis([-1.5 4.5 -3 3]);