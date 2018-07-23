% PlotCircle_1.m:   Plot a circle of radius R centered at the origin
% Input: number of points

n = 100;

angle = 0:2*pi/n:2*pi;            % vector of angles at which points are drawn
R_delta = 0.3;                            % Unit radius
R_epsilon = 0.1;

x1 = R_delta*cos(angle);  y1 = R_delta*sin(angle);   % Coordinates of the circle
hold on;
x2 = R_epsilon*cos(angle);  y2 = R_epsilon*sin(angle);   % Coordinates of the circle
plot(x2,y2,'r',x1,y1,'b');                             % Plot the circles
hold on;
line(0,0,'marker','o','linestyle','none','markerfacecolor','r');
axis([-0.4 0.4 -0.4 0.4]);
h = legend('\epsilon', '\delta');
xlabel('x_1(t)');
ylabel('x_2(t)');
set(gca,'ytick',0);
set(gca,'xtick',0);
grid on;