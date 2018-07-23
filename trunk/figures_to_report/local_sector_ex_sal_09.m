% local sector

x1 = -0.5:0.1:0.5;
a1 = 0;
a2 = 0.2474;
a1x1 = a1.*x1;
a2x1 = a2.*x1;
f_x1 = sin(x1.^2).*x1;
plot(x1, a1x1, 'r', x1, a2x1, 'm', x1, f_x1, 'b');
h = legend('0', '0.2484 x_1(t)','sin(x_1(t)^2) x_1(t)');
xlabel('x_1(t)');