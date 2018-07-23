x1 = -0.5:0.1:0.5;
M11 = sin(x1.^2)/sin(0.5^2);
M12 = 1 - M11;

plot(x1, M11, 'b', x1, M12, 'r');
h = legend('M_{11}', 'M_{12}');
xlabel('x_1(t)');