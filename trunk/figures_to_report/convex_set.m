clear all; clear; clc;

% x = 0:0.01:1.2;
% x = x(1:end-27);
% 
% % non-convex set
% y1 = 0:0.02:0.6;
% y2_aux = 0.3:0.02:0.6;
% y2 = zeros(1,length(y2_aux));
% for i = 1:length(y2_aux)
%     y2(i) = y2_aux(end - i + 1);
% end
% y3 = 0.3:0.02:0.8;
% y4_aux = 0:0.04:0.8;
% y4 = zeros(1,length(y4_aux));
% for i = 1:length(y4_aux)
%     y4(i) = y4_aux(end - i + 1);
% end
% y = [y1 y2 y3 y4];
% plot (x,y);

x = [0 0.3 0.4 0.7 1];
y = [0 0.6 0.4 0.8 0];
plot(x,y);

x = [0 0.1 0.1 0.3 1];
y = [0 0 0.6 0.8 0];
plot(x,y);

% % convex set
% y1 = 0:0.02:0.3;
% y2 = 0.3:0.01:0.5;
% y3 = 0.5:0.02:1;
% y4 = zeros(1,length(y3));
% for i = 1:length(y3)
%     y4(i) = y3(end - i + 1);
% end
% y5 = zeros(1,length(y2));
% for i = 1:length(y2)
%     y5(i) = y2(end - i + 1);
% end
% y6 = zeros(1,length(y1));
% for i = 1:length(y1)
%     y6(i) = y1(end - i + 1);
% end
% 
% y = [y1 y2 y3 y4 y5 y6];
% 
% length_x = length(y);
% 
% x = 0:1/length_x:1;
% x = x(1:end-1);
% 
% plot(x,y);