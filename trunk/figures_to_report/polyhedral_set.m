% line([-1; -1; 1; 1; -1],[-2; 2; 2; -2; -2], 'Color', 'b');
fill([-1; -1; 1; 1; -1],[-2; 2; 2; -2; -2],'b')
axis([-2 2 -3 3]);
line(-1, -2,'marker','o','linestyle','none','markerfacecolor','r');
line(-1, 2,'marker','o','linestyle','none','markerfacecolor','r');
line(1, -2,'marker','o','linestyle','none','markerfacecolor','r');
line(1, 2,'marker','o','linestyle','none','markerfacecolor','r');