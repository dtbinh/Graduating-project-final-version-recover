function attraction_region_nonlinear_system(P, fx, xi, x1, x2, color)

dVx = 2*[x1 x2]*double(P)*fx;

xx1 = xi(1):0.1:xi(2);
xx2 = xi(1):0.1:xi(2);
[X1,X2] = meshgrid(xx1,xx2);

for i=1:length(xx1)
    for j=1:length(xx2)
        x1 = X1(i,j);
        x2 = X2(i,j);
        dV(i,j) = eval(subs(dVx));
    end
end
contour(X1, X2, dV,[0 0], 'LineColor', color);

end

