function x_k = StateVariablesVertices(x1)
    x2 = x1;
    l = 1;
    for i = 1:2
        for j = 1:2
            x_k(1, l) = x1(i);
            x_k(2, l) = x2(j);
            l = l + 1;
        end
    end

end

