function y = func(x)
    % The number of variable n = 2.
    x1 = x(1);
    x2 = x(2);
    s1 = 0;
    s2 = 0;
    for i = 1:5
        s1 = s1 + i * cos((i + 1) * x1 + i);
        s2 = s2 + i * cos((i + 1) * x2 + i);
    end
    y = s1 * s2;
end