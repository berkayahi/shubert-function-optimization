function ypr = gradfunc(x)
    % Define symbolic variables
    syms x1 x2

    % Initialize the Shubert function components
    F1 = 0;
    F2 = 0;

    % Calculate the Shubert function components
    for i = 1:5
        F1 = F1 + i * cos((i + 1) * x1 + i);
        F2 = F2 + i * cos((i + 1) * x2 + i);
    end

    % Combine the components to get the Shubert function
    F = F1 * F2;

    % Calculate the gradient of the Shubert function
    g = gradient(F, [x1, x2]);

    % Substitute the values of x1 and x2 with the input values and convert to double
    ypr = double(subs(g, [x1, x2], [x(1), x(2)]));
end