clear all
close all
clc

% Define the range for x1 and x2 with higher resolution
X1 = -10:0.1:10;
X2 = -10:0.1:10;

% Create the meshgrid
[x1, x2] = meshgrid(X1, X2);

% Initialize the Shubert function components
F1 = zeros(size(x1));
F2 = zeros(size(x2));

% Calculate the Shubert function components
for i = 1:5
    F1 = F1 + i * cos((i + 1) * x1 + i);
    F2 = F2 + i * cos((i + 1) * x2 + i);
end

% Combine the components to get the Shubert function
F = F1 .* F2;

% Find the minimum value of the function
realFMin = min(min(F));
disp(['The minimum value of the Shubert function is ', num2str(realFMin)])

% Plot the contour of the Shubert function
figure
contourf(x1, x2, F)
title('Contour Plot of the Shubert Function')
xlabel('x1')
ylabel('x2')
hold on

% %% Newton-Raphson Algorithm with multiple starts
% colors = {'r', 'm', 'c'}; % Define colors for the 3 starts: red, magenta, cyan
% num_starts = 3; % Number of random starts
% 
% figure;
% contourf(x1, x2, F)
% title('Newton-Raphson Algorithm with Multiple Starts')
% xlabel('x1')
% ylabel('x2')
% hold on
% 
% for start = 1:num_starts
%     fprintf('Newton-Raphson Algorithm\n');
%     x = -10 + 20*rand(2,1);
%     epsilon = 10^(-4);
%     color = colors{start};
% 
%     tic
%     fprintf('Initial point %d: k=1, x1=%f, x2=%f, f(x)=%f\n', start, x(1), x(2), func(x))
%     plot(x(1), x(2), [color, '.'], 'MarkerSize', 15) % Plot initial point with the color
% 
%     x_next = x - inv(hessianfunc(x)) * gradfunc(x); %%Newton-Raphson Formulü
%     x_next = mod(x_next + 10, 20) - 10;
%     fprintf('Initial point %d: k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n', start, x_next(1), x_next(2), func(x_next), abs(realFMin - func(x_next)))
%     plot(x_next(1), x_next(2), [color, '*'], 'MarkerSize', 10)
%     k = 3;
% 
%     while(norm(gradfunc(x_next)) > epsilon)
%         x = x_next;
%         x_next = x - inv(hessianfunc(x)) * gradfunc(x);
%         x_next = mod(x_next + 10, 20) - 10;
%         fprintf('Initial point %d: k=%d, x1=%f, x2=%f, f(x)=%f, error=%f\n', start, k, x_next(1), x_next(2), func(x_next), abs(realFMin - func(x_next)))
%         plot(x_next(1), x_next(2), [color, '*'], 'MarkerSize', 10)
%         k = k + 1;
%     end
% 
%     % Mark the final point with a larger marker
%     plot(x_next(1), x_next(2), [color, 'o'], 'MarkerSize', 12, 'MarkerFaceColor', color)
%     toc
% end
% 
% % Add text annotations for the different initial points
% text(-12, 9, 'Init 1', 'Color', 'r', 'FontSize', 15, 'FontWeight', 'bold');
% text(-12, 8, 'Init 2', 'Color', 'm', 'FontSize', 15, 'FontWeight', 'bold');
% text(-12, 7, 'Init 3', 'Color', 'c', 'FontSize', 15, 'FontWeight', 'bold');
% 
% set(gca, 'fontsize', 35)
% set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
% 
% %% Hestenes-Stiefel Algorithm with multiple starts
% colors = {'r', 'm', 'c'}; % Define colors for the 3 starts: red, magenta, cyan
% num_starts = 3; % Number of random starts
% 
% figure;
% contourf(x1, x2, F)
% title('Hestenes-Stiefel Algorithm with Multiple Starts')
% xlabel('x1')
% ylabel('x2')
% hold on
% 
% for start = 1:num_starts
%     fprintf('Hestenes-Stiefel Algorithm\n');
%     x = -10 + 20*rand(2,1);
%     epsilon = 10^(-4);
%     color = colors{start};
% 
%     tic
%     fprintf('Initial point %d: k=1, x1=%f, x2=%f, f(x)=%f\n', start, x(1), x(2), func(x))
%     plot(x(1), x(2), [color, '.'], 'MarkerSize', 15) % Plot initial point with the color
%     g = gradfunc(x);
%     d = -g;
% 
%     % alpha argmin procedure
%     alpha = 0:0.01:1; % alpha is between 0 and 1 interval.
%     func_alpha = zeros(length(alpha), 1);
%     for i=1:length(alpha)
%         func_alpha(i)=func(x+alpha(i) * d);
%     end
%     [val, ind] = min(func_alpha); %%minimumu kaçıncı indiste bulunuyor?
%     alpha=alpha(ind);
%     % end of alpha argmin procedure.
%     x_next = x + alpha * d;
%     x_next = mod(x_next + 10, 20) - 10;
%     g_next = gradfunc(x_next);
%     beta = (g_next'*(g_next-g)) / (d'*(g_next-g));
%     d_next = -g_next+beta*d;
% 
%     fprintf('Initial point %d: k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n', start, x_next(1), x_next(2), func(x_next), abs(realFMin-func(x_next)))
%     plot(x_next(1), x_next(2), [color, '*'], 'MarkerSize', 10)
%     k = 3;
% 
%     while(norm(gradfunc(x_next)) > epsilon && alpha ~= 0)
%         x = x_next;
%         g = g_next;
%         d = d_next;
% 
%         % alpha argmin procedure
%         alpha = 0:0.01:1; % alpha is between 0 and 1 interval.
%         func_alpha = zeros(length(alpha), 1);
%         for i=1:length(alpha)
%             func_alpha(i)=func(x + alpha(i) * d);
%         end
%         [val, ind] = min(func_alpha); %%minimumu kaçıncı indiste bulunuyor?
%         alpha = alpha(ind);
%         % end of alpha argmin procedure
% 
%         x_next = x + alpha * d;
%         x_next = mod(x_next + 10, 20) - 10;
%         g_next = gradfunc(x_next);
%         beta = (g_next'*(g_next-g)) / (d'*(g_next-g));
%         d_next = -g_next + beta * d;
% 
%         fprintf('Initial point %d: k=%d, x1=%f, x2=%f, f(x)=%f, error=%f\n', start, k, x_next(1), x_next(2), func(x_next), abs(realFMin-func(x_next)))
%         plot(x_next(1), x_next(2), [color, '*'], 'MarkerSize', 10)
%         k = k + 1;
%     end
% 
%     % Mark the final point with a larger marker
%     plot(x_next(1), x_next(2), [color, 'o'], 'MarkerSize', 12, 'MarkerFaceColor', color)
%     toc
% end
% 
% % Add text annotations for the different initial points
% text(-12, 9, 'Init 1', 'Color', 'r', 'FontSize', 15, 'FontWeight', 'bold');
% text(-12, 8, 'Init 2', 'Color', 'm', 'FontSize', 15, 'FontWeight', 'bold');
% text(-12, 7, 'Init 3', 'Color', 'c', 'FontSize', 15, 'FontWeight', 'bold');
% 
% title('Hestenes-Stiefel Algorithm with Multiple Starts')
% set(gca, 'fontsize', 35)
% set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
% 
% %% Polak-Ribiere Algorithm with multiple starts
% colors = {'r', 'm', 'c'}; % Define colors for the 3 starts: red, magenta, cyan
% num_starts = 3; % Number of random starts
% 
% figure;
% contourf(x1, x2, F)
% title('Polak-Ribiere Algorithm with Multiple Starts')
% xlabel('x1')
% ylabel('x2')
% hold on
% 
% for start = 1:num_starts
%     fprintf('Polak-Ribiere Algorithm\n');
%     x = -10 + 20*rand(2,1);
%     epsilon = 10^(-4);
%     color = colors{start};
% 
%     tic
%     fprintf('Initial point %d: k=1, x1=%f, x2=%f, f(x)=%f\n', start, x(1), x(2), func(x))
%     plot(x(1), x(2), [color, '.'], 'MarkerSize', 15) % Plot initial point with the color
%     g = gradfunc(x);
%     d = -g;
% 
%     % alpha argmin procedure
%     alpha = 0:0.01:1; % alpha is between 0 and 1 interval.
%     func_alpha = zeros(length(alpha), 1);
%     for i = 1:length(alpha)
%         func_alpha(i) = func(x + alpha(i) * d);
%     end
%     [val, ind] = min(func_alpha); % minimumu kaçıncı indiste bulunuyor?
%     alpha = alpha(ind);
%     % end of alpha argmin procedure
%     x_next = x + alpha * d;
%     x_next = mod(x_next + 10, 20) - 10;
%     g_next = gradfunc(x_next);
%     beta = (g_next' * (g_next - g)) / (g' * g); % Polak-Ribiere formula
%     d_next = -g_next + beta * d;
% 
%     fprintf('Initial point %d: k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n', start, x_next(1), x_next(2), func(x_next), abs(realFMin - func(x_next)))
%     plot(x_next(1), x_next(2), [color, '*'], 'MarkerSize', 10)
%     k = 3;
% 
%     while norm(gradfunc(x_next)) > epsilon && alpha ~= 0
%         x = x_next;
%         g = g_next;
%         d = d_next;
% 
%         % alpha argmin procedure
%         alpha = 0:0.01:1; % alpha is between 0 and 1 interval.
%         func_alpha = zeros(length(alpha), 1);
%         for i = 1:length(alpha)
%             func_alpha(i) = func(x + alpha(i) * d);
%         end
%         [val, ind] = min(func_alpha); % minimumu kaçıncı indiste bulunuyor?
%         alpha = alpha(ind);
%         % end of alpha argmin procedure
% 
%         x_next = x + alpha * d;
%         x_next = mod(x_next + 10, 20) - 10;
%         g_next = gradfunc(x_next);
%         beta = (g_next' * (g_next - g)) / (g' * g); % Polak-Ribiere formula
%         d_next = -g_next + beta * d;
% 
%         fprintf('Initial point %d: k=%d, x1=%f, x2=%f, f(x)=%f, error=%f\n', start, k, x_next(1), x_next(2), func(x_next), abs(realFMin - func(x_next)))
%         plot(x_next(1), x_next(2), [color, '*'], 'MarkerSize', 10)
%         k = k + 1;
%     end
% 
%     % Mark the final point with a larger marker
%     plot(x_next(1), x_next(2), [color, 'o'], 'MarkerSize', 12, 'MarkerFaceColor', color)
%     toc
% end
% 
% % Add text annotations for the different initial points
% text(-12, 9, 'Init 1', 'Color', 'r', 'FontSize', 15, 'FontWeight', 'bold');
% text(-12, 8, 'Init 2', 'Color', 'm', 'FontSize', 15, 'FontWeight', 'bold');
% text(-12, 7, 'Init 3', 'Color', 'c', 'FontSize', 15, 'FontWeight', 'bold');
% 
% set(gca, 'fontsize', 35)
% set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
% 
% %% Fletcher-Reeves Algorithm with multiple starts
% colors = {'r', 'm', 'c'}; % Define colors for the 3 starts: red, magenta, cyan
% num_starts = 3; % Number of random starts
% 
% figure;
% contourf(x1, x2, F)
% title('Fletcher-Reeves Algorithm with Multiple Starts')
% xlabel('x1')
% ylabel('x2')
% hold on
% 
% for start = 1:num_starts
%     fprintf('Fletcher-Reeves Algorithm\n');
%     x = -10 + 20*rand(2,1);
%     epsilon = 10^(-4);
%     color = colors{start};
% 
%     tic
%     fprintf('Initial point %d: k=1, x1=%f, x2=%f, f(x)=%f\n', start, x(1), x(2), func(x))
%     plot(x(1), x(2), [color, '.'], 'MarkerSize', 15) % Plot initial point with the color
%     g = gradfunc(x);
%     d = -g;
% 
%     % alpha argmin procedure
%     alpha = 0:0.01:1; % alpha is between 0 and 1 interval.
%     func_alpha = zeros(length(alpha), 1);
%     for i = 1:length(alpha)
%         func_alpha(i) = func(x + alpha(i) * d);
%     end
%     [val, ind] = min(func_alpha); % minimumu kaçıncı indiste bulunuyor?
%     alpha = alpha(ind);
%     % end of alpha argmin procedure
%     x_next = x + alpha * d;
%     x_next = mod(x_next + 10, 20) - 10;
%     g_next = gradfunc(x_next);
%     beta = (g_next' * g_next) / (g' * g); % Fletcher-Reeves formula
%     d_next = -g_next + beta * d;
% 
%     fprintf('Initial point %d: k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n', start, x_next(1), x_next(2), func(x_next), abs(realFMin - func(x_next)))
%     plot(x_next(1), x_next(2), [color, '*'], 'MarkerSize', 10)
%     k = 3;
% 
%     while norm(gradfunc(x_next)) > epsilon && alpha ~= 0
%         x = x_next;
%         g = g_next;
%         d = d_next;
% 
%         % alpha argmin procedure
%         alpha = 0:0.01:1; % alpha is between 0 and 1 interval.
%         func_alpha = zeros(length(alpha), 1);
%         for i = 1:length(alpha)
%             func_alpha(i) = func(x + alpha(i) * d);
%         end
%         [val, ind] = min(func_alpha); % minimumu kaçıncı indiste bulunuyor?
%         alpha = alpha(ind);
%         % end of alpha argmin procedure
% 
%         x_next = x + alpha * d;
%         x_next = mod(x_next + 10, 20) - 10;
%         g_next = gradfunc(x_next);
%         beta = (g_next' * g_next) / (g' * g); % Fletcher-Reeves formula
%         d_next = -g_next + beta * d;
% 
%         fprintf('Initial point %d: k=%d, x1=%f, x2=%f, f(x)=%f, error=%f\n', start, k, x_next(1), x_next(2), func(x_next), abs(realFMin - func(x_next)))
%         plot(x_next(1), x_next(2), [color, '*'], 'MarkerSize', 10)
%         k = k + 1;
%     end
% 
%     % Mark the final point with a larger marker
%     plot(x_next(1), x_next(2), [color, 'o'], 'MarkerSize', 12, 'MarkerFaceColor', color)
%     toc
% end
% 
% % Add text annotations for the different initial points
% text(-12, 9, 'Init 1', 'Color', 'r', 'FontSize', 15, 'FontWeight', 'bold');
% text(-12, 8, 'Init 2', 'Color', 'm', 'FontSize', 15, 'FontWeight', 'bold');
% text(-12, 7, 'Init 3', 'Color', 'c', 'FontSize', 15, 'FontWeight', 'bold');
% 
% set(gca, 'fontsize', 35)
% set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Dai-Yuan Algorithm
colors = {'r', 'm', 'c'}; % Define colors for the 3 starts: red, magenta, cyan
num_starts = 3; % Number of random starts

figure;
contourf(x1, x2, F)
title('Dai-Yuan Algorithm with Multiple Starts')
xlabel('x1')
ylabel('x2')
hold on

for start = 1:num_starts
    fprintf('Dai-Yuan Algorithm\n');
    x = -10 + 20*rand(2,1);
    epsilon = 10^(-4);
    color = colors{start};

    tic
    fprintf('Initial point %d: k=1, x1=%f, x2=%f, f(x)=%f\n', start, x(1), x(2), func(x))
    plot(x(1), x(2), [color, '.'], 'MarkerSize', 15) % Plot initial point with the color
    g = gradfunc(x);
    d = -g;

    % alpha argmin procedure
    alpha = 0:0.01:1; % alpha is between 0 and 1 interval.
    func_alpha = zeros(length(alpha), 1);
    for i = 1:length(alpha)
        func_alpha(i) = func(x + alpha(i) * d);
    end
    [val, ind] = min(func_alpha); % minimumu kaçıncı indiste bulunuyor?
    alpha = alpha(ind);
    % end of alpha argmin procedure
    x_next = x + alpha * d;
    x_next = mod(x_next + 10, 20) - 10;
    g_next = gradfunc(x_next);
    % Dai-Yuan formula for beta
    beta = (norm(g_next)^2) / (d' * (g_next - g));
    d_next = -g_next + beta * d;

    fprintf('Initial point %d: k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n', start, x_next(1), x_next(2), func(x_next), abs(realFMin - func(x_next)))
    plot(x_next(1), x_next(2), [color, '*'], 'MarkerSize', 10)
    k = 3;

    while norm(gradfunc(x_next)) > epsilon && alpha ~= 0
        x = x_next;
        g = g_next;
        d = d_next;

        % alpha argmin procedure
        alpha = 0:0.01:1; % alpha is between 0 and 1 interval.
        func_alpha = zeros(length(alpha), 1);
        for i = 1:length(alpha)
            func_alpha(i) = func(x + alpha(i) * d);
        end
        [val, ind] = min(func_alpha); % minimumu kaçıncı indiste bulunuyor?
        alpha = alpha(ind);
        % end of alpha argmin procedure

        x_next = x + alpha * d;
        x_next = mod(x_next + 10, 20) - 10;
        g_next = gradfunc(x_next);
        % Dai-Yuan formula for beta
        beta = (norm(g_next)^2) / (d' * (g_next - g));
        d_next = -g_next + beta * d;

        fprintf('Initial point %d: k=%d, x1=%f, x2=%f, f(x)=%f, error=%f\n', start, k, x_next(1), x_next(2), func(x_next), abs(realFMin - func(x_next)))
        plot(x_next(1), x_next(2), [color, '*'], 'MarkerSize', 10)
        k = k + 1;
    end

    % Mark the final point with a larger marker
    plot(x_next(1), x_next(2), [color, 'o'], 'MarkerSize', 12, 'MarkerFaceColor', color)
    toc
end

% Add text annotations for the different initial points
text(-12, 9, 'Init 1', 'Color', 'r', 'FontSize', 15, 'FontWeight', 'bold');
text(-12, 8, 'Init 2', 'Color', 'm', 'FontSize', 15, 'FontWeight', 'bold');
text(-12, 7, 'Init 3', 'Color', 'c', 'FontSize', 15, 'FontWeight', 'bold');

set(gca, 'fontsize', 35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

