% PlotPhasePlan: phase-plane visualization of regional sets and vector field
% Inputs:
%   A, B   - system matrices
%   Ts     - sampling period
%   sol    - struct with fields K1, K2, P, Phat, H1, H2
%   range  - plotting range [x_max, y_max]
% Outputs:
%   none (plots only)

function  PlotPhasePlan(A, B, Ts, sol, range)

    %% Extract inputs
    range_x = range(1);
    range_y = range(2);
        
    K1   = sol.K1;
    K2   = sol.K2;
    P    = sol.P;
    Phat = sol.Phat;
    H1   = sol.H1;
    H2   = sol.H2;
 
    [n,~] = size(B);
    
    %% Build grids
    [x1_grid, x2_grid] = meshgrid(linspace(-range_x,range_x,800), linspace(-range_y,range_y,800));
    
    V_grid      = zeros(size(x1_grid));
    W_values    = zeros(size(x1_grid));
    H_values    = zeros(size(x1_grid));
    Vhat_values = zeros(size(x1_grid));
    for i = 1:size(x1_grid,1)
        for j = 1:size(x1_grid,2)
            x = [x1_grid(i,j); x2_grid(i,j); zeros(n-2,1)];
    
            [~, ~, dz] = Saturation(x, K1, K2);
    
            V_grid(i,j) = [x; dz]'*P*[x; dz];
            Vhat_values(i,j) = x'*Phat*x;
    
            hx = H1*x + H2*dz;
    
            if norm(hx,inf) <= 1
                H_values(i,j) = norm(hx, inf); 
                W_values(i,j) = min(V_grid(i,j), 1);
            else
                W_values(i,j) = 1;
                H_values(i,j) = 1;
            end
        end
    end

    figure
        S_h_mask   = double(H_values    < 1);
        W_mask     = double(W_values     < 1);
        Vhat_mask  = double(Vhat_values < 1);
        
        hold on
        contourf(x1_grid, x2_grid, S_h_mask,  [0.5 1], 'FaceColor','c',      'EdgeColor','none','FaceAlpha',0.1);
        contourf(x1_grid, x2_grid, W_mask,    [0.5 1], 'FaceColor','k',      'EdgeColor','none','FaceAlpha',0.11);
        contourf(x1_grid, x2_grid, Vhat_mask, [0.5 1], 'FaceColor','#FFB9BA','EdgeColor','none','FaceAlpha',0.4);

        contour(x1_grid, x2_grid, H_values,    [1 1]); hold on;
        contour(x1_grid, x2_grid, W_values,     [1 1], 'Color','#525252', 'LineWidth',0.1); hold on;
        contour(x1_grid, x2_grid, Vhat_values, [1 1], 'r');
        legend('$S_h$', "W(x)", "$\hat{V}$", 'Interpreter','latex')
        xlabel("$x_1$", 'Interpreter','latex', 'FontSize',15); 
        ylabel("$x_2$", 'Interpreter','latex', 'FontSize',15);
        title('Phase Plane');

    %% Plot trajectory
    [x01_grid,x02_grid] = meshgrid(linspace(-range_x,range_x,5), linspace(-range_y,range_y,14));
    
    num_iterations = 100;
    for i = 1:size(x01_grid,1)
        for j = 1:size(x01_grid,2)
            x = [x01_grid(i,j); x02_grid(i,j); zeros(n-2,1)];

            for k = 1:num_iterations/Ts
                [~, u_sat, ~] = Saturation(x, K1, K2);

                x_prev = x;
                x = A*x + B*u_sat;
                dx = x(1) - x_prev(1);
                dy = x(2) - x_prev(2);

                if abs(x(1,:)) > range_x || abs(x(2,:)) > range_y || all(abs(x) < 1e-4)
                    break;
                end
    
                h = quiver(x_prev(1), x_prev(2), dx, dy, 'Color',"#A02BA8", 'LineWidth',0.7, 'MaxHeadSize',0.5); hold on;
                set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
            end     
        end
    end

end

