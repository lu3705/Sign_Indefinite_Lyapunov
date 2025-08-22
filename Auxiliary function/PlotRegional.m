% PlotRegional: plots trajectories and stability regions for regional analysis/synthesis
% Inputs:
%   A, B  - system matrices
%   Ts    - sampling period
%   sol   - structure with controller and Lyapunov matrices (K1, K2, P, Phat, H1, H2)
%   x0    - initial condition
%   r     - optimization parameter (tau in analysis, r in synthesis)
%   range - plotting range [x_max, y_max]
%   debug - if true, generates additional plots
% Outputs:
%   none (plots only)

function PlotRegional(A, B, Ts, sol, x0, r, range, debug)

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
    
    V_grid       = zeros(size(x1_grid));
    W_values     = zeros(size(x1_grid));
    H_values     = zeros(size(x1_grid));
    Vhat_values  = zeros(size(x1_grid));
    for i = 1:size(x1_grid,1)
        for j = 1:size(x1_grid,2)
            x = [x1_grid(i,j); x2_grid(i,j); zeros(n-2,1)];
    
            [~, ~, dz] = Saturation(x, K1, K2);
    
            V_grid(i,j) = [x; dz]'*P*[x; dz];
            Vhat_values(i,j)  = x'*Phat*x;
    
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

    %% Simulate trajectory
    trajectory_data = [];
    V = 0;
    x = x0;
    
    num_iterations = 100;
    for k = 1:num_iterations/Ts
        [u, u_sat, dz] = Saturation(x, K1, K2);
        
        V_prev = V;
        V = [x; dz]'*P*[x; dz];
        dV = V - V_prev;
    
        step_data = [x; V; dV; u; u_sat];
        trajectory_data = [trajectory_data step_data];
    
        x = A*x + B*u_sat;

        if abs(x(1,:)) > range_x
            disp("Trajectory exceeded range limit in X"); break;
        elseif abs(x(2,:)) > range_y
            disp("Trajectory exceeded range limit in Y"); break;
        elseif all(abs(x) < 1e-4)
            break;
        end
    end

    %% Plot regions and trajectory
    figure
        plot(trajectory_data(1,:), trajectory_data(2,:),'.-'); hold on
        plot(trajectory_data(1,1), trajectory_data(2,1),'xr','LineWidth',2); hold on
        contour(x1_grid, x2_grid, H_values, [1 1]); hold on;
        contour(x1_grid, x2_grid, W_values, [1 1], 'k'); hold on;
        contour(x1_grid, x2_grid, Vhat_values, [1 1], 'r');
        label = sprintf('$S_h : %s = %.2f$', inputname(6), r);
        legend("trajectory","start",label,"W(x)","$\hat{V}$","$V_0(x)$",'Interpreter','latex')

        xlabel("$x_1$",'Interpreter','latex','FontSize',15); 
        ylabel("$x_2$",'Interpreter','latex','FontSize',15);

        if inputname(6)=="r"
            label2 = sprintf("Regional Synthesis - $ T_s = %.2f$", Ts);
        else 
            label2 = sprintf("Regional Analysis - $ T_s = %.2f$", Ts);
        end
        title(label2,'Interpreter','latex');

    %% Debug plots
    if debug
        t = 0:Ts:Ts*(size(trajectory_data,2)-1);

        PlotPhasePlanRegional(A, B, Ts, sol, range);

        figure
        subplot 211
            plot(t,trajectory_data(n+1,:),'.-');
            legend("V")
            title(sprintf('$T_s=%.2f$',Ts),'Interpreter','latex')
        subplot 212
            plot(t(2:end),trajectory_data(n+2,2:end),'.-')
            legend("△V")

        figure
            surf(x1_grid,x2_grid,V_grid); hold on
            contour3(x1_grid, x2_grid, V_grid, [1 1], 'k', 'LineWidth', 1)
            colormap(parula)
            shading interp
            legend("V", "V(x) = 1")
            xlabel("$x_1$",'Interpreter','latex','FontSize',14);
            ylabel("$x_2$",'Interpreter','latex','FontSize',14);
            zlabel("V",'FontSize',14)
    
        figure
        for i = 1:n
            subplot(n,1,i)
            plot(t,trajectory_data(i,:),'.-'); hold on;
            plot(t, zeros(1,length(t)), 'Color',[0.7 0.7 0.7],'LineStyle','--');
            xlabel("Time")
            legend(sprintf('x%d',i))
        end
    
        figure
            stairs(t,trajectory_data(n+3,:)); hold on;
            stairs(t,trajectory_data(n+4,:)); hold on;
            plot(t, zeros(1,length(t)), 'Color',[0.7 0.7 0.7],'LineStyle','--');
            legend('$u$', '$u_{\mathrm{sat}}$', 'Interpreter','latex')
            xlabel("Time")
            title("Control signal")
    end
end

