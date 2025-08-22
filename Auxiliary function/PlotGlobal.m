% PlotGlobal: trajectories and stability visualization for global analysis
% Inputs:
%   Ac, Bc - continuous-time system matrices
%   Ts     - sampling period
%   sol    - structure with controller and Lyapunov matrix (K1, K2, P)
%   x0     - initial condition
%   range  - plotting range [x_max, y_max]
% Outputs:
%   none (plots only)

function PlotGlobal(Ac, Bc, Ts, sol, x0, range)

    %% Discretize
    [A, B] = Discretize(Ac, Bc, Ts, 0);
    dt = 0.01;

    range_x = range(1);
    range_y = range(2);

    K1 = sol.K1;
    K2 = sol.K2;
    P  = sol.P;

    [n, ~] = size(B);
    
    %% Build grid
    [x1_grid, x2_grid] = meshgrid(linspace(-range_x, range_x, 800), ...
                                  linspace(-range_y, range_y, 800));
    
    V_grid = zeros(size(x1_grid));
    for i = 1:size(x1_grid,1)
        for j = 1:size(x1_grid,2)
            x = [x1_grid(i,j); x2_grid(i,j); zeros(n-2,1)];
            [~, ~, dz] = Saturation(x, K1, K2);
            V_grid(i,j) = [x; dz]'*P*[x; dz];
        end
    end
    
    %% Discrete-time simulation
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
    
    %% Continuous-time simulation (Euler)
    x = x0;
    continuous_trajectory = [];
    for k = 1:num_iterations/dt
        continuous_trajectory = [continuous_trajectory x];
        [~, u_sat, ~] = Saturation(x, K1, K2);
        dx = Ac*x + Bc*u_sat;
        x = x + dt*dx;
    end

    %% Plot: discrete trajectory
    figure
        plot(trajectory_data(1,:), trajectory_data(2,:),'.-'); hold on
        plot(trajectory_data(1,1), trajectory_data(2,1),'xr','LineWidth',2); hold on
        legend("trajectory","start")
        xlabel("$x_1$",'Interpreter','latex'); 
        ylabel("$x_2$",'Interpreter','latex')
        axis equal

    %% Plot: V and ΔV
    t = 0:Ts:Ts*(size(trajectory_data,2)-1);

    figure
    subplot 211
        plot(t, trajectory_data(n+1,:), '.-');
        legend("V")
        title(sprintf('Discrete - $T_s=%.2f$', Ts), 'Interpreter','latex')
    subplot 212
        plot(t(2:end), trajectory_data(n+2, 2:end), '.-')
        legend("△V")
        xlabel("Time")
    
    %% Plot: V(x) surface
    figure
        surf(x1_grid, x2_grid, V_grid); hold on
        colormap(parula)
        shading interp
        legend("V")
        xlabel("$x_1$",'Interpreter','latex','FontSize',14);
        ylabel("$x_2$",'Interpreter','latex','FontSize',14);
        zlabel("V",'FontSize',14)
    
    %% Plot: states over time
    figure
    for i = 1:n
        subplot(n,1,i)
        plot(t, trajectory_data(i,:), '.-'); hold on;
        plot(t, zeros(1,length(t)), 'Color',[0.7 0.7 0.7], 'LineStyle','--');
        xlabel("Time")
        legend(sprintf('x%d', i))
    end

    %% Plot: control signals
    figure
        stairs(t, trajectory_data(n+3,:)); hold on;
        stairs(t, trajectory_data(n+4,:)); hold on;
        plot(t, zeros(1,length(t)), 'Color',[0.7 0.7 0.7], 'LineStyle','--');
        legend('$u$', '$u_{\mathrm{sat}}$', 'Interpreter','latex')
        xlabel("Time")
        title("Control signal")

end

