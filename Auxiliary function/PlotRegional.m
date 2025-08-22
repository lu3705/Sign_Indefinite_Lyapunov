function PlotRegional(A, B, Ts, sol, x0, r, range, debug)

    range_x=range(1);
    range_y=range(2);
        
    K1=sol.K1;
    K2=sol.K2;
    P=sol.P;
    Phat=sol.Phat;
    H1=sol.H1;
    H2=sol.H2;
 
    [n,~] = size(B);
    
    [x1_grid,x2_grid]=meshgrid(linspace(-range_x,range_x,800),linspace(-range_y,range_y,800));
    
    Vgrid=zeros(size(x1_grid));
    Wvalues=zeros(size(x1_grid));
    H_values=zeros(size(x1_grid));
    Vhat_values=zeros(size(x1_grid));
    for i=1:size(x1_grid,1)
        for j=1:size(x1_grid,2)
            x = [x1_grid(i,j); x2_grid(i,j); zeros(n-2,1)];
    
            [~, ~, dz] = Saturation(x, K1, K2);
    
            Vgrid(i,j) = [x; dz]'*P*[x; dz];

            Vhat_values(i,j)=x'*Phat*x;
    
            hx=H1*x+H2*dz;
    
            if norm(hx,inf)<=1
                H_values(i,j) = norm(hx, inf); 
                Wvalues(i,j)=min(Vgrid(i,j), 1);
            else
                Wvalues(i,j)=1;
                H_values(i,j) = 1;
                % Vhat_values(i,j)=1;
            end
        end
    end

    traj=[];
    V=0;
    x=x0;
    
    N=100;
    for k=1:N/Ts
        
        [u, u_sat, dz] = Saturation(x, K1, K2);
        
        Vant=V;
        V = [x; dz]'*P*[x; dz];
        dV=V-Vant;
    
        data = [x; V; dV; u; u_sat];
        traj = [traj data];
    
        x = A*x + B*u_sat;

        if abs(x(1,:))>range_x
            disp("Trajectory exceeded range limit in X"); break;
        elseif abs(x(2,:))>range_y
            disp("Trajectory exceeded range limit in Y"); break;
        elseif all(abs(x)<1e-4)
            break;
        end

    end

    figure
        plot(traj(1,:)  ,traj(2,:),'.-'); hold on
        plot(traj(1,1)  ,traj(2,1),'xr','LineWidth',2); hold on
        contour(x1_grid, x2_grid, H_values, [1 1]); hold on;
        contour(x1_grid, x2_grid, Wvalues, [1 1], 'k'); hold on;
        contour(x1_grid, x2_grid, Vhat_values, [1 1], 'r');
        % label=sprintf('$S_h$');
        label=sprintf('$S_h : %s = %.2f$', inputname(6), r);
        legend("trajectory","start",label,"W(x)","$\hat{V}$","$V_0(x)$",'Interpreter','latex')

        xlabel("$x_1$",Interpreter='latex',FontSize=15); 
        ylabel("$x_2$",Interpreter='latex',FontSize=15);

        if inputname(6)=="r"
            label2=sprintf("Regional Synthesis - $ T_s = %.2f$", Ts);
        else 
            label2=sprintf("Regional Analysis - $ T_s = %.2f$", Ts);
        end
        title(label2,'Interpreter','latex');

    if debug
        t=0:Ts:Ts*(size(traj,2)-1);

        PlotPhasePlan(A, B, Ts, sol, range);

        figure
        subplot 211
            plot(t,traj(n+1,:),'.-');
            legend("V")
            titulo=sprintf('$T_s=%.2f$',Ts);
            title(titulo,'Interpreter','latex')
        subplot 212
            plot(t(2:end),traj(n+2, 2:end),'.-')
            legend("△V")

        figure
            surf(x1_grid,x2_grid,Vgrid); hold on
            contour3(x1_grid, x2_grid, Vgrid, [1 1], 'k', 'LineWidth', 1)
            colormap(parula)
            shading interp
            % zlim([-20,5]);
            legend("V", "V(x) = 1")
            xlabel("$x_1$",'Interpreter','latex', 'FontSize', 14);
            ylabel("$x_2$",'Interpreter','latex', 'FontSize', 14);
            zlabel("V", 'FontSize', 14)
    
        figure
        for i=1:n
            subplot(n,1,i)
            plot(t,traj(i,:),'.-');hold on;
            plot(t, zeros(1, length(t)), 'Color', [0.7 0.7 0.7], 'LineStyle','--');
            xlabel("Time")
            titulo=sprintf('x%d ',i);
            legend(titulo)
        end
    
        figure
            stairs(t,traj(n+3,:)); hold on;
            stairs(t,traj(n+4,:)); hold on;
            plot(t, zeros(1, length(t)), 'Color', [0.7 0.7 0.7], 'LineStyle','--');
            legend('$u$', '$u_{\mathrm{sat}}$', 'Interpreter', 'latex')
            xlabel("Time")
            title("Control signal")
    end
end