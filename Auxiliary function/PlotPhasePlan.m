function  PlotPhasePlan(A,B,Ts,sol, range)

    N=100;

    range_x=range(1);
    range_y=range(2);
        
    K1=sol.K1;
    K2=sol.K2;
    P0=sol.P0;
    P=sol.P;
    Phat=sol.Phat;
    H1=sol.H1;
    H2=sol.H2;
 
    [n,~] = size(B);

    [x1_grid,x2_grid]=meshgrid(linspace(-range_x,range_x,800),linspace(-range_y,range_y,800));
    
    V0=zeros(size(x1_grid));
    Vgrid=zeros(size(x1_grid));
    Wvalues=zeros(size(x1_grid));
    H_values=zeros(size(x1_grid));
    Vhat_values=zeros(size(x1_grid));
    for i=1:size(x1_grid,1)
        for j=1:size(x1_grid,2)
            x = [x1_grid(i,j); x2_grid(i,j); zeros(n-2,1)];
    
            [~, ~, dz] = Saturation(x, K1, K2);
    
            Vgrid(i,j) = [x; dz]'*P*[x; dz];
            V0(i,j) = x'*P0*x;

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

    figure

        % Máscaras lógicas para as regiões internas
        S_h_mask      = double(H_values    < 1);
        W_mask        = double(Wvalues     < 1);
        Vhat_mask     = double(Vhat_values < 1);
        
        % Plota o preenchimento com opacidade (usar transparência se desejar ver sobreposições)
        hold on
        contourf(x1_grid, x2_grid, S_h_mask, [0.5 1], 'FaceColor','c', 'EdgeColor', 'none','FaceAlpha', 0.1);  % cyan
        contourf(x1_grid, x2_grid, W_mask,   [0.5 1], 'FaceColor','k', 'EdgeColor', 'none','FaceAlpha', 0.11);  % black
        contourf(x1_grid, x2_grid, Vhat_mask,[0.5 1], 'FaceColor','#FFB9BA', 'EdgeColor', 'none','FaceAlpha', 0.4);  % red

    % figure
        contour(x1_grid, x2_grid, H_values   , [1 1]); hold on;
        contour(x1_grid, x2_grid, Wvalues    , [1 1],'Color', '#525252',LineWidth=0.1); hold on;
        contour(x1_grid, x2_grid, Vhat_values, [1 1], 'r');
        legend('$S_h$',"W(x)","$\hat{V}$",Interpreter='latex')
        xlabel("$x_1$",Interpreter='latex',FontSize=15); 
        ylabel("$x_2$",Interpreter='latex',FontSize=15);
        title('Phase Plane');

    [x01_grid,x02_grid]=meshgrid(linspace(-range_x,range_x,5),linspace(-range_y,range_y,14));
    for i=1:size(x01_grid,1)
        for j=1:size(x01_grid,2)
            x = [x01_grid(i,j); x02_grid(i,j); zeros(n-2,1)];

            for k=1:N/Ts
                [~, u_sat, ~] = Saturation(x, K1, K2);

                x_prev = x;
                x = A*x + B*u_sat;
                dx = x(1) - x_prev(1);
                dy = x(2) - x_prev(2);

                if abs(x(1,:))>range_x || abs(x(2,:))>range_y || all(abs(x)<1e-4)
                    break;
                end
    
                h=quiver(x_prev(1), x_prev(2), dx, dy, Color="#A02BA8", LineWidth=0.7,MaxHeadSize=0.5); hold on;
                set(get(get(h,'Annotation'),'LegendInformation'),IconDisplayStyle='off');
            end     
        end
    end

end