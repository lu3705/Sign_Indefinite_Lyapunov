function PlotRegional(A,B,Ts,sol, x0, N,r, range)

    range_x=range(1);
    range_y=range(2);
        
    K1=sol.K1;
    K2=sol.K2;
    P=sol.P;
    Phat=sol.Phat;
    H1=sol.H1;
    H2=sol.H2;
 
    [n,~] = size(B);
    traj=[];
    V=0;
    x=x0;
    
    for k=1:N/Ts
        [u, u_sat, dz] = Saturation(x, K1, K2);
        
        Vant=V;
        V = [x; dz]'*P*[x; dz];
        dV=V-Vant;
    
        data = [x; V; dV; u; u_sat];
        traj = [traj data];
    
        x = A*x + B*u_sat;
    end
    
    [x1_grid,x2_grid]=meshgrid(linspace(-range_x,range_x,800),linspace(-range_y,range_y,800));
    
    Vgrid=zeros(size(x1_grid));
    Wvalues=zeros(size(x1_grid));
    H_values=zeros(size(x1_grid));
    Vhat_values=zeros(size(x1_grid));
    for i=1:size(x1_grid,1)
        for j=1:size(x1_grid,2)
            x = [x1_grid(i,j); x2_grid(i,j)];
    
            [~, ~, dz] = Saturation(x, K1, K2);
    
            Vgrid(i,j) = [x; dz]'*P*[x; dz];
    
            Vgrid(i,j)=[x; dz]'*P*[x; dz];
            Vhat_values(i,j)=x'*Phat*x;
    
            hx=H1*x+H2*dz;
    
            if norm(hx,inf)<=1
                H_values(i,j) = norm(hx, inf); 
                Wvalues(i,j)=min(Vgrid(i,j), 1);
            else
                Wvalues(i,j)=1;
                H_values(i,j) = 1; 
                Vhat_values(i,j)=1;
            end
        end
    end

    figure
    subplot 211
        plot(0:Ts:N-Ts,traj(n+1,:),'.-');
        legend("V")
        titulo=sprintf('Discret-Desconsidera sat no continuo-Ts=%.2f',Ts);
        title(titulo)
    subplot 212
        plot(2*Ts:Ts:N,traj(n+2, 2:end),'.-')
        legend("△V")

    figure
        plot(traj(1,:)  ,traj(2,:),'.-'); hold on
        plot(traj(1,1)  ,traj(2,1),'xr','LineWidth',2); hold on
        contour(x1_grid, x2_grid, H_values, [1 1]); hold on;
        contour(x1_grid, x2_grid, Wvalues, [1 1], 'k'); hold on;
        contour(x1_grid, x2_grid, Vhat_values, [1 1], 'r');
        label=sprintf('$S_h : %s = %.2f$', inputname(7), r);
        legend("trajectory","start",label,"W(x)","$\hat{V}$",'Interpreter','latex')
        xlabel("x1"); ylabel("x2")
        if inputname(7)=="tal"
            label2=sprintf("Análise Regional Discreto - $ T_s = %.2f$", Ts);
        else 
            label2=sprintf("Sintese Regional Discreto - $ T_s = %.2f$", Ts);
        end
        title(label2,'Interpreter','latex');

    figure
        plot3(traj(1,:),traj(2,:),traj(3,:),'.r','MarkerSize',10)
        xlabel("x1"); ylabel("x2"); zlabel("V")
        hold on
        mesh(x1_grid,x2_grid,Vgrid)
        legend("trajet states","V")

    figure
    for i=1:n
        subplot(n,1,i)
        plot(0:Ts:N-Ts,traj(i,:),'.-')
        legend("discret")
        xlabel("time")
        titulo=sprintf('x%d ',i);
        title(titulo)
    end

    figure
        plot(0:Ts:N-Ts,traj(n+3,:),'.-'); hold on;
        plot(0:Ts:N-Ts,traj(n+4,:),'.-'); 
        legend("u","u_sat")
        xlabel("time")
        title("control signal")
end