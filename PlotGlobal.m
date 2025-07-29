function PlotGlobal(Ac, Bc, Ts, sol, x0, N, range)

    [A,B] = Discretize(Ac,Bc,Ts,0);

    dt = 0.01;

    range_x = range(1);
    range_y = range(2);
        
    K1 = sol.K1;
    K2 = sol.K2;
    P  = sol.P ;
 
    [n,~] = size(B);
    
    x = x0;
    traj = [];
    V=0;
    
    for k = 1:N/Ts
        [u, u_sat, dz] = Saturation(x, K1, K2);
        
        Vant=V;
        V = [x; dz]'*P*[x; dz];
        dV=V-Vant;
    
        data = [x; V; dV; u; u_sat];
        traj = [traj data];
    
        x = A*x + B*u_sat;
    
    end
    

    x = x0;
    traj_c = [];
    
    for k = 1:N/dt
        traj_c = [traj_c x];
    
        [~, u_sat, ~] = Saturation(x, K1, K2);
    
        dx = Ac*x + Bc*u_sat;
        x = x + dt*dx;
    end
    
    [x1_grid, x2_grid] = meshgrid(linspace(-range_x,range_x,800),linspace(-range_y,range_y,800));
    
    Vgrid = zeros(size(x1_grid));
    for i = 1:size(x1_grid,1)
        for j = 1:size(x1_grid,2)
    
            x = [x1_grid(i,j); x2_grid(i,j)];
    
            [~, ~, dz] = Saturation(x, K1, K2);
    
            Vgrid(i,j) = [x; dz]'*P*[x; dz];
        end
    end
    
    figure
    subplot 211
        plot(0:Ts:N-Ts,traj(n+1,:),'.-');
        legend("V")
        titulo=sprintf('Discret - Ts=%.2f',Ts);
        title(titulo)
    subplot 212
        plot(2*Ts:Ts:N,traj(n+2, 2:end),'.-')
        legend("△V")
    
    figure
        plot(traj(1,:) ,traj(2,:),'.-'); hold on
        plot(traj(1,1) ,traj(2,1),'xr','LineWidth',2); hold on
        legend("trajectory","start")
        xlabel("x1"); ylabel("x2")
        axis equal
    
    figure
        plot3(traj(1,:),traj(2,:),traj(3,:),'.-r','MarkerSize',10);hold on;
        plot3(traj(1,1),traj(2,1),traj(3,1),'xr','LineWidth',2);hold on;
        mesh(x1_grid,x2_grid,Vgrid)
        xlabel("x1"); ylabel("x2"); zlabel("V")
        legend("trajet states","start","V(x)")
    
    
    figure
    for i=1:n
        subplot(n,1,i)
            plot(0:Ts:N-Ts,traj(i,:),'.'); hold on;
            plot(0:dt:N-dt,traj_c(i,:))
            legend("discret","continu")
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