function [A,B,Ts]=Discretize(Ac,Bc,Ts,do_plot)

    if nargin < 4
        do_plot = 0;
    end
    
    [n,m] = size(Bc);
    sysc=ss(Ac,Bc,eye(n), zeros(1,m));

    sysd=c2d(sysc,Ts);
    
    A=sysd.A;
    B=sysd.B;
    
    if do_plot
        figure;
        bodemag(sysc); hold on;
        bodemag(sysd);
        legend("Continuos", "Discretized");
    end

end