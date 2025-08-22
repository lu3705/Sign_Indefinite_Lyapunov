function [feas,sol]=YALMIP_Analysis_Regional(A,B,K1,K2,tal,do_print,relax)
    
    sol.K1=K1;
    sol.K2=K2;

    [n,m] = size(B);
    Im = eye(m);
    In = eye(n);

    if nargin < 6
        relax = [1e-20,1e-20,1e-20];
    end

    relax_cond1=relax(1);
    relax_cond2=relax(2);
    relax_Vhat =relax(3);


    %% Variables
    
    P0=sdpvar(n, n);
    P1=sdpvar(n, m, 'full');
    P2=sdpvar(m, m);
    
    T1=sdpvar(m, m);
    T2=sdpvar(m, m);
    T3=sdpvar(m, m);
    
    Phat=sdpvar(n, n);
    That=sdpvar(m, m);
    
    Z1  = sdpvar(m, n, 'full');
    Z2  = sdpvar(m, m, 'full');

    M1=A+B*K1;
    M2=K2-Im;
    
    constraints = [];

    %% Condição 1 – V(x) > 0
    cond1_11= P0;
    cond1_21= P1' + Z1 - T1*K1;
    cond1_22= P2  + Z2 + Z2' - T1*M2 - (T1*M2)';
    
    for i = 1:m
        line=zeros(m,1); line(i)=1;
        diag =zeros(m,m); diag(i,i)=1;
        cond1_31= line'*Z1;
        cond1_32= line'*Z2;
        cond1_33= T1*tal*diag;
        
        cond1=[cond1_11, cond1_21', cond1_31';
               cond1_21, cond1_22 , cond1_32';
               cond1_31, cond1_32 , cond1_33];

        constraints = [constraints, cond1 >= relax_cond1*eye(size(cond1))];
    end

    
    %% Condição 2 – △V(x) < 0
    cond2_11 = M1'* P0 *M1 - P0;
    cond2_21 = (B*M2)' * P0 * M1   -  P1'  +   T1*K1   -   Z1  -  T3*K1*(M1-In);
    cond2_31 = P1'*M1  +  T2*K1*(M1-In);
    cond2_22 = (B*M2)' * P0 * (B*M2) - P2  +  T1*M2  + (T1*M2)'  -  Z2   -  Z2'  -   T3*(K1*B-Im)*M2  -  (T3*(K1*B-Im)*M2)';
    cond2_32 = P1'*B*M2   +   T2 * (K1*B-Im)*M2  -  M2'*T3;
    cond2_33 = P2  +  T2*M2 +  (T2*M2)';

    cond2=[cond2_11, cond2_21', cond2_31';
               cond2_21, cond2_22 , cond2_32';
               cond2_31, cond2_32 , cond2_33];
  
    constraints = [constraints, cond2 <= -relax_cond2*eye(size(cond2))];

    %% Condição 3 – Max area (Phat)
    cond3_11= Phat-P0;
    cond3_21=-(P1' + That*K1);
    cond3_22=-(P2  + That*M2 + (That*M2)' );

    cond3=[cond3_11, cond3_21';
           cond3_21, cond3_22];

    constraints = [constraints, cond3 >= relax_Vhat*eye(size(cond3))];

    
    %% Condições de positividade
    constraints = [constraints,
                   P0   >= 1e-4*In,
                   Phat >= 1e-4*In,
                   T1   >= tal *Im,
                   T2   >= 1e-4*Im,
                   T3   >= 1e-4*Im,
                   That >= 1e-4*Im
    ];

    %% Solve
    ops = sdpsettings('solver', 'sdpt3');
    diagnostics = optimize(constraints, trace(Phat), ops);

    feas = diagnostics.problem;
    if feas == 0
        sol.P0 = value(P0);
        sol.P1 = value(P1);
        sol.P2 = value(P2);
        sol.T1 = value(T1);
        sol.T2 = value(T2);
        sol.T3 = value(T3);
        sol.That = value(That);
        sol.Z1 = value(Z1);
        sol.Z2 = value(Z2);
        sol.Phat = value(Phat);
    
        sol.P     = [sol.P0, sol.P1; sol.P1', sol.P2];
    
        sol.H1    = inv(sol.T1)*sol.Z1;
        sol.H2    = inv(sol.T1)*sol.Z2;
    
        %% Check Condition 1:
        cond1_11= sol.P0;
        cond1_21= sol.P1' + sol.Z1 - sol.T1*K1;
        cond1_22= sol.P2  + sol.Z2 + sol.Z2' -  sol.T1*M2  -  (sol.T1*M2)';
        
        for i = 1:m
            line=zeros(m,1); line(i)=1;
            diag =zeros(m,m); diag(i,i)=1;
            cond1_31= line'*sol.Z1;
            cond1_32= line'*sol.Z2;
            cond1_33= sol.T1*tal*diag;
            
            sol.cond1{i}=[cond1_11, cond1_21', cond1_31';
                          cond1_21, cond1_22 , cond1_32';
                          cond1_31, cond1_32 , cond1_33];
            sol.eig1{i} =eig(sol.cond1{i});
    
        end

        %% Check Condition 2:
        cond2_11 = M1'* sol.P0 *M1 - sol.P0;
        cond2_21 = (B*M2)' * sol.P0 * M1   -  sol.P1'  +   sol.T1*K1   -   sol.Z1  -  sol.T3*K1*(M1-In);
        cond2_31 = sol.P1'*M1  +  sol.T2*K1*(M1-In);
        cond2_22 = (B*M2)' * sol.P0 * (B*M2) - sol.P2  +  sol.T1*M2  + (sol.T1*M2)'  -  sol.Z2   -  sol.Z2'  -   sol.T3*(K1*B-Im)*M2  -  (sol.T3*(K1*B-Im)*M2)';
        cond2_32 = sol.P1'*B*M2   +   sol.T2 * (K1*B-Im)*M2  -  M2'*sol.T3;
        cond2_33 = sol.P2  +  sol.T2*M2 +  (sol.T2*M2)';

        sol.cond2=[cond2_11, cond2_21', cond2_31';
                   cond2_21, cond2_22 , cond2_32';
                   cond2_31, cond2_32 , cond2_33];
        sol.eig2=eig(sol.cond2);
    
        %% Check Condition 3:
        cond3_11= sol.Phat-sol.P0;
        cond3_21=-(sol.P1' + sol.That*K1);
        cond3_22=-(sol.P2  + sol.That*M2 + (sol.That*M2)' );
    
        sol.cond3=[cond3_11, cond3_21';
                   cond3_21, cond3_22];
        sol.eig3=eig(sol.cond3);


        %% Print conditions

        if do_print
            for i=1:m
                disp("Eigenvalue Condition 1 (V>0)"); disp(sol.eig1{i});
            end
            disp("Eigenvalue Condition 2 (△V)"); disp(sol.eig2);
            disp("Eigenvalue Condition 3 (Vhat)"); disp(sol.eig3);
        end

    else
        disp("Unfeasible system. Feas:");
        disp(feas);

        sol = 0;
    end

end
