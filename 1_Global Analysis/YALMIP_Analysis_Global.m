function [feas,sol]=YALMIP_Analysis_Global(A,B,K1,K2,do_print,relax)
    
    sol.K1 = K1;
    sol.K2 = K2;

    [n,m] = size(B);
    Im = eye(m);
    In = eye(n);

    if nargin < 6
        relax = [1e-20,1e-20];
    end

    relax_cond1 = relax(1);
    relax_cond2 = relax(2);

    %% Variables
    
    setlmis([]);

    P0=sdpvar(n, n);
    P1=sdpvar(n, m, 'full');
    P2=sdpvar(m, m);
    
    T =sdpvar(m, m);
    T1=sdpvar(m, m);
    T2=sdpvar(m, m);
    T3=sdpvar(m, m);

    M1 = A+B*K1;
    M2 = K2-Im;

    constraints = [];

    %% Condição 1 – V(x) > 0
    cond1_11 = P0;
    cond1_21 = P1' - T*K1;
    cond1_22 = P2  - T*M2  - (T*M2)';
    
    cond1 = [cond1_11, cond1_21';
             cond1_21, cond1_22];

    constraints = [constraints, cond1 >= relax_cond1*eye(size(cond1))];
    
    %% Condition 2 (△V<0)

    cond2_11 = M1'*P0*M1 - P0;
    cond2_21 = (B*M2)'*P0*M1   -  P1'  +   T1*K1  -  T3*K1*(M1-In);
    cond2_31 = P1'*M1    +  T2*K1*(M1-In);
    cond2_22 = (B*M2)'*P0*(B*M2) - P2  +  T1*M2   +  (T1*M2)'  -  T3*(K1*B-Im)*M2  -  (T3*(K1*B-Im)*M2)';
    cond2_32 = P1'*B*M2   +   T2 * (K1*B-Im)*M2  -  M2'*T3;
    cond2_33 = P2  +  T2*M2 +  (T2*M2)';
    
    cond2 = [cond2_11, cond2_21', cond2_31';
             cond2_21, cond2_22 , cond2_32';
             cond2_31, cond2_32 , cond2_33];
    
    constraints = [constraints, cond2 <= -relax_cond2*eye(size(cond2))];
    
    %% Condições de positividade
    constraints = [constraints,
                   P0 >= 1e-4*In,
                   T  >= 1e-4*Im,
                   T1 >= 1e-4*Im,
                   T2 >= 1e-4*Im,
                   T3 >= 1e-4*Im];

    %% Solve
    ops = sdpsettings('solver', 'sdpt3');
    diagnostics = optimize(constraints, [], ops);

    feas = diagnostics.problem;
    if feas == 0

        sol.P0 = value(P0);
        sol.P1 = value(P1);
        sol.P2 = value(P2);
        sol.T  = value(T );
        sol.T1 = value(T1);
        sol.T2 = value(T2);
        sol.T3 = value(T3);

    
        sol.P  = [sol.P0, sol.P1; sol.P1', sol.P2];

        sol.H1=0;
        sol.H2=0;
    
        %% Check Condition 1:
        
        cond1_11 = sol.P0;
        cond1_21 = sol.P1' - sol.T*K1;
        cond1_22 = sol.P2  - sol.T*M2  - (sol.T*M2)';
        
        sol.cond1 = [cond1_11, cond1_21';
                     cond1_21, cond1_22];
        sol.eig1 = eig(sol.cond1);

        %% Check Condition 2:
        cond2_11 = M1'* sol.P0 *M1 - sol.P0;
        cond2_21 = (B*M2)' * sol.P0 * M1   -  sol.P1'  +   sol.T1*K1  -  sol.T3*K1*(M1-In);
        cond2_31 = sol.P1'*M1  +  sol.T2*K1*(M1-In);
        cond2_22 = (B*M2)' * sol.P0 * (B*M2) - sol.P2  +  sol.T1*M2   +  (sol.T1*M2)'  -  sol.T3*(K1*B-Im)*M2  -  (sol.T3*(K1*B-Im)*M2)';
        cond2_32 = sol.P1'*B*M2   +   sol.T2 * (K1*B-Im)*M2  -  M2'*sol.T3;
        cond2_33 = sol.P2  +  sol.T2*M2 +  (sol.T2*M2)';
        
        sol.cond2 = [cond2_11, cond2_21', cond2_31';
                     cond2_21, cond2_22 , cond2_32';
                     cond2_31, cond2_32 , cond2_33];
        sol.eig2 = eig(sol.cond2);

        %% Print conditions

        if do_print
            disp("Eigenvalue Condition 1 (V>0)"); disp(sol.eig1);
            disp("Eigenvalue Condition 2 (△V)"); disp(sol.eig2);
        end

    else
        disp("Sistema infactível. Feas:");
        disp(feas);

        sol=0;
    end

end