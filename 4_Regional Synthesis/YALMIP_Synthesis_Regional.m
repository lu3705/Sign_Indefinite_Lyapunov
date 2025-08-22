function [feas, sol] = YALMIP_Synthesis_Regional(A, B, r, do_print,relax)

    [n, m] = size(B);
    In = eye(n);
    Im = eye(m);

    if nargin < 5
        relax = [1e-20, 1e-20, 1e-20, 1e-20];
    end

    relax_cond1=relax(1);
    relax_cond2=relax(2);
    relax_r =relax(3);
    relax_Vhat=relax(4);

    %% Variables

    R   = sdpvar(n, n);

    Q0  = sdpvar(n, n);
    Q1  = sdpvar(n, m, 'full');
    Q2  = sdpvar(m, m);

    S   = sdpvar(m, m);
    M   = sdpvar(n, n, 'full');

    Y1  = sdpvar(m, n, 'full');
    Y2  = sdpvar(m, m, 'full');

    Phat= sdpvar(n, n);

    Z1  = sdpvar(m, n, 'full');
    Z2  = sdpvar(m, m, 'full');
    
    constraints = [];

    %% Condição 1 – V(x) > 0

    cond1_11= Q0;
    cond1_21= Q1' + Z1 - Y1;
    cond1_22= Q2  + Z2 + Z2' -  Y2  -  Y2' + 2*S;

    for i = 1:m
        line=zeros(m,1); line(i)=1;

        cond1_31= line'*Z1;
        cond1_32= line'*Z2;
        cond1_33= 1;
        
        cond1=[cond1_11, cond1_21', cond1_31';
               cond1_21, cond1_22 , cond1_32';
               cond1_31, cond1_32 , cond1_33];
        
        constraints = [constraints, cond1 >= relax_cond1*eye(size(cond1))];
    end

    %% Condição 2 – △V(x) < 0
    cond2_11 = -Q0;
    cond2_21 = -Q1' +  2*Y1 -  Z1; 
    cond2_31 = A*M  +  B*Y1;
    cond2_41 = -Y1;
    cond2_22 = -Q2  +  2*Y2 + 2*Y2' - 4*S  - Z2 - Z2';
    cond2_32 = -Y1' +  B*Y2 - B*S;
    cond2_42 = -Y2  -  Y2'  + 2*S;
    cond2_33 =  Q0  -  M    -  M';
    cond2_43 =  Q1' +  Y1;
    cond2_44 =  Q2  +  Y2   +  Y2'  - 2*S;
    
    cond2=[cond2_11, cond2_21', cond2_31' , cond2_41';
           cond2_21, cond2_22 , cond2_32' , cond2_42';
           cond2_31, cond2_32 , cond2_33  , cond2_43';
           cond2_41, cond2_42 , cond2_43  , cond2_44];

    constraints = [constraints, cond2 <= -relax_cond2*eye(size(cond2))];

    %% Condição 3 – Alocação de polos (decay rate)
    cond3_11 = -r^2*(M  + M'- R);
    cond3_21 = A*M + B*Y1; 
    cond3_22 = -R;
    
    cond3=[cond3_11, cond3_21';
           cond3_21, cond3_22];

    constraints = [constraints, cond3 <= -relax_r*eye(size(cond3))];

    %% Condição 4 – Maximizar região (Phat)
    cond4_11 =  M   + M'  - Q0;
    cond4_21 = In;
    cond4_31 = -Q1' - Y1 ; %+ Z1
    cond4_22 = Phat;
    cond4_32 = zeros(m,n);
    cond4_33 =-Q2   + 2*S - Y2 - Y2' ; %+ Z2+ Z2'
    
    cond4=[cond4_11, cond4_21', cond4_31';
           cond4_21, cond4_22 , cond4_32';
           cond4_31, cond4_32 , cond4_33];

    constraints = [constraints, cond4 >= relax_Vhat*eye(size(cond4))];

    % Condições de positividade
    constraints = [constraints,
                   Q0   >= 1e-4*In,
                   S    >= 1e-4*Im,
                   R    >= 1e-4*In,
                   Phat >= 1e-4*In
    ];

    %% Solve
    ops = sdpsettings('solver', 'sdpt3');
    diagnostics = optimize(constraints, trace(Phat), ops);

    feas = diagnostics.problem;
    if feas == 0
        sol.Q0 = value(Q0);
        sol.Q1 = value(Q1);
        sol.Q2 = value(Q2);
        sol.S  = value(S);
        sol.M  = value(M);
        sol.Y1 = value(Y1);
        sol.Y2 = value(Y2);
        sol.Z1 = value(Z1);
        sol.Z2 = value(Z2);
        sol.R  = value(R);
        sol.Phat = value(Phat);

        sol.P0  = inv(sol.M')*sol.Q0*inv(sol.M);
        sol.P1  = inv(sol.M')*sol.Q1*inv(sol.S);
        sol.P2  = inv(sol.S )*sol.Q2*inv(sol.S);

        sol.P   = [sol.P0, sol.P1; sol.P1', sol.P2];
    
        sol.H1  = sol.Z1*inv(sol.M);
        sol.H2  = sol.Z2*inv(sol.S);

        sol.K1  = sol.Y1*inv(sol.M);
        sol.K2  = sol.Y2*inv(sol.S);

        %% Check Condition 1:
        cond1_11= sol.Q0;
        cond1_21= sol.Q1' + sol.Z1 - sol.Y1;
        cond1_22= sol.Q2  + sol.Z2 + sol.Z2' -  sol.Y2  -  sol.Y2' + 2*sol.S;

        for i = 1:m
            line=zeros(m,1); line(i)=1;

            cond1_31= line'*sol.Z1;
            cond1_32= line'*sol.Z2;
            cond1_33= 1;
            
            sol.cond1{i}=[cond1_11, cond1_21', cond1_31';
                          cond1_21, cond1_22 , cond1_32';
                          cond1_31, cond1_32 , cond1_33];
            sol.eig1{i} =eig(sol.cond1{i});
        end

        %% Check Condition 2:
        cond2_11 = -sol.Q0;
        cond2_21 = -sol.Q1' +  2*sol.Y1 -  sol.Z1; 
        cond2_31 = A*sol.M  +  B*sol.Y1;
        cond2_41 = -sol.Y1;
        cond2_22 = -sol.Q2  +  2*sol.Y2 + 2*sol.Y2' - 4*sol.S  - sol.Z2 - sol.Z2';
        cond2_32 = -sol.Y1' +  B*sol.Y2 - B*sol.S;
        cond2_42 = -sol.Y2  -  sol.Y2'  + 2*sol.S;
        cond2_33 =  sol.Q0  -  sol.M    -  sol.M';
        cond2_43 =  sol.Q1' +  sol.Y1;
        cond2_44 =  sol.Q2  +  sol.Y2   +  sol.Y2'  - 2*sol.S;
        
        sol.cond2=[cond2_11, cond2_21', cond2_31' , cond2_41';
                   cond2_21, cond2_22 , cond2_32' , cond2_42';
                   cond2_31, cond2_32 , cond2_33  , cond2_43';
                   cond2_41, cond2_42 , cond2_43  , cond2_44];
        sol.eig2=eig(sol.cond2);
        
        %% Check Condition 3 (r):
        cond3_11 = -r^2*(sol.M  + sol.M'- sol.R);
        cond3_21 = A*sol.M + B*sol.Y1; 
        cond3_22 = -sol.R;
        
        sol.cond3=[cond3_11, cond3_21';
                   cond3_21, cond3_22];
        sol.eig3=eig(sol.cond3);

        %% Check Condition 4 (Vhat)
        cond4_11 =  sol.M   + sol.M'  - sol.Q0;
        cond4_21 = In;
        cond4_31 = -sol.Q1' - sol.Y1;  %+ sol.Z1
        cond4_22 = sol.Phat;
        cond4_32 = zeros(m,n);
        cond4_33 =-sol.Q2   + 2*sol.S - sol.Y2 - sol.Y2'; %+ sol.Z2+ sol.Z2'
        
        sol.cond4=[cond4_11, cond4_21', cond4_31';
                   cond4_21, cond4_22 , cond4_32';
                   cond4_31, cond4_32 , cond4_33];
        sol.eig4=eig(sol.cond4);

        %% Print conditions

        if do_print
            for i=1:m
                disp("Eigenvalue Condition 1 (V>0)"); disp(sol.eig1{i});
            end
            disp("Eigenvalue Condition 2 (△V)");  disp(sol.eig2);
            disp("Eigenvalue Condition 3 (r)");    disp(sol.eig3);
            disp("Eigenvalue Condition 4 (Vhat)"); disp(sol.eig4);
        end

    else
        disp("Sistema infactível. Feas:");
        disp(feas);
    end

end