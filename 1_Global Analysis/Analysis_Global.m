function [feas,sol]=Analysis_Global(A,B,K1,K2,do_print,relax,options)
    
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

    if nargin < 7
        options = [0,0,0,0,0];
    end

    %% Variables
    
    setlmis([]);

    P0 = lmivar(1,[n,1]);
    P1 = lmivar(2,[n,m]);
    P2 = lmivar(1,[m,1]);
    
    T  = lmivar(1,[m,1]);
    T3 = lmivar(1,[m,1]);
    T2 = lmivar(1,[m,1]);
    T1 = lmivar(1,[m,1]);

    %% Conditions 

    M1 = A+B*K1;
    M2 = K2-Im;
    
    %Condition 1 (V > 0)
    %|  P0             *       | >relax_cond1*I>0
    %|P1'-TK1   P2-He(T(K2-Im))|    
    lmiterm([-1 1 1 P0],1,1);
    lmiterm([-1 1 1 0 ],-relax_cond1*In);
    
    lmiterm([-1 2 1 -P1], 1,1);
    lmiterm([-1 2 1  T ],-1,K1);
    
    lmiterm([-1 2 2 P2], 1,1);
    lmiterm([-1 2 2 T ],-1,M2,'s');
    lmiterm([-1 2 2 0 ],-relax_cond1*Im);
    
    %Condition 2 (△V<0)
    lmiterm([2 1 1 P0],M1',M1);  
    lmiterm([2 1 1 P0],-1,1);
    lmiterm([2 1 1 0 ],relax_cond2*In);
    
    lmiterm([2 2 1  P0],(B*M2)',M1);
    lmiterm([2 2 1 -P1],-1,1);
    lmiterm([2 2 1  T1], 1,K1);
    lmiterm([2 2 1  T3],-1,K1*(M1-In));
    
    lmiterm([2 3 1 -P1],1,M1);
    lmiterm([2 3 1  T2],1,K1*(M1-In));
    
    lmiterm([2 2 2 P0],(B*M2)',B*M2);
    lmiterm([2 2 2 P2],-1,1);
    lmiterm([2 2 2 T1], 1,M2,'s');
    lmiterm([2 2 2 T3],-1,(K1*B-Im)*M2,'s');
    lmiterm([2 2 2 0 ],relax_cond2*Im);
    
    lmiterm([2 3 2 -P1],1,B*M2);
    lmiterm([2 3 2  T2],1,(K1*B-Im)*M2);
    lmiterm([2 3 2  T3],-M2',1);
    
    lmiterm([2 3 3 P2],1,1);
    lmiterm([2 3 3 T2],1,M2,'s');
    lmiterm([2 3 3 0 ],relax_cond2*Im);
    
    %P0>0
    lmiterm([-3 1 1 P0],1,1);
    
    %T's>0
    lmiterm([-4 1 1 T ],1,1);
    lmiterm([-5 1 1 T3],1,1);
    lmiterm([-6 1 1 T2],1,1);
    lmiterm([-7 1 1 T1],1,1);

    %% Solve
    LMIsys = getlmis;
    feas = feasp(LMIsys);  % without optimisation

    if round(feas,3) < 0

        [~,LMIsol] = feasp(LMIsys,options);
        
        sol.P0 = dec2mat(LMIsys, LMIsol, P0);
        sol.P1 = dec2mat(LMIsys, LMIsol, P1);
        sol.P2 = dec2mat(LMIsys, LMIsol, P2);
        sol.T  = dec2mat(LMIsys, LMIsol, T );
        sol.T3 = dec2mat(LMIsys, LMIsol, T3);
        sol.T2 = dec2mat(LMIsys, LMIsol, T2);
        sol.T1 = dec2mat(LMIsys, LMIsol, T1);
    
        sol.P  = [sol.P0, sol.P1; sol.P1', sol.P2];

        sol.H1 = 0;
        sol.H2 = 0;
    
        %% Check Condition 1:
        
        cond1_11 = sol.P0;
        cond1_21 = sol.P1' - sol.T*K1;
        cond1_22 = sol.P2  - sol.T*M2  - (sol.T*M2)';
        
        sol.cond1 = [cond1_11, cond1_21';
                     cond1_21, cond1_22];
        sol.eig1  = eig(sol.cond1);

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
        sol.eig2  = eig(sol.cond2);

        %% Print conditions

        if do_print
            disp("Eigenvalue Condition 1 (V>0)"); disp(sol.eig1);
            disp("Eigenvalue Condition 2 (△V)"); disp(sol.eig2);
        end

    else
        disp("Unfeasible system. Feas:");
        disp(feas);

        sol = 0;
    end

end
