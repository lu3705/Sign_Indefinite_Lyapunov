function [feas,sol]=Analysis_Regional(A,B,K1,K2,tal,do_print,relax,options)
    
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

    if nargin < 8
        options = [0,0,0,0,0];
    end

    %% Variables

    setlmis([])
    
    P0=lmivar(1,[n,1]);
    P1=lmivar(2,[n,m]);
    P2=lmivar(1,[m,1]);
    
    T1=lmivar(1,[m,1]);
    T2=lmivar(1,[m,1]);
    T3=lmivar(1,[m,1]);
    
    Phat=lmivar(1,[n,1]);
    That=lmivar(1,[m,1]);
    
    Z1=lmivar(2,[m,n]);
    Z2=lmivar(2,[m,m]);

    %% Conditions 

    M1=A+B*K1;
    M2=K2-Im;
    
    %Condição 1 (V>0 e inclusão)     (mesmo Sophie)
    % |     P0                *               *  | 
    % |P1'+Z1-T1*K1   P2+He(Z2-T1(K2-Im))     *  | >relax_cond1*I>0
    % |     Z1i              Z2i           tal*ti| 
    for i = 1:m
        idx =-8-i;
        diag=zeros(m,m);  diag(i,i)=1;
        line=zeros(m,1); line(i)=1;
    
        lmiterm([idx 1 1 P0],1,1);
        lmiterm([idx 1 1 0],-relax_cond1*In);
    
    
        lmiterm([idx 2 1 -P1], 1,1);
        lmiterm([idx 2 1  Z1], 1,1);
        lmiterm([idx 2 1  T1],-1,K1);
    
        lmiterm([idx 3 1 Z1],line',1);
    
        lmiterm([idx 2 2 P2], 1, 1);
        lmiterm([idx 2 2 Z2], 1, 1,'s');
        lmiterm([idx 2 2 T1],-1,M2,'s');
        lmiterm([idx 1 1 0],-relax_cond1*Im);
    
        lmiterm([idx 3 2 Z2], line', 1);
    
        lmiterm([idx 3 3 T1], tal*diag, 1);
        lmiterm([idx 1 1 0],-relax_cond1*In);
    end
    
    %Condition 2 (△V<0)
    lmiterm([1 1 1 P0],M1',M1);  
    lmiterm([1 1 1 P0],-1,1); 
    lmiterm([1 1 1 0],relax_cond2*In)
     
    lmiterm([1 2 1  P0],(B*M2)',M1);
    lmiterm([1 2 1  T3],-1,K1*(M1-In));
    lmiterm([1 2 1  T1], 1,K1);
    lmiterm([1 2 1  Z1],-1,1);
    lmiterm([1 2 1 -P1],-1,1);
    
    lmiterm([1 3 1 -P1],1,M1);
    lmiterm([1 3 1  T2],1,K1*(M1-In));
    
    lmiterm([1 2 2 P0],(B*M2)',B*M2);
    lmiterm([1 2 2 T3],-1,(K1*B-Im)*M2,'s');
    lmiterm([1 2 2 T1], 1,M2,'s');
    lmiterm([1 2 2 Z2],-1, 1,'s');
    lmiterm([1 2 2 P2],-1, 1);
    lmiterm([1 2 2 0],relax_cond2*Im)
    
    lmiterm([1 3 2 -P1],1,B*M2);
    lmiterm([1 3 2  T2],1,(K1*B-Im)*M2);
    lmiterm([1 3 2  T3],M2',-1);
    
    lmiterm([1 3 3 P2],1,1);
    lmiterm([1 3 3 T2],1,M2,'s');
    lmiterm([1 3 3 0],relax_cond2*Im)
  
    
    % Condição 3 (maximiza área)
    % |Phat-P0       -(P1+K1'That)     | >relax_Vhat*I>0
    % |   *         -P2-He(That(K2-Im))|
    lmiterm([-2 1 1 Phat], 1,1);
    lmiterm([-2 1 1  P0 ],-1,1);
    lmiterm([-2 1 1  0 ],-relax_Vhat*In);
    
    lmiterm([-2 2 1  -P1],-1, 1);
    lmiterm([-2 2 1 That],-1,K1);
    
    lmiterm([-2 2 2  P2 ],-1, 1);
    lmiterm([-2 2 2 That],-1,M2,'s');
    lmiterm([-2 2 2  0 ],-relax_Vhat*Im);
    
    % P0>0
    lmiterm([-3 1 1 P0], 1,1);
    
    %Phat>0
    lmiterm([-4 1 1 Phat],1,1);
    
    %T>0
    lmiterm([-5 1 1 T2],1,1);
    lmiterm([-6 1 1 T3],1,1); 
    lmiterm([-7 1 1 That],1,1);
    
    %Condição 0 (T-tal>0)            (mesmo Sophie)
    lmiterm([-8 1 1 T1], 1 , 1);  
    lmiterm([-8 1 1 0 ],-tal*Im);

    %% Solve
    LMIsys=getlmis;
    feas = feasp(LMIsys,options);

    if round(feas,3) < 0

        nx = decnbr(LMIsys); 
        c = zeros(nx,1);
        
        for i = 1:nx 
            [Phati] = defcx(LMIsys,i,Phat); 
            c(i) = trace(Phati);
        end
        
        [~, LMIsol] = mincx(LMIsys, c,options);
        
        sol.P0    = dec2mat(LMIsys, LMIsol, P0);
        sol.P1    = dec2mat(LMIsys, LMIsol, P1);
        sol.P2    = dec2mat(LMIsys, LMIsol, P2);
        sol.T1    = dec2mat(LMIsys, LMIsol, T1);
        sol.T2    = dec2mat(LMIsys, LMIsol, T2);
        sol.T3    = dec2mat(LMIsys, LMIsol, T3);
        sol.Z1    = dec2mat(LMIsys, LMIsol, Z1);
        sol.Z2    = dec2mat(LMIsys, LMIsol, Z2);
        sol.Phat  = dec2mat(LMIsys, LMIsol, Phat);
        sol.That  = dec2mat(LMIsys, LMIsol, That);
    
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
        disp("Sistema infactível. Feas:");
        disp(feas);

        sol=0;
    end

end