function [feas,sol]=Synthesis_Regional(A,B,r,relax,do_print,options)
    
    [n,m] = size(B);
    Im = eye(m);
    In = eye(n);

    if nargin < 4
        relax = [1e-20, 1e-20, 1e-20, 1e-20];
    end

    relax_cond1=relax(1);
    relax_cond2=relax(2);
    relax_r =relax(3);
    relax_Vhat=relax(4);

    if nargin < 6
        options = [0,0,0,0,0];
    end

    %% Variables

    setlmis([]);

    R =lmivar(1,[n,1]);
    
    Q0=lmivar(1,[n,1]);
    Q1=lmivar(2,[n,m]);
    Q2=lmivar(1,[m,1]); 
    
    S =lmivar(1,[m,1]);
    M =lmivar(2,[n,n]);
    
    Y1=lmivar(2,[m,n]);
    Y2=lmivar(2,[m,m]);
    
    Phat=lmivar(1,[n,1]);
    
    Z1=lmivar(2,[m,n]);
    Z2=lmivar(2,[m,m]);
    
    %% Conditions 
    
    %Condição 1 (V>0)  (mesmo Sophie)
    % |Q0     Q1+Z1'-Y1'   Z1i'|
    % | *  Q2+He(Z2-Y2+S)  Z2i'| >0
    % | *        *          1  |
    for i = 1:m
        idx =-8-i;
        linha=zeros(m,1); linha(i)=1;
    
        lmiterm([idx 1 1 Q0], 1,1);
        lmiterm([idx 1 1 0],-relax_cond1*In);
    
        lmiterm([idx 2 1 -Q1], 1,1);
        lmiterm([idx 2 1 Z1], 1,1);
        lmiterm([idx 2 1 Y1],-1,1);
    
        lmiterm([idx 3 1 Z1],linha',1);
    
        lmiterm([idx 2 2 Q2], 1,1);  
        lmiterm([idx 2 2 Z2], 1,1,'s'); 
        lmiterm([idx 2 2 Y2],-1,1,'s'); 
        lmiterm([idx 2 2 S ], 1,1,'s');
        lmiterm([idx 2 2 0],-relax_cond1*Im);
    
        lmiterm([idx 3 2 Z2],linha',1);
    
        lmiterm([idx 3 3 0 ], 1);
        lmiterm([idx 3 3 0],-relax_cond1*Im);
    end
    
    %Condição 2 (deltaV<0)
    % lmiterm([2 1 1 M ], A,1,'s');
    % lmiterm([2 1 1 Y1], B,1,'s');
    lmiterm([2 1 1 Q0],-1,1);
    lmiterm([2 1 1 0],relax_cond2*In);

    lmiterm([2 2 1 -Q1],-1,1);
    lmiterm([2 2 1  Y1], 2,1);
    lmiterm([2 2 1  Z1],-1,1);
    % lmiterm([2 2 1 -Y2], 1,B');
    % lmiterm([2 2 1  S ],-1,B');
    
    lmiterm([2 3 1 M ], A,1);
    lmiterm([2 3 1 Y1], B,1);
    % lmiterm([2 3 1 -M ],-1,1);
    
    lmiterm([2 4 1 Y1],-1,1);
    
    lmiterm([2 2 2 Q2],-1,1);
    lmiterm([2 2 2 Y2], 2,1,'s');
    lmiterm([2 2 2 S ],-2,1,'s');
    lmiterm([2 2 2 Z2],-1,1,'s');
    lmiterm([2 2 2 0],relax_cond2*Im);
    
    lmiterm([2 3 2 -Y1],1,-1);
    lmiterm([2 3 2  Y2],B, 1);
    lmiterm([2 3 2  S ],B,-1);
    
    lmiterm([2 4 2 Y2],-1,1,'s');
    lmiterm([2 4 2 S ], 2,1);
    
    lmiterm([2 3 3 M ],-1,1,'s');
    lmiterm([2 3 3 Q0], 1,1);
    lmiterm([2 3 3 0],relax_cond2*In);
    
    lmiterm([2 4 3  Y1], 1,1);
    lmiterm([2 4 3 -Q1], 1,1);
    
    lmiterm([2 4 4 Y2], 1,1,'s');
    lmiterm([2 4 4 S ],-1,1,'s');
    lmiterm([2 4 4 Q2], 1,1);
    lmiterm([2 4 4 0],relax_cond2*Im);
    
    %Condition 3 (Aumentar velocidade resposta linear) - alocação polos (autovalores)
    % |-alpha(M+M'-R)   * |<-relax_r*I<0
    % |  AM+BY1        -R |
    lmiterm([3 1 1 M ],r^2,-1,'s');
    lmiterm([3 1 1 R ],r^2, 1);
    lmiterm([3 1 1 0],relax_r*In);

    lmiterm([3 2 1 M ], A,1);
    lmiterm([3 2 1 Y1], B,1);

    lmiterm([3 2 2 R ],-1,1);
    lmiterm([3 2 2 0],relax_r*In);
    
    % Condição 4 (maximiza área)
    % |M+M'-Q0     *           *     | >tal*I>0
    % |  In       Phat         *     |
    % |-Q1'-Y1     0    -Q2+2S-Y2-Y2'|
    lmiterm([-4 1 1 M ], 1,1,'s');
    lmiterm([-4 1 1 Q0],-1,1);
    lmiterm([-4 1 1 0],-relax_Vhat*In);

    lmiterm([-4 2 1 0],In);

    lmiterm([-4 3 1 -Q1],-1,1);
    lmiterm([-4 3 1  Y1],-1,1);

    lmiterm([-4 2 2 Phat],1,1);
    lmiterm([-4 2 2 0],-relax_Vhat*In);

    lmiterm([-4 3 2 0],zeros(m,n));

    lmiterm([-4 3 3 Q2],-1,1);
    lmiterm([-4 3 3 S ], 1,1,'s');
    lmiterm([-4 3 3 Y2],-1,1,'s');
    lmiterm([-4 3 3 0],-relax_Vhat*Im);
    
    %Q0>0
    lmiterm([-5 1 1 Q0], 1,1);
    
    %S>0
    lmiterm([-6 1 1 S ],1,1);
    
    %R>0
    lmiterm([-7 1 1 R ],1,1);
    
    %Phat>0
    lmiterm([-8 1 1 Phat],1,1);

    %% Solve
    LMIsys = getlmis;
    feas = feasp(LMIsys,options);

    if round(feas,3) < 0

        nx = decnbr(LMIsys); 
        c = zeros(nx,1);
        
        for i = 1:nx 
            [Phati] = defcx(LMIsys,i,Phat); 
            c(i) = trace(Phati);
        end
        
        [~, LMIsol] = mincx(LMIsys, c,options);
        
        sol.Q0  = dec2mat(LMIsys, LMIsol, Q0);
        sol.Q1  = dec2mat(LMIsys, LMIsol, Q1);
        sol.Q2  = dec2mat(LMIsys, LMIsol, Q2);
        sol.S   = dec2mat(LMIsys, LMIsol, S );
        sol.M   = dec2mat(LMIsys, LMIsol, M );
        sol.Y1  = dec2mat(LMIsys, LMIsol, Y1);
        sol.Y2  = dec2mat(LMIsys, LMIsol, Y2);
        sol.Z1  = dec2mat(LMIsys, LMIsol, Z1);
        sol.Z2  = dec2mat(LMIsys, LMIsol, Z2);
        sol.R   = dec2mat(LMIsys, LMIsol, R );
        sol.Phat= dec2mat(LMIsys, LMIsol, Phat);
    
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
        cond4_31 = -sol.Q1' - sol.Y1; 
        cond4_22 = sol.Phat;
        cond4_32 = zeros(m,n);
        cond4_33 =-sol.Q2   + 2*sol.S - sol.Y2 - sol.Y2'; 
        
        sol.cond4=[cond4_11, cond4_21', cond4_31';
                   cond4_21, cond4_22 , cond4_32';
                   cond4_31, cond4_32 , cond4_33];
        sol.eig4=eig(sol.cond4);

        %% Print conditions

        if do_print
            for i=1:m
                disp("Eigenvalue Condition 1 (V>0)"); disp(sol.eig1{i});
            end
            disp("Eigenvalue Condition 2 (△V)"); disp(sol.eig2);
            disp("Eigenvalue Condition 3 (r)"); disp(sol.eig3);
            disp("Eigenvalue Condition 4 (Vhat)"); disp(sol.eig4);
        end

    else
        disp("Sistema infactível. Feas:");
        disp(feas);
    end

end