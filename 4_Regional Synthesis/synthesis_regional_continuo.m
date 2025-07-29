clc
clear all

%% 
alphas=[];
for alpha=0.1:0.1:20
dt = 0.01;
x0=[-2.3; 4];
% alpha=1.6;
pause(0.1)

A=[0 1 ; 1 0];
B=[0;-5];

%B=nxm
[n,m] = size(B);
Im = eye(m);
In = eye(n);

%%

setlmis([]) %Initialization

%type = 1,P square and symmetrical
%     = 2,P rectangular [m,n]
%     = 3,P is of other type
%tipo,dim[l,c]
R =lmivar(1,[n,1]);         %quadrado e simétrico

Q0=lmivar(1,[n,1]);         %quadrado e simétrico
Q1=lmivar(2,[n,m]);         %retangular
Q2=lmivar(1,[m,1]);         %quadrado e simetrico 

S =lmivar(1,[m,1]);         %quadrado e simétrico

M =lmivar(2,[n,n]);         %quadrado 

Y1=lmivar(2,[m,n]);         %quadrado e simétrico
Y2=lmivar(2,[m,m]);         %quadrado

Phat=lmivar(1,[n,1]);        %quadrado e simétrico

Z1 = zeros(m, 1);
Z2 = zeros(m, 1);

for i=1:m
    Z1(i)=lmivar(2, [1 n]);
    Z2(i)=lmivar(2, [1 m]);
end

%Condição 1 (V>0)
% |Q0     Q1+Z1'-Y1'   Z1i'| >0
% | *  Q2+He(Z2-Y2+S)  Z2i'|
% | *        *          1  |
for i = 1:m
    idx =-9-i;
    lmiterm([idx 1 1 Q0], 1,1);  
    
    lmiterm([idx 1 2 Q1], 1,1);
    lmiterm([idx 2 1 Z1], 1,1);
    lmiterm([idx 2 1 Y1],-1,1);
    
    lmiterm([idx 3 1 Z1(i)], 1,1);
    
    lmiterm([idx 2 2 Q2], 1,1);  
    lmiterm([idx 2 2 Z2], 1,1,"s"); 
    lmiterm([idx 2 2 Y2],-1,1,"s"); 
    lmiterm([idx 2 2 S ], 1,1,"s"); 
    
    lmiterm([idx 3 2 Z2(i)], 1,1);
    
    lmiterm([idx 3 3 0 ], 1);
end

%Condição 2 (V.<0)
% |He(AM+BY1)   BY2-BS+Y1'   Q0-M+(AM+BY1)'      Q1    |
% |    *        2(Y2-S)     Q1'+Y1+(BY2-BS)'   Q2+Y2-S |<0
% |    *          *             -2M              Y1'   |
% |    *          *              *             2(Y2-S) |
lmiterm([2 1 1 M ],A,1,"s");
lmiterm([2 1 1 Y1],B,1,"s");

lmiterm([2 1 2 Y2],B, 1);
lmiterm([2 1 2 S ],B,-1);
lmiterm([2 2 1 Y1],1, 1);
lmiterm([2 2 1 Z1],1,-1);

lmiterm([2 1 3 Q0],1, 1);  
lmiterm([2 1 3 M ],1,-1);
lmiterm([2 3 1 M ],A, 1);
lmiterm([2 3 1 Y1],B, 1);

lmiterm([2 1 4 Q1],1,1);

lmiterm([2 2 2 Y2], 1,1,"s");
lmiterm([2 2 2 S ],-1,1,"s");
lmiterm([2 2 2 Z2],-1,1,"s");

lmiterm([2 3 2 Q1],1, 1);
lmiterm([2 2 3 Y1],1, 1);
lmiterm([2 3 2 Y2],B, 1);
lmiterm([2 3 2 S ],B,-1);

lmiterm([2 2 4 Q2],1, 1);
lmiterm([2 2 4 Y2],1, 1);
lmiterm([2 2 4 S ],1,-1);

lmiterm([2 3 3 M ],1,-1,"s");

lmiterm([2 4 3 Y1],1, 1);

lmiterm([2 4 4 Y2], 1,1,"s");
lmiterm([2 4 4 S ],-1,1,"s");


%Condition 3 (Aumentar velocidade resposta linear) 
% |He(alphaM+AM+BY1)    -M   |<0
% |alphaM+AM+BY1+R     He(-M)|
lmiterm([3 1 1 M ],alpha*In+A,1,"s");
lmiterm([3 1 1 Y1],B,1,"s");

lmiterm([3 2 1 M ],alpha*In+A,1);
lmiterm([3 2 1 Y1], B,1);
lmiterm([3 2 1 R ], 1,1);
lmiterm([3 1 2 M ],-1,1);               %transpose

lmiterm([3 2 2 M ],-1,1,"s");

% Condição 4 (maximiza área)
% |M+M'-Q0     In       -Q1-Y1'   |
% |   *       Phat         0      | >0
% |   *         *    -Q2+2S-Y2-Y2'|

lmiterm([-4 1 1  M  ], 1,1,"s");
lmiterm([-4 1 1  Q0 ],-1,1);

lmiterm([-4 1 2  0  ],In);

lmiterm([-4 1 3  Q1 ],-1,1);
lmiterm([-4 3 1  Y1 ],-1,1);

lmiterm([-4 2 2 Phat],1,1);

lmiterm([-4 3 3  Q2 ],-1,1);
lmiterm([-4 3 3  S  ], 2,1);
lmiterm([-4 3 3  Y2 ],-1,1,"s");


%Q0>0
lmiterm([-5 1 1 Q0],1,1);

%S>0
lmiterm([-6 1 1 S ],1,1);

%R>0
lmiterm([-7 1 1 R ],1,1);

%Phat>0
lmiterm([-8 1 1 Phat],1,1);


%% declara sistema
LMIsys=getlmis;
feas = feasp(LMIsys);

if round(feas,9) < 0

nx = decnbr(LMIsys); 
c = zeros(nx,1);

for i = 1:nx 
    Phati = defcx(LMIsys,i,Phat); 
    c(i) = trace(Phati);
end

options=[1e-5 0 0 0 1e-3];
% options=[0 0 0 0 0];
[copt, sol] = mincx(LMIsys, c,options);

%extrai valores da solução
Q0opt=dec2mat(LMIsys,sol,Q0);
Q1opt=dec2mat(LMIsys,sol,Q1);
Q2opt=dec2mat(LMIsys,sol,Q2);
Sopt =dec2mat(LMIsys,sol,S );
Mopt =dec2mat(LMIsys,sol,M );
Y1opt=dec2mat(LMIsys,sol,Y1);
Y2opt=dec2mat(LMIsys,sol,Y2);
Z1opt=dec2mat(LMIsys,sol,Z1);
Z2opt=dec2mat(LMIsys,sol,Z2);
Ropt =dec2mat(LMIsys,sol,R );
Phatopt =dec2mat(LMIsys,sol,Phat);


P0=inv(Mopt')*Q0opt*inv(Mopt);
P1=inv(Mopt')*Q1opt*inv(Sopt);
P2=inv(Sopt )*Q2opt*inv(Sopt);
P=[P0 P1; P1' P2];

K1=Y1opt*inv(Mopt);
K2=Y2opt*inv(Sopt);

H1=Z1opt*inv(Mopt);
H2=Z2opt*inv(Sopt);

eigP=eig(P)

if round(eigP(3),4) ~= 0 && round(eigP(2),4) ~= 0 && round(eigP(1),4) ~= 0
    alphas=[alphas;alpha];
end

%% It is Hurwitz?

eigvalues1=eig(A);
eigvalues2=eig(A+B*K1);
eigvalues3=eig(A+B*K1+alpha*In);
% figure
% 
% for i=1:n
%     plot(real(eigvalues1(i)),imag(eigvalues1(i)),'r.','MarkerSize', 10); hold on;
%     plot(real(eigvalues2(i)),imag(eigvalues2(i)),'b.','MarkerSize', 10); hold on;
%     plot(real(eigvalues3(i)),imag(eigvalues3(i)),'k.','MarkerSize', 10); hold on;
% end
% 
% xline(-alpha, '--m', ['\alpha = ' num2str(alpha)]);
% 
% 
% axis equal;
% legend('eigvalues(A)','eigvalues(A+BK1)','eigvalues(A+BK1+aphaI)');
% xlabel('Real'); ylabel('Imaginary');
% title("Eigvalues");

%% Plot traject
x=x0;
N=10;
traj=[];
dV=0;

for k=1:N/dt
    u=K1*x;
    u_sat=min(1,max(-1,u));
    dz=u-u_sat;

    u=K1*x+K2*dz;
    u_sat=min(1, max(-1, u));
    dz=u-u_sat;

    V=[x;dz]'*P*[x;dz]; 

    data=[x;V;dV];
    traj=[traj data];
    
    if k~=1
        dV=(traj(n+1,k)-traj(n+1,k-1))/dt;
    end

    dx=A*x+B*u_sat;
    x=x+dt*dx;
end

[x1_grid,x2_grid]=meshgrid(linspace(-20,20,700),linspace(-20,20,700));

Vgrid=zeros(size(x1_grid));
Wvalues=zeros(size(x1_grid));
H_values=zeros(size(x1_grid));
Vhat_values=zeros(size(x1_grid));

for i=1:size(x1_grid,1)
    for j=1:size(x1_grid,2)

        x=[x1_grid(i,j); x2_grid(i,j)];

        u=K1*x;
        u_sat=min(1,max(-1,u));
        dz=u-u_sat;

        u=K1*x+K2*dz;
        u_sat=min(1,max(-1,u));
        dz=u-u_sat;

        Vgrid(i,j)=[x; dz]'*P*[x; dz];
        Vhat_values(i,j)=x'*Phatopt*x;

        hx=H1*x+H2*dz;
        H_values(i,j) = norm(hx, inf); 
        
        if norm(hx,inf)<=1
            Wvalues(i,j)=min([x; dz]'*P*[x; dz], 1);
        else
            Wvalues(i,j)=1;
        end
    end
end

% figure
% subplot 211
%     plot(0:dt:N-dt,traj(n+1,:));
%     legend("V")
%     titulo=sprintf('Continu');
%     title(titulo)
% 
% subplot 212
%     plot(2*dt:dt:N,traj(n+2, 2:end))
%     legend("dV")

figure
    plot(traj(1,:)  ,traj(2,:)); hold on
    plot(traj(1,1)  ,traj(2,1),'xr','LineWidth',2); hold on
    contour(x1_grid, x2_grid, H_values, [1 1]); hold on;
    contour(x1_grid, x2_grid, Wvalues, [1 1], 'k'); hold on;
    contour(x1_grid, x2_grid, Vhat_values, [1 1], 'r');
    label=sprintf("Synthesis Regional - Continu - $\alpha$= %.2f",alpha);
    legend("trajectory","start","$S_h$","W(x)","$\hat{V}$",'Interpreter','latex')
    title(label,'Interpreter','latex')
    xlabel("x1"); ylabel("x2")




% figure
%     plot3(traj(1,:),traj(2,:),traj(3,:),'r','LineWidth',2); hold on;
%     plot3(traj(1,1),traj(2,1),traj(3,1),'xr','LineWidth',2); hold on;
%     mesh(x1_grid,x2_grid,Vgrid)
%     xlabel("x1"); ylabel("x2"); zlabel("V")
%     legend("trajet states","start","V(x)")


% figure
% for i=1:n
%     subplot(n,1,i)
%     plot(0:dt:N-dt,traj(i,:))
%     xlabel("time")
%     titulo=sprintf('x%d ',i);
%     title(titulo)
% end

end
end
