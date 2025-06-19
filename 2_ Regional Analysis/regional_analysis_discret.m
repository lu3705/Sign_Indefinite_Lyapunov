clc
clear all

%% 
Ts=0.25;
dt = 0.01;
x0=[-2;4];
tal=0.01;

% System discret 
% Ex1(Pag 118 livro)
% A=[0.8 0.5;-0.4 1.2];
% B=[0;1];                    %Problema: limite saturação [-7,7]
% K1=[0.2888 -1.835];
% K2=0;

% % %Ex2
% A=[1.0048 0.0953;0.0953 0.9095];
% B=[0.1050;0.1002];                   
% K1=[-1.4551 -0.8925];
% K2=0;

%Ex3 (Discretização do Exemplo 2 de [0]Sophie)
Ac=[0 1 ; 1 0];
Bc=[0;-5];
K1=[2.6 1.4];
K2=0;


%B=nxmy
[n,m] = size(Bc);
Im = eye(m);
In = eye(n);

Cc=zeros(1,n);
Dc=zeros(1,m);

%% Discretizar

sysc=ss(Ac,Bc,Cc,Dc);                   

sysd=c2d(sysc,Ts);

A=sysd.A;
B=sysd.B;

%% Check if is Schur Cohm

angles=linspace(0,2*pi,100);
x = cos(angles);
y = sin(angles);

eigvalues1=eig(A);
eigvalues2=eig(A+B*K1);

figure
for i=1:n
    plot(real(eigvalues1(i)),imag(eigvalues1(i)),'r.','MarkerSize', 10); hold on;
    plot(real(eigvalues2(i)),imag(eigvalues2(i)),'b.','MarkerSize', 10); hold on;
end

plot(x, y, 'Color',[.8 .8 .8]); hold on
axis equal;
legend('eigvalues(A)','eigvalues(A+BK1)');
xlabel('Real'); ylabel('Imaginary');
title("Eigvalues - Schur Cohm");


%%

setlmis([]) %Initialization


%type = 1,square and symmetrical
%     = 2,rectangular [m,n]
P0=lmivar(1,[n,1]);          %quadrado e simétrico
P1=lmivar(2,[n,m]);          %retangular
P2=lmivar(1,[m,1]);          %quadrado e simetrico 

T1=lmivar(1,[m,1]);          %quadrado e simétrico
T2=lmivar(1,[m,1]);          %quadrado e simétrico
T3=lmivar(1,[m,1]);          %quadrado e simétrico

Z1=lmivar(2,[n,m]);          %quadrado e simétrico  (ESTOU DECLARANDO Z1=(TH1)'
Z2=lmivar(1,[m,1]);          %quadrado e simétrico

Phat=lmivar(1,[n,1]);         %quadrado e simétrico
That=lmivar(1,[m,1]);        %quadrado e simétrico

%Condição 0 (T-tal>0)
lmiterm([-1 1 1 T1], 1 , 1);  
lmiterm([-1 1 1 0 ],-1*tal*Im);  

M1=A+B*K1;
M2=K2-Im;

%Condição 1 (V>0 e inclusão)
% |P0   2[P1-(K1-H1)'T]  2H1'T| 
% | 0   P2-2T(K2-Im-H2)  2H2'T| >0
% | 0           0        tal*T| 
lmiterm([-2 1 1 P0],1, 1);  

lmiterm([-2 1 2 P1],1, 1);
lmiterm([-2 1 2 T1],-1*K1',1);
lmiterm([-2 1 2 Z1],1,1);           %transpose

lmiterm([-2 1 3 Z1],1,1);           %transpose

lmiterm([-2 2 2 P2],1,1);  
lmiterm([-2 2 2 Z2],2,1); 
lmiterm([-2 2 2 T1],-2*M2',1);

lmiterm([-2 2 3 Z2],1,1);           %transpose

lmiterm([-2 3 3 T1],tal,1);


%Condition 2 (△V<0)
%|(A+BK1)'P0(A+BK1)-P0    (A+BK1)'P0B(K2-Im)-P1-[K1(A+BK1-In)]'T3+(K1-H1)'T1              (A+BK1)'P1+[K1(A+BK1-In)]'T2           | 
%|       *             [B(K2-Im)]'P0[B(K2-Im)]-P2-2T3[(K1B-Im)(K2-Im)]+2T1(K2-Im-H2   [B(K2-Im)]'P1+[(K1B-Im)(K2-Im)]'T2-T3(K2-Im)] | <0
%|       *                                       *                                                     P2+2T2(K2-Im)                   |

lmiterm([3 1 1 P0],M1',M1);  
lmiterm([3 1 1 P0],-1,1); 

lmiterm([3 1 2 P0],M1',B*M2);
lmiterm([3 1 2 P1],-1,1);
lmiterm([3 1 2 T3],-1*(K1*(M1-In))',1);
lmiterm([3 1 2 T1],K1',1);
lmiterm([3 1 2 Z1],1,1);

lmiterm([3 1 3 P1],M1',1);
lmiterm([3 1 3 T2],(K1*(M1-In))',1);

lmiterm([3 2 2 P0],(B*M2)',B*M2);
lmiterm([3 2 2 P2],-1,1);
lmiterm([3 2 2 T3],-2,(K1*B-Im)*M2);
lmiterm([3 2 2 T1],2*M2',1);
lmiterm([3 2 2 Z2],-2,1);

lmiterm([3 2 3 P1],(B*M2)',1);
lmiterm([3 2 3 T2],((K1*B-Im)*M2)',1);
lmiterm([3 2 3 T3],-1,M2);

lmiterm([3 3 3 P2],1,1);
lmiterm([3 3 3 T2],2,M2);

% Condição 3 (maximiza área)
% | Pchap-P0   -2(P1+K1'Tchap)    | >0
% | 0        -(P2+2Tchap(K2-Im))  |

lmiterm([-4 1 1 Phat],1,1);
lmiterm([-4 1 1 P0],-1,1);

lmiterm([-4 1 2 P1],-1,1);
lmiterm([-4 1 2 That],-1*K1',1);

lmiterm([-4 2 2 P2],-1,1);
lmiterm([-4 2 2 That],-2,M2);


% P0>0
lmiterm([-5 1 1 P0],1,1);

%Phat>0
lmiterm([-6 1 1 Phat],1,1);

%T>0
lmiterm([-7 1 1 T2],1,1);
lmiterm([-8 1 1 T3],1,1); 
lmiterm([-9 1 1 That],1,1);


%% Solve 

LMIsys=getlmis;
feas = feasp(LMIsys);

nx = decnbr(LMIsys); 
c = zeros(nx,1);

for i = 1:nx 
    [Phati] = defcx(LMIsys,i,Phat); 
    c(i) = trace(Phati);
end

options=[1e-5 0 0 0 1e-3];

% c = mat2dec(LMIsys, Phat);
[copt, sol] = mincx(LMIsys, c,options);

P0opt=dec2mat(LMIsys,sol,P0);
P1opt=dec2mat(LMIsys,sol,P1);
P2opt=dec2mat(LMIsys,sol,P2);

T1opt=dec2mat(LMIsys,sol,T1);
T2opt=dec2mat(LMIsys,sol,T2);
T3opt=dec2mat(LMIsys,sol,T3);

Z1opt=dec2mat(LMIsys,sol,Z1);
Z2opt=dec2mat(LMIsys,sol,Z2);

Phatopt=dec2mat(LMIsys,sol,Phat);
Thatopt=dec2mat(LMIsys,sol,That);

H1=inv(T1opt)*Z1opt';
H2=inv(T1opt)*Z2opt';

P=[P0opt P1opt ; P1opt' P2opt];


%% Plot traject
x=x0;
N=5;
traj=[];
dV=0;
dz=0;

for k=1:N/Ts

    u = K1*x + K2*dz;                   %sempre vai ta usando o dz antigo
    u_sat = min(1, max(-1, u));
    dz = u - u_sat;

    V=[x; dz]'*P*[x; dz]; 

    data=[x; V; dV];
    traj=[traj data];
    
    if k~=1
        dV=(traj(n+1,k)-traj(n+1,k-1));
    end

    x=A*x+B*u_sat;
end


[x1_grid,x2_grid]=meshgrid(linspace(-20,20,40/Ts),linspace(-20,20,40/Ts));

Vgrid=zeros(size(x1_grid));
Wvalues=zeros(size(x1_grid));
H_values=zeros(size(x1_grid));
Vhat_values=zeros(size(x1_grid));
dz=0;
for i=1:size(x1_grid,1)
    for j=1:size(x1_grid,2)
        x=[x1_grid(i,j); x2_grid(i,j)];

        u=K1*x+K2*dz;           %sempre vai ta usando dz antigo
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
%     plot(0:Ts:N-Ts,traj(n+1,:),'.');
%     legend("V")
%     titulo=sprintf('Discret-Desconsidera sat no continuo-Ts=%.2f',Ts);
%     title(titulo)
% subplot 212
%     plot(2*Ts:Ts:N,traj(n+2, 2:end),'.')
%     legend("△V")

figure
    plot(traj(1,:)  ,traj(2,:),'.'); hold on
    plot(traj(1,1)  ,traj(2,1),'xr','LineWidth',2); hold on
    contour(x1_grid, x2_grid, H_values, [1 1]); hold on;
    contour(x1_grid, x2_grid, Wvalues, [1 1], 'k'); hold on;
    contour(x1_grid, x2_grid, Vhat_values, [1 1], 'r');
    label=sprintf("$ S_h : \\tau = %.2f $", tal);
    legend("trajectory","start",label,"W(x)","$\hat{V}$",'Interpreter','latex')
    xlabel("x1"); ylabel("x2")
    label2=sprintf("$ T_s = %.2f$", Ts);
    title(label2,'Interpreter','latex');

% figure
%     plot3(traj(1,:),traj(2,:),traj(3,:),'.r','MarkerSize',10)
%     xlabel("x1"); ylabel("x2"); zlabel("V")
%     hold on
%     mesh(x1_grid,x2_grid,Vgrid)
%     legend("trajet states","V")
% 
% figure
% for i=1:n
%     subplot(n,1,i)
%     plot(0:Ts:N-Ts,traj(i,:),'.')
%     legend("discret","continu")
%     xlabel("time")
%     titulo=sprintf('x%d ',i);
%     title(titulo)
% end