clc
clear all

%% 
dt = 0.01;
x0 = [4; -0.82];
tal=0.1;

%System
%Ex1 (Exemplo 2 - Artigo Sophie)
A=[0 1 ; 1 0];
B=[0;-5];
K1=[2.6 1.4];
K2=0;

% % Ex2 (Exemplo 3 - Artigo Sophie - Instável)   (RETIRAR LMI 2 E 5)
% A=[0 1 ; 0 1];
% B=[-1;1];
% K1=[-0.1 2];
% K2=0.9;

%B=nxmy
[n,m] = size(B);
Im = eye(m);
In = eye(n);

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

Phat=lmivar(1,[n,1]);        %quadrado e simétrico
That=lmivar(1,[m,1]);        %quadrado e simétrico

%Condição 0 (T-tal>0)
lmiterm([-1 1 1 T1], 1 , 1);  
lmiterm([-1 1 1 0 ],-1*tal*Im);  

M1=A+B*K1;
M2=K2-Im;

%Condição 1 (V>0 e inclusão)
% |P0    [P1-(K1-H1)'T]  H1'T| 
% | *   P2-2T(K2-Im-H2)  H2'T| >0
% | *           *       tal*T| 
lmiterm([-2 1 1 P0],1, 1);  

lmiterm([-2 1 2 P1],1, 1);
lmiterm([-2 1 2 T1],-1*K1',1);
lmiterm([-2 1 2 Z1],1,1);           %transpose

lmiterm([-2 1 3 Z1],1,1);           %transpose

lmiterm([-2 2 2 P2],1,1);  
lmiterm([-2 2 2 T1],-2*M2',1);
lmiterm([-2 2 2 Z2],2,1); 

lmiterm([-2 2 3 Z2],1,1);           %transpose

lmiterm([-2 3 3 T1],tal,1);


%Condição 2 (V.<0 somente em Sh)
% |2P0*(A+BK1)      P0B(K2-Im)+(A+B*K1)'P1+(K1-H1)'T1+[K1(A+B*K1)]'T3             P1+[K1(A+BK1)]'T2         | 
% |     *           2[B*(K2-Im)]'*P1+2(K2-Im-H2)'T1+2(K1*B(K2-Im))'T3       P2+[K1*B(K2-Im)]'T2+T3(K2-Im) | <0
% |     *                                *                                           2T2*(K2-Im)            |

lmiterm([3 1 1 P0],2, M1);

lmiterm([3 1 2 P0],1,B*M2);
lmiterm([3 1 2 P1],M1', 1);
lmiterm([3 1 2 T1],K1',1);
lmiterm([3 1 2 Z1],-1,1);            %transpose
lmiterm([3 1 2 T3],(K1*M1)',1);

lmiterm([3 1 3 P1],1,1);
lmiterm([3 1 3 T2],(K1*M1)',1);

lmiterm([3 2 2 P1],2*(B*M2)',1);
lmiterm([3 2 2 T1],2*M2',1);
lmiterm([3 2 2 Z2],-2,1);
lmiterm([3 2 2 T3],2*(K1*B*M2)',1);

lmiterm([3 2 3 P2],1,1);
lmiterm([3 2 3 T2],(K1*B*M2)',1);
lmiterm([3 2 3 T3],1,M2);

lmiterm([3 3 3 T2],2,M2);

% Condição 3 (maximiza área)
% | Phat-P0   -(P1+K1'That)    | >0
% |    *     -(P2+2That(K2-Im))|

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
% lmiterm([-8 1 1 T3],1,1); (NAO PRECISA POIS ENTRA EM dzT3sat=0)
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

%extrai valores da solução
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

eig(P)

%% Plot traject
x=x0;
N=10;
traj=[];
dV=0;
dz=0;

for k=1:N/dt

    u=K1*x+K2*dz;                   %sempre vai ta usando o dz antigo
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

[x1_grid,x2_grid]=meshgrid(linspace(-12,12,700),linspace(-12,12,700));

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

figure
subplot 211
    plot(0:dt:N-dt,traj(n+1,:));
    legend("V")
    titulo=sprintf('Continu');
    title(titulo)

subplot 212
    plot(2*dt:dt:N,traj(n+2, 2:end))
    legend("dV")

figure
plot(traj(1,:)  ,traj(2,:)); hold on
plot(traj(1,1)  ,traj(2,1),'xr','LineWidth',2); hold on
contour(x1_grid, x2_grid, H_values, [1 1]); hold on;
contour(x1_grid, x2_grid, Wvalues, [1 1], 'k'); hold on;
contour(x1_grid, x2_grid, Vhat_values, [1 1], 'r');
label=sprintf("$ S_h : \\tau = %.1f $", tal);
legend("trajectory","start",label,"W(x)","$\hat{V}$",'Interpreter','latex')
xlabel("x1"); ylabel("x2")
axis equal



figure
% plot3(traj(1,:),traj(2,:),traj(3,:),'r','LineWidth',2); hold on;
% plot3(traj(1,1),traj(2,1),traj(3,1),'xr','LineWidth',2); hold on;
mesh(x1_grid,x2_grid,Vgrid)
xlabel("x1"); ylabel("x2"); zlabel("V")
legend("trajet states","start","V(x)")

%%
figure
for i=1:n
    subplot(n,1,i)
    plot(0:dt:N-dt,traj(i,:))
    xlabel("time")
    titulo=sprintf('x%d ',i);
    title(titulo)
end
