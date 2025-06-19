clear all;
clc;

%% 
Ts=0.5;
dt=0.01;
x0=[10;10];

%System continu
Ac=[-0.05 1; -10 -0.5];
Bc=[0; 1];
K1=[9.9 0.495];
K2=0;

[n,m]=size(Bc);
Im=eye(m);
In=eye(n);

Cc=zeros(1,n);
Dc=zeros(1,m);

%caso global
H2=zeros(m,m);
H1=zeros(m,n);

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

%% Conditions 

setlmis([]);

%type=1,square and symmetrical
%    =2,rectangular [m,n]
P0=lmivar(1,[n,1]);
P1=lmivar(2,[n,m]);
P2=lmivar(1,[m,1]);

T =lmivar(1,[m,1]);
T1=lmivar(1,[m,1]);
T2=lmivar(1,[m,1]);
T3=lmivar(1,[m,1]);

M1=A+B*K1;
M2=K2-Im;

%Condition 1 (V > 0) – same for continuous and discrete
%|P0   2[P1-(K1-H1)'T] | >0
%| 0   P2-2T(K2-Im-H2) |    
lmiterm([-1 1 1 P0],1,1);  

lmiterm([-1 1 2 P1],1,1);
lmiterm([-1 1 2 T],-1*(K1-H1)',1);

lmiterm([-1 2 2 P2],1,1);  
lmiterm([-1 2 2 T],-2,(M2-H2)); 

%Condition 2 (△V<0)
%|(A+BK1)'P0(A+BK1)-P0    2(A+BK1)'P0B(K2-Im)-2P1-2[K1(A+BK1-In)]'T3+2(K1-H1)'T1              2(A+BK1)'P1+2[K1(A+BK1-In)]'T2           | 
%|       0             [B(K2-Im)]'P0[B(K2-Im)]-P2-2T3[(K1B-Im)(K2-Im)]+2T1(K2-Im-H2)   2[B(K2-Im)]'P1+2[(K1B-Im)(K2-Im)]'T2-2T3(K2-Im)] | <0
%|       0                                       0                                                     P2+2T2(K2-Im)                   |

lmiterm([2 1 1 P0],M1',M1);  
lmiterm([2 1 1 P0],-1,1); 

lmiterm([2 1 2 P0],M1',B*M2);
lmiterm([2 1 2 P1],-1,1);
lmiterm([2 1 2 T3],-1*(K1*(M1-In))',1);
lmiterm([2 1 2 T1],(K1-H1)',1);

lmiterm([2 1 3 P1],(M1)',1);
lmiterm([2 1 3 T2],(K1*(M1-In))',1);

lmiterm([2 2 2 P0],(B*M2)',B*M2);
lmiterm([2 2 2 P2],-1,1);
lmiterm([2 2 2 T3],-2,(K1*B-Im)*M2);
lmiterm([2 2 2 T1],2,M2-H2);

lmiterm([2 2 3 P1],(B*(M2))',1);
lmiterm([2 2 3 T2],((K1*B-Im)*(M2))',1);
lmiterm([2 2 3 T3],-1,(M2));

lmiterm([2 3 3 P2],1,1);
lmiterm([2 3 3 T2],2,M2);

%P0>0
lmiterm([-3 1 1 P0],1,1);

%Ts>0
lmiterm([-4 1 1 T ],1,1);
lmiterm([-5 1 1 T1],1,1);
lmiterm([-6 1 1 T2],1,1);
lmiterm([-7 1 1 T3],1,1);

%% Solve
LMIsys=getlmis;

[feas,sol]=feasp(LMIsys);      

P0=dec2mat(LMIsys,sol,P0);
P1=dec2mat(LMIsys,sol,P1);
P2=dec2mat(LMIsys,sol,P2);

P=[P0 P1; P1' P2];

%% Plot traject
x=x0;
N=10;
traj=[];
V_values=[];
dV_values=[];

for k=1:N/Ts
    traj=[traj x];

    u=K1*x;
    u_sat=min(1,max(-1,u));
    dz=u-u_sat;

    u=K1*x+K2*dz;
    dz=u+min(1,max(-1,u));

    V=[x; dz]'*P*[x; dz];
    V_values=[V_values V];

    if k~=1
        dV=V_values(k)-V_values(k-1);
        dV_values=[dV_values dV];
    end

    x=A*x+B*u_sat;

end

x_c =x0;
traj_c=[];

for k=1:N/dt
    traj_c=[traj_c x_c];

    u_c=K1*x_c;
    u_sat_c=min(1,max(-1,u_c));
    dz_c=u_c-u_sat_c;

    u_c=K1*x_c+K2*dz_c;
    u_sat_c=min(1,max(-1,u_c));

    dx_c =Ac*x_c+Bc*u_sat_c;
    x_c=x_c+dt*dx_c;
end

figure
subplot 311
    plot(0:Ts:N-Ts,V_values,'.');
    legend("V")
    titulo=sprintf('Discret-Desconsidera sat no continuo-Ts=%.2f',Ts);
    title(titulo)

subplot 312
    plot(2*Ts:Ts:N,dV_values,'.')
    legend("△V")

subplot 313
    plot(traj(1,:)  ,traj(2,:),'.'); hold on
    plot(traj_c(1,:),traj_c(2,:))  ; hold on
    plot(traj(1,1)  ,traj(2,1),'xr','LineWidth',2)
    legend("discret","continu","start")
    xlabel("x1"); ylabel("x2")


%%
[x1_grid,x2_grid]=meshgrid(linspace(-15,15,30/Ts),linspace(-35,30,65/Ts));

Vgrid=zeros(size(x1_grid));
for i=1:size(x1_grid,1)
    for j=1:size(x1_grid,2)
        x=[x1_grid(i,j); x2_grid(i,j)];

        u=K1*x;
        u_sat=min(1,max(-1,u));
        dz=u-u_sat;

        u=K1*x+K2*dz;
        dz=u+min(1,max(-1,u));

        Vgrid(i,j)=[x; dz]'*P*[x; dz];
    end
end

figure
plot3(traj(1,:),traj(2,:),V_values,'.r','MarkerSize',10)
xlabel("x1"); ylabel("x2"); zlabel("V")
hold on
mesh(x1_grid,x2_grid,Vgrid)
legend("trajet states","V")

%%
figure
for i=1:n
    subplot(n,1,i)
    plot(0:Ts:N-Ts,traj(i,:),'.')
    hold on;
    plot(0:dt:N-dt,traj_c(i,:))
    legend("discret","continu")
    xlabel("time")
    titulo=sprintf('x%d ',i);
    title(titulo)
end