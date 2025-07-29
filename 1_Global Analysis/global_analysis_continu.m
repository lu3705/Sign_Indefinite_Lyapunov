clc
clear all

%%
dt = 0.01;
x0 = [4; -0.82];

%System
%Exemplo 1 - Artigo Sophie
A=[-0.05 1 ; -10 -0.5];
B=[0;1];
% K1=[9.9 0.495];
% K2=0;

K1=[-0.4322   -4.0410];
K2= -0.2861;

%B=nxm
[n,m]=size(B);
Im=eye(m);
In=eye(n);

H2=zeros(m,m);
H1=zeros(m,n);

%%

setlmis([]) %Initialization

%type = 1,P square and symmetrical
%     = 2,P rectangular [m,n]
%     = 3,P is of other type
%tipo,dim[l,c]
P0=lmivar(1,[n,1]);          %quadrado e simétrico
P1=lmivar(2,[n,m]);          %retangular
P2=lmivar(1,[m,1]);          %quadrado e simetrico 

T =lmivar(1,[m,1]);          %quadrado e simétrico
T1=lmivar(1,[m,1]);          %quadrado e simétrico
T2=lmivar(1,[m,1]);          %quadrado e simétrico
T3=lmivar(1,[m,1]);          %quadrado e simétrico

M1=A+B*K1;
M2=K2-Im;

%Condição 1 (V>0)
% |P0   [P1-(K1-H1)'T] | >0
% | *   P2-2T(K2-Im-H2) |    
lmiterm([-1 1 1 P0],1,1);  

lmiterm([-1 1 2 P1],1,1);
lmiterm([-1 1 2 T ],-1*(K1-H1)',1);

lmiterm([-1 2 2 P2],1,1);  
lmiterm([-1 2 2 T ],-2,M2-H2); 

%Condição 2 (V.<0)
% |P0*(A+BK1)  (A+B*K1)'*P1+P0*B*(K2-Im)+[K1(A+B*K1)]'*T3+(K1-H1)'*T1      P1+[K1(A+BK1)]'T2           | 
% | 0               T3*K1*B(K2-Im)+T1*(K2-Im-H2)+[B*(K2-Im)]'*P1       T3*(K2-Im)+P2+[K1*B(K2-Im)]'*T2 | <0
% | 0                             0                                           T2*(K2-Im)               |
lmiterm([2 1 1 P0],1, M1);  

lmiterm([2 1 2 P1],M1', 1);
lmiterm([2 1 2 P0],1,B*M2);
lmiterm([2 1 2 T3],(K1*M1)',1);
lmiterm([2 1 2 T1],(K1-H1)',1);

lmiterm([2 2 2 P1],(B*M2)',1); 
lmiterm([2 2 2 T3],1,K1*B*M2);  
lmiterm([2 2 2 T1],1,M2-H2); 

lmiterm([2 1 3 P1],1,1);  
lmiterm([2 1 3 T2],(K1*M1)',1); 

lmiterm([2 2 3 T3],1,M2);  
lmiterm([2 2 3 P2],1,1);
lmiterm([2 2 3 T2],(K1*B*M2)',1); 

lmiterm([2 3 3 T2],1,M2);  

%P0>0
lmiterm([-3 1 1 P0],1,1);

%T>0
lmiterm([-4 1 1 T ],1,1);
lmiterm([-5 1 1 T1],1,1);
lmiterm([-6 1 1 T2],1,1);
% lmiterm([-7 1 1 T3],1,1);   NAO PRECISA 


%% Solve 

LMIsys=getlmis;
options=[1e-5 0 0 0 1e-3];
[feas,sol]=feasp(LMIsys, options); 
%feas<0 (factivel)
%sol conjunto de soluções

%extrai valores da solução
P0opt=dec2mat(LMIsys,sol,P0);
P1opt=dec2mat(LMIsys,sol,P1);
P2opt=dec2mat(LMIsys,sol,P2);

Topt =dec2mat(LMIsys,sol,T );
T1opt=dec2mat(LMIsys,sol,T1);
T2opt=dec2mat(LMIsys,sol,T2);
T3opt=dec2mat(LMIsys,sol,T3);

P=[P0opt P1opt; P1opt' P2opt];

eig(P)

%% It is Hurwitz?

eigvalues1=eig(A);
eigvalues2=eig(A+B*K1);

figure
for i=1:n
    plot(real(eigvalues1(i)),imag(eigvalues1(i)),'r.','MarkerSize', 10); hold on;
    plot(real(eigvalues2(i)),imag(eigvalues2(i)),'b.','MarkerSize', 10); hold on;
end


xline(0, '--','Color',[.8 .8 .8]);
axis equal;
legend('eigvalues(A)','eigvalues(A+BK1)');
xlabel('Real'); ylabel('Imaginary');
title("Eigvalues");


%% Plot traject
x=x0;
N = 10;
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

[x1_grid,x2_grid]=meshgrid(linspace(-5,5,50),linspace(-12,12,50));

Vgrid=zeros(size(x1_grid));

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
plot(traj(1,1)  ,traj(2,1),'xr','LineWidth',2);
legend("trajectory","start")
xlabel("x1"); ylabel("x2")
axis equal


figure
plot3(traj(1,:),traj(2,:),traj(3,:),'r','LineWidth',2); hold on;
plot3(traj(1,1),traj(2,1),traj(3,1),'xr','LineWidth',2); hold on;
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