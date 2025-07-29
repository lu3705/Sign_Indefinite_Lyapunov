clc
clear all

%% 
dt = 0.01;
x0 = [4; -4];
alpha=2;

A=[-0.05 1 ; -10 -0.5];
B=[0;1];

%B=nxm
[n,m] = size(B);
Im = eye(m);
In = eye(n);

%%

setlmis([]) 

R =lmivar(1,[n,1]);

Q0=lmivar(1,[n,1]);
Q1=lmivar(2,[n,m]);
Q2=lmivar(1,[m,1]);

S =lmivar(1,[m,1]);

M =lmivar(2,[n,n]);

Y1=lmivar(2,[m,n]);
Y2=lmivar(2,[m,m]);


%Condição 1 (V>0)
% |Q0      Q1-Y1'    | >0
% | *   Q2-He(Y2)+2S |    
lmiterm([-1 1 1 Q0], 1,1);  

lmiterm([-1 1 2 Q1], 1,1);
lmiterm([-1 2 1 Y1],-1,1);      %transpose

lmiterm([-1 2 2 Q2], 1,1);
lmiterm([-1 2 2 Y2],-1,1,"s"); 
lmiterm([-1 2 2 S ], 2,1); 

%Condição 2 (V.<0)
% |He(AM+BY1)   BY2-BS+Y1'-Z1'    Q0-M+(AM+BY1)'        Q1    |
% |    *        He(Y2-Z2-S)      Q1'+Y1+(BY2-BS)'    Q2+Y2-S  |<0
% |    *          *                  He(-M)             Y1'   |
% |    *          *                    *             He(Y2-S) |
lmiterm([2 1 1 M ],A,1,"s");
lmiterm([2 1 1 Y1],B,1,"s");

lmiterm([2 1 2 Y2], B,1);
lmiterm([2 1 2 S ],-1*B,1);
lmiterm([2 2 1 Y1], 1,1);           %transpose

lmiterm([2 1 3 Q0],1,1);  
lmiterm([2 1 3 M ],-1,1);
lmiterm([2 3 1 M ],A,1);           %transpose
lmiterm([2 3 1 Y1],B,1);           %transpose

lmiterm([2 1 4 Q1],1,1);

lmiterm([2 2 2 Y2], 1,1,"s");
lmiterm([2 2 2 S ],-1,1,"s");

lmiterm([2 3 2 Q1],1,1);           %transpose
lmiterm([2 2 3 Y1],1,1);
lmiterm([2 3 2 Y2],B,1);           %transpose
lmiterm([2 3 2 S ],-1*B,1);        %transpose

lmiterm([2 2 4 Q2],1,1);
lmiterm([2 2 4 Y2],1,1);
lmiterm([2 2 4 S ],-1,1);

lmiterm([2 3 3 M ],-1,1,"s");

lmiterm([2 4 3 Y1],1,1);           %transpose

lmiterm([2 4 4 Y2], 1,1,"s");
lmiterm([2 4 4 S ],-1,1,"s");

%Condition 3 (Aumentar velocidade resposta linear) 
% |He(alphaM+AM+BY1)     *   |<0
% |alphaM+AM+BY1+R-M'  He(-M)|
lmiterm([3 1 1 M ],alpha*In+A,1,"s");
lmiterm([3 1 1 Y1],B,1,"s");

lmiterm([3 2 1 M ],alpha*In+A,1);
lmiterm([3 2 1 Y1], B,1);
lmiterm([3 2 1 R ], 1,1);
lmiterm([3 1 2 M ],-1,1);               %transpose

lmiterm([3 2 2 M ],-1,1,"s");

%Q0>0
lmiterm([-4 1 1 Q0],1,1);

%S>0
lmiterm([-5 1 1 S ],1,1);

%R>0
lmiterm([-6 1 1 R ],1,1);

%% declara sistema
LMIsys=getlmis;

options=[1e-5 0 0 0 1e-3];
[feas,sol]=feasp(LMIsys,options);                
%feas<0 (factivel)
%sol conjunto de soluções

%extrai valores da solução
Q0opt=dec2mat(LMIsys,sol,Q0);
Q1opt=dec2mat(LMIsys,sol,Q1);
Q2opt=dec2mat(LMIsys,sol,Q2);
Sopt =dec2mat(LMIsys,sol,S );
Mopt =dec2mat(LMIsys,sol,M );
Y1opt=dec2mat(LMIsys,sol,Y1);
Y2opt=dec2mat(LMIsys,sol,Y2);
Ropt =dec2mat(LMIsys,sol,R );

P0=inv(Mopt')*Q0opt*inv(Mopt);
P1=inv(Mopt')*Q1opt*inv(Sopt);
P2=inv(Sopt )*Q2opt*inv(Sopt);
P=[P0 P1; P1' P2];

K1=Y1opt*inv(Mopt);
K2=Y2opt*inv(Sopt);

eig(P)


%% It is Hurwitz?

eigvalues1=eig(A);
eigvalues2=eig(A+B*K1);
eigvalues3=eig(A+B*K1+alpha*In);
figure

for i=1:n
    plot(real(eigvalues1(i)),imag(eigvalues1(i)),'r.','MarkerSize', 10); hold on;
    plot(real(eigvalues2(i)),imag(eigvalues2(i)),'b.','MarkerSize', 10); hold on;
    plot(real(eigvalues3(i)),imag(eigvalues3(i)),'k.','MarkerSize', 10); hold on;
end

xline(-alpha, '--m', ['\alpha = ' num2str(alpha)]);
xline(0, '--','Color',[.8 .8 .8]);


axis equal;
legend('eigvalues(A)','eigvalues(A+BK1)','eigvalues(A+BK1+aphaI)');
xlabel('Real'); ylabel('Imaginary');
title("Eigvalues");

%% Plot traject
x=x0;
N=10;
traj=[];
dV=0;

for k=1:N/dt
    u=K1*x;
    u_sat=min(1,max(-1,u));
    dz=u-u_sat;

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

[x1_grid,x2_grid]=meshgrid(linspace(-5,5,50),linspace(-12,12,50));

Vgrid=zeros(size(x1_grid));
dz=0;
for i=1:size(x1_grid,1)
    for j=1:size(x1_grid,2)

        x=[x1_grid(i,j); x2_grid(i,j)];

        u=K1*x+K2*dz;           %sempre vai ta usando dz antigo
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
plot(traj(1,1)  ,traj(2,1),'xr','LineWidth',2); hold on
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
