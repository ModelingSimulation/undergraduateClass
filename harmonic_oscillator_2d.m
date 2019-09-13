%9/12 Recitation
clear all
%%Initialization

X=zeros(2,2);
X(1,:)=[0 1];
X(2,:)=[0 0];
U=zeros(2,2);
mass=[1;1];
links=[1 2];
dt=0.01;
T=10;
clockmax=T/dt;

for timer=1:clockmax
    
    t_sim=timer*dt
%%  Forward Euler Part  
    f=force_spr(X,links,0.8,5); %Force Calculation
    a=f./mass; %Acceleration  
    U=U+dt*a; %Update Time
    X=X+U*dt; %Update Position
%%  Plot   
    figure (1)
    subplot(1,2,1)
    plot (X(:,1),X(:,2),'.','MarkerSize',50,'Color',[0.6350 0.0780 0.1840]);
    name_title=sprintf('time=%g',t_sim);
    title (name_title);
    xlim([-1 1])
    ylim([-1 2])
    axis equal
    
    subplot(1,2,2)
    plot (t_sim,X(1,2),'b-o');
    hold on
    plot (t_sim,X(2,2),'r-*');
    xlim([0 5])
    ylim([-1 2])
end