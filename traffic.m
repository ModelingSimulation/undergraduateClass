%traffic
%meters/seconds units
global d_min d_max v_max
l_road=pi*20; %length of the road
l_car=2; %length of the car
d_min=l_car/2; %minimum distance
d_max=5*l_car; %maximum distance
car_max=floor(l_road/(d_min+l_car)); %maximum number of cars
density_car=1/2; %car density we choose to put
num_car=floor(car_max*density_car); %number of cars based on the density
v_max=30; %maximum velocity
X=zeros(num_car,1);
for k = 1:num_car
    X(k,1)=(k-1)*(l_car+d_min);
end
% plot (X,0,'o');
% axis ([0 l_road -1 1])
% drawnow
dt=0.0005;
t=0;
t_max=10;
clockmax=t_max/dt;
v_car=zeros(num_car,1);
for l = 1:num_car
    v_car(l)=randi(v_max); %choose a random velocity for each car
end
figure (1)
Y=zeros(num_car,2); %to plot in a circle

for i = [1:clockmax]
    X=mod(X,l_road); %periodic domain
    d_car=X([2:num_car,1])-X;
    d_car=mod(d_car,l_road);
     for j = 1:num_car
         v_car(j,1)=vel_car(d_car(j,1));
     end
     
    %accl_car=a_car(d_car,v_car,num_car);
    %v_car=v_car+dt*accl_car;
    
    X=X+dt*v_car; %forward euler
%     plot (X,0,'o'); %plot a line
%     axis ([0 l_road -1 1])
%     drawnow

    %you may want to figure out how below works by yourself
    for m = 1:num_car
        arg=(X(m,1)/l_road)*2*pi;
        Y(m,1)=10*cos(arg);
        Y(m,2)=10*sin(arg);
    end
    plot (Y(:,1),Y(:,2),'o'); %plot a circle
    axis equal
    axis ([-15 15 -15 15])
    drawnow
end

