function [accel] = a_car(distance,velocity,num_car)
global v_max
accel=zeros(num_car,1);
for k = 1:num_car-1
    %"Follow the leader"
    accel(k)=(velocity(k+1)-velocity(k))/distance(k);
    %"Optimal Velocity"
    %accel(k)=v_max*distance(k)-velocity(k);
end
accel(num_car)=(velocity(1)-velocity(num_car))/distance(num_car);
end

