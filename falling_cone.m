close all; clear all; clc;
% building a flat bouncing ball that is rigid
NDIM = 3;

% set angle for cone edges
THETA = 60; % in degrees
THETA_PERT = 0.0;

% initialize arrays
num_nodes = 3;
num_links = 3;
jj = zeros(num_links,1);
kk = zeros(num_links,1);
X = zeros(num_nodes, NDIM);
U = zeros(num_nodes, NDIM);

% set physical parameters
M = ones(num_nodes,1);
M(1:num_nodes) = 5e2*ones(num_nodes,1); % mass of each node
G = 5e3; % magnitude of the gravitational force
Sd = 1e6; % drag force constant

% set numerical parameters
dt = 1e-4;
end_time = 2e-1;
timevec = 0:dt:end_time;

% set positions of nodes around the circumference of the ball
X(1,:) = [0,0,0];
X(2,:) = [cos(THETA*pi/180),0,sin(THETA*pi/180)];
X(3,:) = [-cos(THETA*pi/180),0,sin(THETA*pi/180)];

% perturb initial configuration
PERT = [cos(THETA_PERT*pi/180) 0 -sin(THETA_PERT*pi/180);
        0                      0 0
	sin(THETA_PERT*pi/180) 0 cos(THETA_PERT*pi/180);];
X = (PERT*X')';

% build structure by creating links
jj(1) = 1; kk(1) = 2;
jj(2) = 2; kk(2) = 3;
jj(3) = 3; kk(3) = 1;

% plot initial structure
figure(1);
x = [X(jj,1) X(kk,1)];
y = [X(jj,2) X(kk,2)];
z = [X(jj,3) X(kk,3)];
plot3(x',y',z','linewidth',4)
axis equal

% make a gravitational force
F_gravity = zeros(num_nodes, NDIM);
F_gravity(:,1) = 0;
F_gravity(:,2) = 0;
F_gravity(:,3) = -G.*M; 

% initialize the position of the center of mass
Xcm = (sum((M.*X))./sum(M))';

% initialize velocity center of mass
Ucm = (sum((M.*U))./sum(M))';

% initialize Xtwiddle
Xtwiddle = zeros(num_nodes,NDIM);
Xtwiddle = X - Xcm';

% initialize the angular velocity
L = zeros(NDIM,1);

% here is the timestep loop
centroid_position = zeros(length(timevec), NDIM);
centroid_velocity = zeros(length(timevec), NDIM);
bottom_position = zeros(length(timevec), NDIM);
for t = 1:length(timevec)

    % this just displays the current simulation time in the matlab window
    disp(timevec(t));	

    % compute moment of inertia tensor
    I = zeros(NDIM, NDIM);
    for k = 1:num_nodes
	I = I + M(k).*( (norm(Xtwiddle(k,:))^2).*eye(NDIM) - Xtwiddle(k,:)'*Xtwiddle(k,:) );
    end     

    % compute the angular velocity
    Omega = I\L;

    % normalize angular velocity if it is nonzero
    if(norm(Omega) > 100*eps)
         unit_Omega = Omega/norm(Omega);
         Omega_cross = [0 -Omega(3) Omega(2); Omega(3) 0 -Omega(1); -Omega(2) Omega(1) 0];
         P_Omega = unit_Omega*unit_Omega';
	 Xtwiddle = ( P_Omega*(Xtwiddle') + cos(norm(Omega)*dt).*(eye(NDIM) - P_Omega)*(Xtwiddle') + sin(norm(Omega)*dt).*(Omega_cross*(Xtwiddle'))./norm(Omega) )';
    end

    % compute normal vectors to triangle
    V = X(kk,:) - X(jj,:);
    unit_V = zeros(size(V));
    N = zeros(size(V));    
    for nn = 1:3
       unit_V(nn,:) = V(nn,:)./norm(V(nn,:));
       normal = [unit_V(nn,3) -unit_V(nn,1); unit_V(nn,1) unit_V(nn,3)] \ [1; 0];
       N(nn,1) = normal(1);
       N(nn,3) = normal(2);
    end

    % compute net force and net torque
    net_force = zeros(NDIM,1);
    net_torque = zeros(NDIM,1);
    drag_force = F_drag(X,V,N,kk,jj,Sd,Ucm);
    for l = 1:num_nodes
    	net_force = net_force + F_gravity(l,:)' + drag_force(l,:)';
        net_torque = net_torque + cross(Xtwiddle(l,:)', F_gravity(l,:)' + drag_force(l,:)');
    end	

    % update the position and velocity for the center of mass
    Ucm = Ucm + (dt/sum(M)).*net_force;
    Xcm = Xcm + dt.*Ucm;

    % update the angular momentum
    L = L + dt.*net_torque;    

    % update positions of individual masses
    X = Xtwiddle + Xcm';

    % store the position of centroid
    bottom_position(t,:) = X(1,:);	
    centroid_position(t,:) = Xcm';
    centroid_velocity(t,:) = Ucm';

    % plot the current position of the ball
    figure(2);
    x = [X(jj,1) X(kk,1)];
    y = [X(jj,2) X(kk,2)];
    z = [X(jj,3) X(kk,3)];
    plot3(x',y',z','linewidth',4)
    view(180,0)
    xlim([-0.7 0.7]);
    zlim([-10 1]);	
    pause(0.0001)

    % stop simulation if it blows up.
    if(norm(X)> 1e10)
      break
    end

end

figure(3); hold on
plot(timevec, centroid_position,'linewidth',2)
legend('x','y','z')

figure(4); hold on
plot(timevec, centroid_velocity,'linewidth',2)
legend('x','y','z')

figure(5); hold on
plot(timevec, bottom_position,'linewidth',2)
legend('x','y','z')

% function which defines a drag force
function Fd = F_drag(X,V,N,kk,jj,Sd,Ucm)
    nodes = [1:3];
    Fd = zeros(size(X));
    for l = 1:3 % loop over triangle sides
    	if(N(l,:)*Ucm > 0)
            orthographic = norm( (eye(3) - (Ucm*Ucm'./(norm(Ucm)^2)) )*(V(l,:)') );
            Fd(kk(l),:) = Fd(kk(l),:) - Sd*orthographic.*Ucm';
            Fd(jj(l),:) = Fd(jj(l),:) - Sd*orthographic.*Ucm';
	end
    end
end

