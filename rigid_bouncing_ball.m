clear all; clc;
% building a bouncing ball that is rigid
NDIM = 2;

% number of nodes around the circumference of the ball
%n = 30;
n = 10; 

% initialize arrays
num_nodes = n+1;
num_links = 2*n;
jj = zeros(num_links,1);
kk = zeros(num_links,1);
X = zeros(num_nodes, NDIM);
U = zeros(num_nodes, NDIM);

% set physical parameters
RB = 0.8; % radius of ball
center = [1, 10.5]; % initial location for ball center
M = 5e-7*ones(num_nodes,1); % mass of each node 
G = 50000; % magnitude of the gravitational force
Sg = 500000; % strength of the force exerted by the ground

% set numerical parameters
dt = 1e-4;
end_time = 1e0;
timevec = 0:dt:end_time;

% set positions of nodes around the circumference of the ball
for k = 1:n
    theta = 2*pi*k/n;
    X(k,:) = center + RB*[cos(theta), sin(theta)];
end
% set position of center
X(n+1,:) = center;

% naming index sets for the links
spokes = 1:n;
rimlinks = (n+1):2*n;

% build structure by creating links
jj(spokes) = 1:n;
kk(spokes) = n+1;
jj(rimlinks) = 1:n;
kk(rimlinks) = [2:n 1];

% plot initial structure
figure(1);
x = [X(jj,1) X(kk,1)];
y = [X(jj,2) X(kk,2)];
plot(x',y','linewidth',4)

% initialize rest length
DX = X(jj,:) - X(kk,:);
Rzero = sqrt(sum(DX.^2,2));

% make a gravitational force
F_gravity = zeros(num_nodes, NDIM);
F_gravity(:,1) = 0;
F_gravity(:,2) = -G; 

% set the angular velocity initially to zero, since the velocity is also
% initially set to zero
L = zeros(3,1);

% initialize velocity center of mass
Ucm = zeros(3,1);

% here is the timestep loop
centroid_position = zeros(length(timevec), NDIM);
for t = 1:length(timevec)

    % this just displays the current simulation time in the matlab window
    disp(timevec(t));	

    % compute center of mass and Xtwiddle
    Xcm = mean((M.*X));
    Xtwiddle = X - Xcm;	       

    % compute moment of inertia tensor
    I = zeros(3,3);
    for kk = 1:num_nodes
	I = M(kk).*( norm(X(kk,:))^2.*eye(3) - X(kk,:)'*X(kk,:) );
    end     

    % compute the angular velocity
    Omega = I\L;

    % normalize angular velocity if it is nonzero
    if(Omega > 100*eps)
         unit_Omega = Omega/norm(Omega);
         Omega_cross = [0 -Omega(3) Omega(2); Omega(3) 0 -Omega(1); -Omega(2) Omega(1) 0];
         P_Omega = unit_Omega*unit_Omega';
	 Xtwiddle = ( P_Omega*(Xtwiddle') + cos(norm(Omega)*dt).*(eye(3) - P_Omega)*(Xtwidde') + sin(norm(Omega)*dt).*(Omega_cross*(Xtwiddle'))./norm(Omega) )'
    end      

    % update the position and velocity for each node
    sum(M)*U = U + (dt*F./[M M]) + dt*F_gravity + dt*F_ground(X,Sg);
    X = X + dt*U;

    % store the position of centroid
    centroid_position(t,:) = X(n+1,:);

    % plot the current position of the ball
    figure(2);
    x = [X(jj,1) X(kk,1)];
    y = [X(jj,2) X(kk,2)];
    plot(x',y','linewidth',4)
    xlim([0 2])
    ylim([-3 12])
    axis equal
    pause(0.001)

    % stop simulation if it blows up.
    if(abs(X)> 1e16)
      break
    end

end

figure(3); hold on
plot(timevec, centroid_position,'linewidth',2)

% function which defines the force on the ground, existing on the plane X_2 = 0
function Fg = F_ground(X,Sg)
    Fg = zeros(size(X));
    Fg(:,2) = (X(:,2) < 0).*(-Sg*X(:,2));
end

