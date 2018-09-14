close all; clear all; clc;
% building a string of springs
NDIM = 2;

num_nodes = 3;
num_links = 2;
jj = zeros(num_links,1);
kk = zeros(num_links,1);

% set numerical parameters
dt = 1e-3;
end_time = 5e-1;
timevec = 0:dt:end_time;

% set physical parameters
T = 5e-2;
omega = 2*pi/T;
M = zeros(num_nodes,1); % mass of nodes
S = zeros(num_links,1); % link spring constant
D = zeros(num_links,1); % link dashpot constant
M(1) = 5e-7; 
M(2) = 5e-7; 
M(3) = 5e-7;
S(1) = 0.001;
S(2) = 0.000001;

% set positions of nodes
X = zeros(num_nodes, NDIM);
X(1,1) = 0; % node 1
X(1,2) = 0;
X(2,1) = 1; % node 2
X(2,2) = 0;
X(3,1) = 2; % node 3
X(3,2) = 0;

% build structure by creating links
jj(1) = 1; % link 1
kk(1) = 2;
jj(2) = 2; % link 2
kk(2) = 3;

figure(1);
x = [X(jj,1) X(kk,1)];
y = [X(jj,2) X(kk,2)];
plot(x',y','linewidth',4)

% initialize velocities to zero
U = zeros(num_nodes, NDIM);

% specify free nodes
free_nodes = 2;

% specify bc nodes with prescribed displacement
bc_nodes = [1 3];

% initialize rest length
DX = X(jj,:) - X(kk,:);
Rzero = sqrt(sum(DX.^2,2));

for t = timevec

    % prescribe boundary conditions on nodes
    %X(bc_nodes(1),:) = [0.1*sin(omega*t); 0];
    %U(bc_nodes(1),:) = [omega*0.1*cos(omega*t); 0];
    X(bc_nodes(1),:) = [0; 0.4*sin(omega*t)];
    U(bc_nodes(1),:) = [0; omega*0.4*cos(omega*t)];

    disp(t);
    
    DX = X(jj,:) - X(kk,:);
    DU = U(jj,:) - U(kk,:);  
    R = sqrt(sum(DX.^2,2));
    T = S.*(R - Rzero) + (D./R).*sum(DX.*DU,2);
    TR = T./R;
    FF = [TR TR].*DX;

    F = zeros(num_nodes, NDIM);
    % compute forces on the nodes
    for l = 1:num_links
    	F(kk(l),:) = F(kk(l),:) + FF(l,:);
       	F(jj(l),:) = F(jj(l),:) - FF(l,:);	
    end

    % update the free nodes
    U(free_nodes,:) = U(free_nodes,:) + dt*F(free_nodes,:)./[M(free_nodes) M(free_nodes)];
    X(free_nodes,:) = X(free_nodes,:) + dt*U(free_nodes,:);
   
    figure(2);
    x = [X(jj,1) X(kk,1)];
    y = [X(jj,2) X(kk,2)];
    plot(x',y','linewidth',4)
    xlim([-0.2 2.2])
    ylim([-0.5 0.5])	
    pause(0.1)

end