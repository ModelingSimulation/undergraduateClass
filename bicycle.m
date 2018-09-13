close all; clear all; clc;
% building a bicycle

% physical parameters
L = 10;
LA = 4;
LBC = 2;
HAB = 2;
HC = 0.2;
n = 30; % number of nodes on rim
W = 1.0;
RW = 1.5;

% initialize arrays
num_nodes = 2*n + 7;
num_links = 6*n + 17;
jj = zeros(num_links,1);
kk = zeros(num_links,1);
X = zeros(num_nodes,3);

% making the nodes on the rims
for k = 1:n
    theta = 2*pi*k/n;
    X(k,:) = RW*[cos(theta), 0, sin(theta)];
    X(k+n,:) = X(k,:) + [L,0,0];
end

% naming indices for the different nodes
iA = 2*n + 1;
iB = 2*n + 2;
iC = 2*n + 3;
iD = 2*n + 4;
iDpr = 2*n + 5;
iE = 2*n + 6;
iEpr = 2*n + 7;

% setting the positions of the other nodes
X(iA,:) = [LA, 0, HAB];
X(iB,:) = [L-LBC, 0, HAB];
X(iC,:) = [L-LBC, 0, HC];
X(iD,:) = [0, -W/2, 0];
X(iDpr,:) = [0, W/2, 0];
X(iE,:) = [L, -W/2, 0];
X(iEpr,:) = [L, W/2, 0];

% naming indices for links
spokesD = 1:n;
spokesDpr = (n+1):2*n;
spokesE = (2*n+1):3*n;
spokesEpr = (3*n+1):4*n;
rimlinksD = (4*n+1):5*n;
rimlinksE = (5*n+1):6*n;
framelinks = (6*n+1):(6*n+15);
axels = (6*n+16):(6*n+17);

% links for spokes
jj(spokesD) = 1:n; kk(spokesD) = iD;
jj(spokesDpr) = 1:n; kk(spokesDpr) = iDpr;
jj(spokesE) = (n+1):2*n; kk(spokesE) = iE;
jj(spokesEpr) = (n+1):2*n; kk(spokesEpr) = iEpr;

% links for frame
jj(framelinks) = [iD iD iD iE iE iE iDpr iDpr iDpr iEpr iEpr iEpr iA iB iC];
kk(framelinks) = [iA iB iC iA iB iC iA iB iC iA iB iC iB iC iA];

% links for the rims
jj(rimlinksD) = [1:n]; kk(rimlinksD) = [2:n 1];
jj(rimlinksE) = [n+1:2*n]; kk(rimlinksE) = [n+2:2*n n+1];

% links for the axels
jj(axels) = [iD iE]; kk(axels) = [iDpr iEpr];

figure(1);
x = [X(jj,1) X(kk,1)];
y = [X(jj,2) X(kk,2)];
z = [X(jj,3) X(kk,3)];
plot3(x',y',z','linewidth',4)

% make a movie
v = VideoWriter('test');
open(v);
num_angles = 500;
frames(num_angles) = struct('cdata',[],'colormap',[]);
for tt = 1:num_angles
    angle = tt*360/num_angles;
    plot3(x',y',z','linewidth',4)
    axis equal
    view(angle, 30);
    frame = getframe(gcf);
    frames(tt) = getframe;
    writeVideo(v,frame);
end
close(v);

