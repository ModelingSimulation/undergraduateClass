close all; clear all; clc;
% building a simple structure in 2d

num_nodes = 4;
num_links = 6;
jj = zeros(num_links,1);
kk = zeros(num_links,1);

% set positions of nodes
X = zeros(num_nodes,2);
X(1,1) = 0; % node 1
X(1,2) = 0;
X(2,1) = 0; % node 2
X(2,2) = 1;
X(3,1) = 1; % node 3
X(3,2) = 1;
X(4,1) = 1; % node 4
X(4,2) = 0;

% build structure by creating links
jj(1) = 1; % link 1
kk(1) = 2;
jj(2) = 2; % link 2
kk(2) = 3;
jj(3) = 3; % link 3
kk(3) = 4;
jj(4) = 4; % link 4
kk(4) = 1;
jj(5) = 1; % link 5
kk(5) = 3;
jj(6) = 2; % link 6
kk(6) = 4;

figure(1);
x = [X(jj,1) X(kk,1)];
y = [X(jj,2) X(kk,2)];
plot(x',y','linewidth',4)



