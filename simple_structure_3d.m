close all; clear all; clc;
% building a simple structure in 3d

num_nodes = 8;
num_links = 12;
jj = zeros(num_links,1);
kk = zeros(num_links,1);

% set positions of nodes
X = zeros(num_nodes,3);
X(1,1) = 0; % node 1
X(1,2) = 0;
X(1,3) = 0;
X(2,1) = 0; % node 2
X(2,2) = 1;
X(2,3) = 0;
X(3,1) = 1; % node 3
X(3,2) = 1;
X(3,3) = 0;
X(4,1) = 1; % node 4
X(4,2) = 0;
X(4,3) = 0;
X(5,1) = 0; % node 5
X(5,2) = 0;
X(5,3) = 1;
X(6,1) = 0; % node 6
X(6,2) = 1;
X(6,3) = 1;
X(7,1) = 1; % node 7
X(7,2) = 1;
X(7,3) = 1;
X(8,1) = 1; % node 8
X(8,2) = 0;
X(8,3) = 1;

% build structure by creating links
jj(1) = 1; % link 1
kk(1) = 2;
jj(2) = 2; % link 2
kk(2) = 3;
jj(3) = 3; % link 3
kk(3) = 4;
jj(4) = 4; % link 4
kk(4) = 1;
jj(5) = 5; % link 5
kk(5) = 6;
jj(6) = 6; % link 6
kk(6) = 7;
jj(7) = 7; % link 7
kk(7) = 8;
jj(8) = 8; % link 8
kk(8) = 5;
jj(9) = 2; % link 9
kk(9) = 6;
jj(10) = 3; % link 10
kk(10) = 7;
jj(11) = 1; % link 11
kk(11) = 5;
jj(12) = 4; % link 12
kk(12) = 8;

figure(1);
x = [X(jj,1) X(kk,1)];
y = [X(jj,2) X(kk,2)];
z = [X(jj,3) X(kk,3)];
plot3(x',y',z','linewidth',4)



