% VEM code for 2D piezoelasticity
% Our paper: Virtual element method for piezoelasticity with contact
% Copyright (C) Bing-Bing Xu
clear;
load("matlab.mat");

sumNode = size(node,1);
p = 100; % press

nodeL = find(node(:,1)<0.001);
sL = size(nodeL,1);
nodeR  = find(node(:,1)>48-0.001);
sR = size(nodeR,1);
% find face
pface = findFace(node,elem,nodeR);
press = [pface,ones(length(pface),1)*p];
fixNode = [nodeL,ones(sL,1),zeros(sL,1);nodeL,2*ones(sL,1),zeros(sL,1);nodeL,3*ones(sL,1),zeros(sL,1)];
fixMes = [fixNode(:,1)+(fixNode(:,2)-1)*sumNode,fixNode(:,3)];

% calculate nodal force
nodeForce = getForce(node,sumNode,elem,press,'y');

d = [126e9,743e8,0;743e8,115e9,0;0,0,256e8];
e = [0,-5.2;0,15.1;12.7,0];
rmt = [-6.463e-9,0;0,-5.622e-9];

GK = globalK(node,elem,d,e,rmt);

[GK,F] = boundaryCondition(GK,nodeForce,fixMes(:,1),fixMes(:,2),1);

uh = GK\F;

uh = full(uh);

ux = uh(1:sumNode);
uy = uh(sumNode+1:2*sumNode);
fai = uh(2*sumNode+1:end);

% figure;
% showsolution(node,elem,ux);

% figure;
% showsolution(node,elem,uy);

figure;
showsolution(node,elem,fai);











