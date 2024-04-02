clear; clc;
addpath ./output/polyDisp_VF20Cap75_circPost/

file = 'tStep';
nsteps = 4000;
nskip = 100;

% Load the wall matrices
load([file num2str(1)])

extWallx = [XwallsExt(1:end/2);XwallsExt(1)];
extWally = [XwallsExt(end/2+1:end);XwallsExt(end/2+1)];

intWallx = [XwallsInt(1:end/2,:);XwallsInt(1,:)];
intWally = [XwallsInt(end/2+1:end,:);XwallsInt(end/2+1,:)];

for it = nskip : nskip : nsteps
load([file num2str(it)])
rigx = [Xhist(1:end/2,1);Xhist(1,1)];
rigy = [Xhist(end/2+1:end,1);Xhist(end/2+1,1)];

vesx = [Xhist(1:end/2,2:end);Xhist(1,2:end)];
vesy = [Xhist(end/2+1:end,2:end);Xhist(end/2+1,2:end)];


figure(1);clf;hold on;
plot(extWallx, extWally, 'k','linewidth',2)
plot(intWallx, intWally, 'k','linewidth',2)
plot(rigx, rigy, 'b','linewidth',2)
plot(vesx, vesy, 'r','linewidth',2)

axis equal

box on
pause(0.1)





end