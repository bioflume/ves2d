clear; clc;

set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

load shanSim_matlab_version.mat


nsteps = numel(X(:,1,1));

for k = 1 :  nsteps

figure(1); clf;


xvec = [X(k,1:end/2,1)'; X(k,1,1)'];
yvec = [X(k,end/2+1:end,1)';X(k,end/2+1,1)'];
plot(xvec, yvec, 'Color',[202,0,32]/255,'linewidth',2)
hold on
hFill = fill(xvec, yvec, [202,0,32]/255);
set(hFill,'EdgeColor', [202,0,32]/255);


xvec = [X(k,1:end/2,2)'; X(k,1,2)'];
yvec = [X(k,end/2+1:end,2)';X(k,end/2+1,2)'];
plot(xvec, yvec, 'Color',[5,113,176]/255,'linewidth',2)
hold on
hFill = fill(xvec, yvec, [5,113,176]/255);
hFill.FaceAlpha = 0.5;
set(hFill,'EdgeColor', [5,113,176]/255);

box on

axis equal
xlim([-0.75 1])
ylim([-0.35 0.35])

title(k)

pause(0.1)
end


%%
addpath ../../src/

istep = 3;

load shanSim_matlab_version.mat
tsteps = [1;2;105;106;150;151];
shan_tracJump1 = tracJump1;
shan_tracJump2 = tracJump2;
shan_vTrac1 = farFieldtracJump1;
shan_vTrac2 = farFieldtracJump2;
shan_MVinf = MVinf;

Xold = reshape(X(tsteps(istep),:,:),256,2);
tenOld = reshape(ten(tsteps(istep),:,:),128,2);
vesicle = capsules(Xold,[],[],1,1,0);
tracJump = vesicle.tracJump(Xold, tenOld);


shan_Xnew = X(tsteps(istep)+1,:,:);
shan_tenNew = ten(tsteps(istep)+1,:,:);
load(['CheckingShansNet_istep' num2str(istep) '.mat']);
% corresponds to tsteps(istep)+1 in Shans

% TRACJump1
figure(1);clf;
plot(Xold(1:end/2,1),Xold(end/2+1:end,1),'k','linewidth',2)
hold on
quiver(Xold(1:end/2,1),Xold(end/2+1:end,1),tracJump1(1:end/2,1),tracJump1(end/2+1:end,1),0.7)
quiver(Xold(1:end/2,1),Xold(end/2+1:end,1),shan_tracJump1(tsteps(istep)+1,1:end/2,1)',shan_tracJump1(tsteps(istep)+1,end/2+1:end,1)',0.7)
axis equal
legend('vesicle','matlab','python')
title('TracJump1-vesicle1')

figure(2);clf;
plot(Xold(1:end/2,2),Xold(end/2+1:end,2),'k','linewidth',2)
hold on
quiver(Xold(1:end/2,2),Xold(end/2+1:end,2),tracJump1(1:end/2,2),tracJump1(end/2+1:end,2),0.7)
quiver(Xold(1:end/2,2),Xold(end/2+1:end,2),shan_tracJump1(tsteps(istep)+1,1:end/2,2)',shan_tracJump1(tsteps(istep)+1,end/2+1:end,2)',0.7)
axis equal
title('TracJump1-vesicle2')
legend('vesicle','matlab','python')
errTJ1_V1 = norm(tracJump1(:,1)-shan_tracJump1(tsteps(istep)+1,:,1)');
errTJ1_V2 = norm(tracJump1(:,2)-shan_tracJump1(tsteps(istep)+1,:,2)');


% farFieldTRACJump1
figure(3);clf;
plot(Xold(1:end/2,1),Xold(end/2+1:end,1),'k','linewidth',2)
hold on
quiver(Xold(1:end/2,1),Xold(end/2+1:end,1),farFieldtracJump1(1:end/2,1),farFieldtracJump1(end/2+1:end,1),0.7)
quiver(Xold(1:end/2,1),Xold(end/2+1:end,1),shan_vTrac1(tsteps(istep)+1,1:end/2,1)',shan_vTrac1(tsteps(istep)+1,end/2+1:end,1)',0.7)
axis equal
legend('vesicle','matlab','python')
title('farFieldTracJump1-vesicle1')


plot(Xold(1:end/2,2),Xold(end/2+1:end,2),'k','linewidth',2)
hold on
quiver(Xold(1:end/2,2),Xold(end/2+1:end,2),farFieldtracJump1(1:end/2,2),farFieldtracJump1(end/2+1:end,2),0.7)
quiver(Xold(1:end/2,2),Xold(end/2+1:end,2),shan_vTrac1(tsteps(istep)+1,1:end/2,2)',shan_vTrac1(tsteps(istep)+1,end/2+1:end,2)',0.7)
axis equal
title('farFieldTracJump1-vesicle2')
legend('vesicle','matlab','python','vesicle','matlatb','python')

errvT1_V1 = norm(farFieldtracJump1(:,1)-shan_vTrac1(tsteps(istep)+1,:,1)');
errvT1_V2 = norm(farFieldtracJump1(:,2)-shan_vTrac1(tsteps(istep)+1,:,2)');


% TRACJump2
figure(4);clf;
plot(Xold(1:end/2,1),Xold(end/2+1:end,1),'k','linewidth',2)
hold on
quiver(Xold(1:end/2,1),Xold(end/2+1:end,1),tracJump2(1:end/2,1),tracJump2(end/2+1:end,1),0.7)
quiver(Xold(1:end/2,1),Xold(end/2+1:end,1),shan_tracJump2(tsteps(istep)+1,1:end/2,1)',shan_tracJump2(tsteps(istep)+1,end/2+1:end,1)',0.7)
axis equal
legend('vesicle','matlab','python')
title('TracJump2-vesicle1')

figure(5);clf;
plot(Xold(1:end/2,2),Xold(end/2+1:end,2),'k','linewidth',2)
hold on
quiver(Xold(1:end/2,2),Xold(end/2+1:end,2),tracJump2(1:end/2,2),tracJump2(end/2+1:end,2),0.7)
quiver(Xold(1:end/2,2),Xold(end/2+1:end,2),shan_tracJump2(tsteps(istep)+1,1:end/2,2)',shan_tracJump2(tsteps(istep)+1,end/2+1:end,2)',0.7)
axis equal
title('TracJump2-vesicle2')
legend('vesicle','matlab','python')

errTJ2_V1 = norm(tracJump2(:,1)-shan_tracJump2(tsteps(istep)+1,:,1)');
errTJ2_V2 = norm(tracJump2(:,2)-shan_tracJump2(tsteps(istep)+1,:,2)');


% farFieldTRACJump2
figure(6);clf;
plot(Xold(1:end/2,1),Xold(end/2+1:end,1),'k','linewidth',2)
hold on
quiver(Xold(1:end/2,1),Xold(end/2+1:end,1),farFieldtracJump2(1:end/2,1),farFieldtracJump2(end/2+1:end,1),0.7)
quiver(Xold(1:end/2,1),Xold(end/2+1:end,1),shan_vTrac2(tsteps(istep)+1,1:end/2,1)',shan_vTrac2(tsteps(istep)+1,end/2+1:end,1)',0.7)
axis equal
legend('vesicle','matlab','python')
title('farFieldTracJump2')

plot(Xold(1:end/2,2),Xold(end/2+1:end,2),'k','linewidth',2)
hold on
quiver(Xold(1:end/2,2),Xold(end/2+1:end,2),farFieldtracJump2(1:end/2,2),farFieldtracJump2(end/2+1:end,2),0.7)
quiver(Xold(1:end/2,2),Xold(end/2+1:end,2),shan_vTrac2(tsteps(istep)+1,1:end/2,2)',shan_vTrac2(tsteps(istep)+1,end/2+1:end,2)',0.7)
axis equal
title('farFieldTracJump2')
legend('vesicle','matlab','python','vesicle','matlatb','python')

errvT2_V1 = norm(farFieldtracJump2(:,1)-shan_vTrac2(tsteps(istep)+1,:,1)');
errvT2_V2 = norm(farFieldtracJump2(:,2)-shan_vTrac2(tsteps(istep)+1,:,2)');


