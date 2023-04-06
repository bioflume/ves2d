clear; clc;
addpath ../src/
% Three rigid bodies

N = 32;
theta = [0:N-1]'/N*2*pi;
X1 = [cos(theta)-4; sin(theta)]; % at (-4,0) with r = 1;
X2 = [2*cos(theta)+4; 2*sin(theta)+4]; % at (4,4) with r = 2;
X3 = [cos(theta)+5; sin(theta)-4]; % at (5, -4) with r = 1;

bb = body([X1 X2 X3], [1;2;1], [-4 0; 4 4; 5 -4],[],[]);
bb1 = body(X1,[],[],[],[]);
bb2 = body(X2,[],[],[],[]);
bb3 = body(X3,[],[],[],[]);

density = [rand(2*N,1) rand(2*N,1) rand(2*N,1)];


ker = kernels(N);

SLP_12 = ker.stokesSLmatrixInteractionWeightless(bb1,bb2,1);
SLP_13 = ker.stokesSLmatrixInteractionWeightless(bb1,bb3,1);

SLP_21 = ker.stokesSLmatrixInteractionWeightless(bb2,bb1,1);
SLP_23 = ker.stokesSLmatrixInteractionWeightless(bb2,bb3,1);

SLP_31 = ker.stokesSLmatrixInteractionWeightless(bb3,bb1,1);
SLP_32 = ker.stokesSLmatrixInteractionWeightless(bb3,bb2,1);


DLP_12 = ker.stokesDLmatrixInteractionWeightless(bb1,bb2);
DLP_13 = ker.stokesDLmatrixInteractionWeightless(bb1,bb3);

DLP_21 = ker.stokesDLmatrixInteractionWeightless(bb2,bb1);
DLP_23 = ker.stokesDLmatrixInteractionWeightless(bb2,bb3);

DLP_31 = ker.stokesDLmatrixInteractionWeightless(bb3,bb1);
DLP_32 = ker.stokesDLmatrixInteractionWeightless(bb3,bb2);


DLPT_12 = ker.stokesDLTmatrixInteractionWeightless(bb1,bb2);
DLPT_13 = ker.stokesDLTmatrixInteractionWeightless(bb1,bb3);

DLPT_21 = ker.stokesDLTmatrixInteractionWeightless(bb2,bb1);
DLPT_23 = ker.stokesDLTmatrixInteractionWeightless(bb2,bb3);

DLPT_31 = ker.stokesDLTmatrixInteractionWeightless(bb3,bb1);
DLPT_32 = ker.stokesDLTmatrixInteractionWeightless(bb3,bb2);

SL_on1 = SLP_21*density(:,2) + SLP_31*density(:,3);
SL_on2 = SLP_12*density(:,1) + SLP_32*density(:,3);
SL_on3 = SLP_13*density(:,1) + SLP_23*density(:,2);

DL_on1 = DLP_21*density(:,2) + DLP_31*density(:,3);
DL_on2 = DLP_12*density(:,1) + DLP_32*density(:,3);
DL_on3 = DLP_13*density(:,1) + DLP_23*density(:,2);

DLT_on1 = DLPT_21*density(:,2) + DLPT_31*density(:,3);
DLT_on2 = DLPT_12*density(:,1) + DLPT_32*density(:,3);
DLT_on3 = DLPT_13*density(:,1) + DLPT_23*density(:,2);

DL_on_all = ker.stokesDL_times_density(bb, bb, density,1);
DLT_on_all = ker.stokesDLT_times_density(bb,bb,density,1);
SL_on_all = ker.stokesSL_times_density(bb,bb,1,density,1);

errInDL1 = norm(DL_on1-DL_on_all(:,1))/norm(DL_on1);
errInDL2 = norm(DL_on2-DL_on_all(:,2))/norm(DL_on2);
errInDL3 = norm(DL_on3-DL_on_all(:,3))/norm(DL_on3);

errInDLT1 = norm(DLT_on1-DLT_on_all(:,1))/norm(DLT_on1);
errInDLT2 = norm(DLT_on2-DLT_on_all(:,2))/norm(DLT_on2);
errInDLT3 = norm(DLT_on3-DLT_on_all(:,3))/norm(DLT_on3);

errInSL1 = norm(SL_on1-SL_on_all(:,1))/norm(SL_on1);
errInSL2 = norm(SL_on2-SL_on_all(:,2))/norm(SL_on2);
errInSL3 = norm(SL_on3-SL_on_all(:,3))/norm(SL_on3);

[stokesDLP,~] = ker.exactStokesDL(bb,density,[]);


