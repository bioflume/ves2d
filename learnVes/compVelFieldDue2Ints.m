clear; clc;
addpath ../src/
load ./output/32modes_TaylorGreen_50Ves_NearNet

vinf = @(X) 2000*[sin(X(1:end/2,:)/2.5 * pi).*cos(X(end/2+1:end,:)/2.5 * pi);-...
    cos(X(1:end/2,:)/2.5 * pi).*sin(X(end/2+1:end,:)/2.5 * pi)]; % Taylor-Green

oc = curve;
N = 32;
nv = 48;
allIdcs = [1:46];

op = poten(N);
vInfs = zeros(2*N,2,251);
%%
for ij = 1 : 251
X = [vesx(:,:,449+ij);vesy(:,:,449+ij)];
vback = vinf(X);
vesicle = capsules(X,[],[],1,1,0);

fBend = vesicle.tracJump(X,zeros(N,nv));
tracJump = fBend;

SLPnoCorr = []; % it would SLPdiag if fmm is on
G = op.stokesSLmatrix(vesicle);

% Get the near structure (this will be done using NN in Python)
NearV2V = vesicle.getZone([],1);    

kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;

SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,NearV2V,...
    kernel,kernelDirect,vesicle,true,false);

tenNew = zeros(N,nv);
G = op.stokesSLmatrix(vesicle);
[~,Ten,Div] = vesicle.computeDerivs;
for k = 1 : nv
  LHS = (Div(:,:,k)*G(:,:,k)*Ten(:,:,k));
  selfBend = G(:,:,k)*fBend(:,k);
  RHS = -Div(:,:,k)*(vback(:,k)+farFieldtracJump(:,k)+selfBend);
  tenNew(:,k) = LHS\RHS;
end % k = 1 : nv

% update the traction jump calculation
fTen = vesicle.tracJump(zeros(2*N,nv), tenNew); 
tracJump = fBend + fTen;

vesicleS = capsules(X(:,allIdcs),[],[],1,1,0);
vesicleT = capsules(X(:,[47;48]),[],[],1,1,0);

[~,NearVS2VT] = vesicleS.getZone(vesicleT,2);  

GS = op.stokesSLmatrix(vesicleS);
SLP = @(X) op.exactStokesSLdiag(vesicleS,GS,X);
vInfs(:,:,ij) = op.nearSingInt(vesicleS,tracJump(:,allIdcs),SLP,SLPnoCorr,NearVS2VT,...
    kernel,kernelDirect,vesicleT,false,false);

end

