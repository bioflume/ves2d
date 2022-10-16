function [appVself,appVnear,appVfar,tt] = computeSLPbend(X,prams,options,tt)
op = tt.op;

N = prams.N;   % number of points
nv = prams.nv; % number of vesicles

vesicle = capsules(X,zeros(N,nv),[],prams.kappa,prams.viscCont,options.antiAlias);
vesicle.setUpRate(op);

tt.SLPnoCorr = op.stokesSLmatrixNoCorr(vesicle);

% fApp = vesicle.tracJump(X,zeros(N,nv)); % traction jump due to bending
% force is the second derivative of the shape
IK = fft1.modes(N,nv);
Dx = fft1.diffFT(X(1:end/2,:),IK);
Dy = fft1.diffFT(X(end/2+1:end,:),IK);
DDx = curve.arcDeriv(Dx,1,ones(N,nv),IK);
DDy = curve.arcDeriv(Dy,1,ones(N,nv),IK);
fApp = [DDx;DDy];

% Self interaction 
appVself = op.exactStokesSLdiag(vesicle,tt.SLPnoCorr,fApp);

% Single layer potential (far+near) w/ FMM (near singular integrals are wrong)
tt.NearV2V = vesicle.getZone([],1);
if tt.fmm
  kernel = @op.exactStokesSLfmm;
else
  kernel = @op.exactStokesSL;
end
[~,appVnear,appVfar] = op.divideNearFarSLP(vesicle,fApp,tt.SLPnoCorr,...
  tt.NearV2V,kernel,@op.exactStokesSL,vesicle,true);

