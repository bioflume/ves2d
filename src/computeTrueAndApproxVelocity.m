function [Vinf,trueVnear,trueVfar,trueVself,sigmaTrue,appVnear,appVfar,...
    appVself,sigmaApp,sigDen] = computeTrueAndApproxVelocity...
    (X,prams,options,tt,iAppOnly,iExOnly,sigDenOld)

op = tt.op;

N = prams.N;   % number of points
nv = prams.nv; % number of vesicles

vesicle = capsules(X,zeros(N,nv),[],prams.kappa,prams.viscCont,options.antiAlias);
vesicle.setUpRate(op);

if iExOnly
% Compute true tension
tt.gmresTol = 1e-8;
[sigmaTrue,~,~,~,~,~,iter,tt,sigDen] = vesicle.computeSigAndEtaTToutInEx(tt,[],[],[]);
disp(['GMRES takes ' num2str(iter) ' iterations.'])

% Self-interaction with singular quadrature
fTrue = vesicle.tracJump(X,sigmaTrue); % traction jump
trueVself = op.exactStokesSLdiag(vesicle,tt.Galpert,fTrue);

% Compute true near and far velocity 
[~,trueVnear,trueVfar,~] = op.nearSingIntAll(vesicle,fTrue,...
    @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X),...
    tt.SLPnoCorr,tt.NearV2V,@op.exactStokesSLfmm,@op.exactStokesSL,...
    vesicle,true);

else
trueVnear = [];
trueVfar = [];
trueVself = [];
sigmaTrue = [];
end

if iAppOnly
tt.NearV2V = vesicle.getZone([],1);
% Compute approximate tension
[sigmaApp,~,iter,tt,sigDen] = vesicle.computeSigmaWRegul(tt,sigDenOld);
disp(['GMRES takes ' num2str(iter) ' iterations.'])

% Self interaction with regularized kernel
fApp = vesicle.tracJump(X,sigmaApp); % traction jump
appVself = op.exactStokesSLdiag(vesicle,tt.SLPnoCorr,fApp);

% Single layer potential (far+near) w/ FMM (near singular integrals are wrong)
[~,appVnear,appVfar] = op.divideNearFarSLP(...
        vesicle,fApp,tt.SLPnoCorr,tt.NearV2V,@op.exactStokesSLfmm,...
        @op.exactStokesSL,vesicle,true);
else
appVnear = [];
appVfar = [];
appVself = [];
sigmaApp = [];
end
% Background velocity
Vinf = tt.farField(X,[]);    
