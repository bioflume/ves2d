function couette_driver(runName, VF, VC1, VC2, VC1ratio, Vinf, resume_file, VC_file)

fprintf('Multiple elliptical vesicles in a couette apparatus.\n');
options.farField = 'couette'; % background velocity
options.farFieldSpeed = Vinf; % scaling of background velocity
prams.nv = 125;    % number of vesicles
options.confined = true; % confined or unbounded geometry

% Spatial-resolution
prams.N = 96;    % points per vesicle
prams.Nbd = 512; % points on a wall
prams.nvbd = 2;  % number of solid walls

% physical properties of vesicles
prams.kappa = 1e-1;   % bending coefficient
prams.errorTol = 1e-2;

% parameters for numerics
options.fmm = ~false;  % fmm for single-layer potentials
options.fmmDLP = ~false; % fmm for double-layer potentials
options.matFreeWalls = false; % W2W interactions are done without a matrix

prams.gmresTol = 1e-10;  % tolerance for gmres

% Low-Resolution Correction Algorithms
options.antiAlias = true; % upsampling for anti-aliasing
options.reparameterization = true; % reparametrization of membrane
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA
options.repulsion = true; % repulsion between vesicles
prams.minDist = 1; %0.3
options.usePlot = false;
% Temporal Resolution (parameters for new implementation)
prams.T = 100/Vinf;   % time horizon (or 100*lenScale/options.farFieldSpeed) etc.
prams.m = (prams.T/1e-3);
options.timeAdap = true; % time adaptivity
prams.dtMax = 1e-3;
prams.dtMin = 1e-5;
prams.betaInc = 1e-2;
prams.betaDec = 5e-2;
prams.betaUp = 1.5;
prams.betaDown = 0.5;
prams.alpha = 0.9;
prams.rtolArea = 1e-3;
prams.rtolLength = 1e-3;

% Save vesicle information and create a log file
% Name of log file for saving messages
options.logFile = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/couette/' runName '.log'];
% Name of binary data file for storing vesicle information
options.dataFile = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/couette/' runName 'Data.bin'];

% Set options and parameters that the user doesn't
% Also add src to path
[options,prams] = initVes2D(options,prams);

% Initial configuration of vesicles
oc = curve;
% Build solid walls
Xwalls = oc.initConfig(prams.Nbd,...
    options.farField,'nv',prams.nvbd,'center',[[0;0] [0;0]]);

if isempty(resume_file)
omVesGen = monitor([],options,prams);
ttVesGen = tstep(options,prams,omVesGen);
[X,prams.nv] = oc.initConfig(prams.N,'scale',1,'volFrac',VF,...
     Xwalls,ttVesGen);
 
nVC1 = ceil(VC1ratio*prams.nv); 
prams.viscCont = VC2*ones(prams.nv, 1);
idVC1 = sort(randperm(prams.nv, nVC1));
prams.viscCont(idVC1) = VC1;
VCs = prams.viscCont;
VCfileName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/couette/' runName '_VCvalues.mat'];
save(VCfileName,'VCs', 'VC1', 'VC2','VF','VC1ratio','nVC1','idVC1','Vinf')
else
VCfileName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/couette/' VC_file '_VCvalues.mat'];    
load(VCfileName)
load(['/mnt/home/gkabacaoglu/codes/ves2dn/examples/couetteResumeFiles/' resume_file '.mat'])
prams.nv = numel(X(1,:));
prams.viscCont = VC2*ones(prams.nv,1);
prams.viscCont(idVC1) = VC1;
end
% Run vesicle code
Xfinal = Ves2D(X,Xwalls,[],[],prams,options,[],[]);
end



