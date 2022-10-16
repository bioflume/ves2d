function DLDStreamRunsWpost(runName,folderName,Xpost,wbox,Dx,Dy,theta,...
    Next,Nint,Nves,Ufar,VCmaj,VCmin,seedRate,useFMM,totnv)

options.farField = 'DLD'; % background velocity
options.farFieldSpeed = Ufar; % scaling of background velocity
prams.nv = 1; % number of vesicles
options.confined = true; % confined or unbounded geometry
options.diffDiscWalls = true; % walls are discretized by different Nbd

% Descriptive run name
prams.runName = runName; 
% This is the name given to all related files 

options.usePlot = false; % Plot on-the-fly
options.track = false; % trackers on membrane

% Keep seeding vesicles in the entrance
options.streaming = true;
% Freeze vesicles coming to the end of the device
options.freezing = true;

if options.streaming
  % 1/scaling for area to check before streaming
  prams.streamRate  = seedRate;
  % Number of vesicles to be streamed
  prams.nSeed = 2;
  % Number of vesicles will be streamed (including the initial ones)
  prams.totnv = totnv;
  % Viscosity contrast of the streamed vesicles (2 different)
  % the majority of the vesicles haS VC of prams.vesViscCont
  % the minority has VC of prams.ves2ViscCont
  prams.ves2ViscCont = VCmin;
  % concentration of the other vesicle; this is the percentage of the
  % total number of vesicles with prams.ves2ViscCont streamed. 
  prams.ves2Concent = 0.1; 
  % VCsXaY means that we want to separate one with VC = X from the other
  % one (so VC=X is the minority)
%   !mkdir ./output/streamSim_VCs5a1_rateP4/
%   prams.folderName = './output/streamSim_VCs5a1_rateP4/';
  prams.folderName = folderName;
  options.saveData = false; % do not save data in the usual way
  options.freezing = true; % we have to freeze if we keep streaming
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Spatial Resolution
prams.N      = Nves; % # of points per vesicle

prams.NbdInt = Nint; % # of points on posts % 64
prams.NbdExt = Next; % # of points on exterior wall % 3072 for 6 rows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% GEOMETRY
% Sizes
scale       = 0.1;
if isempty(Xpost)
    
  prams.Dpost = 15*scale; % diameter of a post in micron * scale
  Dpostx = prams.Dpost/2;
  Dposty = prams.Dpost/2;
  prams.Dy    = 10*scale; % vertical separation in micron * scale
  prams.Dx    = 10*scale; % horizontal separation in micron * scale
else
  Xpost = [interpft(Xpost(1:end/2),Nint);interpft(Xpost(end/2+1:end),Nint)];
  Dpostx = max(Xpost(1:end/2))-min(Xpost(1:end/2)); 
  Dposty = max(Xpost(end/2+1:end))-min(Xpost(end/2+1:end)); 
  if isempty(Dx)
    prams.Dx = wbox-Dpostx;
    prams.Dy = wbox-Dposty;
  else
    prams.Dx = Dx;
    prams.Dy = Dy;
  end
  prams.Dpost = [];
end


% DLD geometry
% Period of a device
prams.periodN = (Dposty+Dy)/(tan(theta)*(Dpostx+Dx)); 
% Row-shift fraction
prams.epsilon = 1/prams.periodN; 

% number of obstacles in x direction  (number of columns)
prams.ncol = ceil(1.5*prams.periodN);
% number of obstacles in y direction (number of rows)
prams.nrow    = 6; 

% Vesicle geometry
DeffVes     = 3; % short-axis thickness in microns
reducedArea = 0.65; % vesicle reduced area
scaleVes    = scale*DeffVes/2.23;
prams.scaleVes = scaleVes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Vesicle
prams.kappa    = 1e-2; % bending coefficient
prams.vesViscCont = VCmaj; % Vesicle viscosity contrast

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for numerics
prams.gmresTol = 1e-7;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

options.fmm = useFMM;  % fmm for single-layer potentials
options.fmmDLP = useFMM; % fmm for double-layer potentials
options.matFreeWalls = false; % W2W interactions are done without a matrix

options.saveWallMat  = false; % save exact factorization on disk
options.haveWallMats = false; % Do we have wall matrices computed already?

options.outOfCore    = false; % Compute inexact factorization out-of-core
                              % Also block applying exact inverse
% if exact factorization and out-of-core, then we block applying exact inverse
options.memsize      = 1;     % memory usage in gb for out-of-core

options.fastDirect   = false;  % Use fast direct solver
options.HODLRforW2W  = false;   % Use HODLR to compress wall2wall interactions
prams.lev_max        = 3;     % Maximum level to go in HODLR
prams.etolFD         = 1e-8; % tolerance in hodlr to decide l

%if 1
%  !mkdir /workspace/gokberk/DLD_Matrices/DLDN20/
%  options.wallMatFile = '/workspace/gokberk/DLD_Matrices/DLDN20/wallMat';
%else
%  !mkdir ./output/LargeMatDLD/
%  options.wallMatFile = './output/LargeMatDLD/wallMat';
%end
% directory where we save large wall matrices
options.wallMatFile = '/workspace/gokberk/wallMatNew';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.antiAlias = true; % upsampling for anti-aliasing
options.reparameterization = true; % reparametrization of membrane
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA
options.repulsion = true; % repulsion between vesicles, and vesicles-walls
% Minimum distance (ratio of max. spacing) after which repulsion is active
prams.minDist      = 1; 

options.timeAdap = true; % time adaptivity
options.nsdc = 0;  % number of sdc corrections
prams.areaLenTol = 5e-4; % tolerance for errors in area-length for time-stepping

if ~options.streaming
  options.saveData = true;    
end
% Name of log file
options.logFile  = ['output/' prams.runName '.log'];
% Name of binary data file for storing vesicle information
options.dataFile = ['output/' prams.runName '_Data.bin'];

% Set options and parameters that the user doesn't
[options,prams] = initVes2D(options,prams);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE WALLS
oc = curve;
% Generate the DLD device

[XwallsInt,XwallsExt,L,H,prams.xfreeze,radiusRight,radiusUp,...
    prams.Dpostx,prams.Dposty] = ...
    oc.initConfigDLD('wholeDLD',prams.NbdInt,prams.NbdExt,Xpost,prams.Dpost,...
    prams.epsilon,prams.periodN,prams.Dx,prams.Dy,prams.nrow,prams.ncol);

prams.nvbdInt = size(XwallsInt,2);
prams.nvbdExt = 1;
prams.nvbd = prams.nvbdInt + 1;

[extWallx,extWally] = oc.getXY(XwallsExt);

xmin = min(extWallx(:)); xmax = max(extWallx(:));
ymin = min(extWally(:)); ymax = max(extWally(:));
options.axis = [xmin-0.5 xmax+0.5 ymin-0.5 ymax+0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE VESICLES


% randomly place vesicles and rigid particles
omVesGen = monitor([],options,prams);
ttVesGen = tstep(options,prams,omVesGen);

% range for x and y coordinates where we generate vesicles
prams.xrange = xmin + radiusRight + [0 0.75*prams.Dy];
prams.yrange = [-0.45 0.1]*prams.Dy;

load relaxed64.mat
X = [interpft(X(1:end/2),prams.N);interpft(X(end/2+1:end),prams.N)];
IA = oc.getIncAngle(X);
Xnew = zeros(size(X));
for i = 1 : prams.N
Xnew(i) = cos(-IA)*X(i)-sin(-IA)*X(i+prams.N);
Xnew(i+prams.N) = sin(-IA)*X(i)+cos(-IA)*X(i+prams.N);
end

prams.nv = 2;
[X,volFrac] = oc.initConfig(prams.N,'createVesforDLD',...
'nv',prams.nv,'scale',scaleVes,...
'reducedArea',reducedArea,'refVesicle',Xnew,...
'randomlyPlaceVes',XwallsExt,XwallsInt,...
prams.xrange,prams.yrange,ttVesGen,'DLDEntrance',prams,0.05);  
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


% Get the length and velocity scale of simulation to decide on params.
[~,~,arcLen] = oc.geomProp(X); 
lenScale = max(arcLen); % length scale of the simulation

% Temporal Resolution
prams.T = 5e+4;   % time horizon
% time scale for simulation
% # of time steps (if constant time stepping) or Dt for the first step
% (adaptive)
prams.m = ceil(prams.T/(0.01*lenScale/options.farFieldSpeed)); 
prams.dtMax = 0.1*lenScale/options.farFieldSpeed; % maximum allowable time-step size
prams.dtMin = 1e-4*lenScale/options.farFieldSpeed; % minimum allowable time-step size  

% first compute velocity at the middle and scale the maximum velocity
gapX = -L/2+[2;3;4;5;6;7;8]*(prams.Dpostx+prams.Dx);
gapY = [1;2;3;4;5;6;7]*tan(theta)*(prams.Dpostx+prams.Dx);
gapY1 = gapY+prams.Dposty+prams.Dy;
gapY2 = gapY-prams.Dposty-prams.Dy;


% 4) Construct tstep structure for both DLD models and get BCs
om = monitor(X,options,prams);
tt = tstep(options,prams,om);
  
[vel,~,~] = tt.computeVelFieldNoVes(prams,[],XwallsInt,XwallsExt,...
  [gapX;gapX;gapX;gapY1;gapY;gapY2],0);
scaling = options.farFieldSpeed/mean(vel(1:end/2));
options.farFieldSpeed = options.farFieldSpeed*scaling; 

tt.farField = @(X,Xint) tt.bgFlow(X,options.farField,...
  'Speed',options.farFieldSpeed,'intWalls',XwallsInt,'nrow',prams.nrow,...
  'ncol',prams.ncol,'Dpostx',prams.Dpostx,'Dposty',prams.Dposty,'GapX',prams.Dx,...
  'GapY',prams.Dy,'epsilon',prams.epsilon);



% Run vesicle code  
Xfinal = Ves2D(X,[],XwallsInt,XwallsExt,prams,options,[],[],[]);
  
