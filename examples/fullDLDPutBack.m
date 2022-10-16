function fullDLDPutBack(runName,postFile,whichPost,Dx,Dy,Dpost,whichCell,...
    theta,nRepeatPeriod,VCs,rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)
Xpost = [];
addpath ../src/
if ~isempty(postFile)  
  load(postFile)
  oc = curve;
  if isempty(Xpost)
    if whichPost <= numel(bestIds)
      if size(bsplinePoints,1) == 16
        [Xpost,Dpostx,Dposty] = oc.samplePost(NbdPill,reshape(...
          bsplinePoints(:,bestIds(whichPost)),8,2));
      else
        [Xpost,Dpostx,Dposty] = oc.samplePost(NbdPill,bsplinePoints);
      end
    else
      [Xpost,Dpostx,Dposty] = oc.samplePost(NbdPill,reshape(...
          jiggBsplines(:,bestIds(whichPost-numel(bestIds))),8,2));
    end
  else
    Xpost = [interpft(Xpost(1:end/2),NbdPill);interpft(Xpost(end/2+1:end),NbdPill)];
    Dpostx = max(Xpost(1:end/2))-min(Xpost(1:end/2));
    Dposty = max(Xpost(end/2+1:end))-min(Xpost(end/2+1:end));
  end
else
  % a circular pillar  
  Xpost = [];
  Dpostx = Dpost;
  Dposty = Dpost;
end

if isempty(Dx)
  Dx = 2.5-Dpostx;
  Dy = 2.5-Dposty;
end

options.farField = 'DLD'; % background velocity
options.farFieldSpeed = 1.2; % scaling of background velocity
prams.nv = 1; % number of vesicles
options.confined = true; % confined or unbounded geometry
options.diffDiscWalls = true; % walls are discretized by different Nbd
options.putBackDLD = true; % put back into the 1st gap when it reaches xfreeze

% Descriptive run name
prams.runName = runName; 
% This is the name given to all related files 

options.usePlot = ~true; % Plot on-the-fly
options.track = false; % trackers on membrane

% Keep seeding vesicles in the entrance
options.streaming = false;
% Freeze vesicles coming to the end of the device
options.freezing = false;
if options.putBackDLD
  options.freezing = false;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Vesicle
prams.kappa    = 1e-2; % bending coefficient
prams.vesViscCont = 5; % Vesicle viscosity contrast

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% GEOMETRY
scale       = 0.1;
% DLD geometry
% Period of a device
% periodN = floor(1/epsilon);
prams.periodN = (Dposty+Dy)/(tan(theta)*(Dpostx+Dx));

% number of obstacles in x direction  (number of columns)
prams.ncol    = 9; 

% number of obstacles in y direction (number of rows) 
prams.nrow    = 12; 
% So that away from top and bottom walls

% Row-shift fraction
prams.epsilon = tan(theta)*(Dpostx+Dx)/(Dposty+Dy); 

% Sizes
prams.Dpost = Dpost; % diameter of a post in micron * scale
prams.Dy    = Dy; % vertical separation in micron * scale
prams.Dx    = Dx; % horizontal separation in micron * scale

% Vesicle geometry
DeffVes     = 3; % short-axis thickness in microns
scaleVes    = scale*DeffVes/2.23;
prams.scaleVes = scaleVes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Spatial Resolution
prams.N      = Nves; % # of points per vesicle

prams.NbdInt = NbdPill; % # of points on posts
prams.NbdExt = NbdOut; % # of points on exterior wall 
% (1/4 arcLengthSpacing for pillars)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for numerics
prams.gmresTol = 1e-7;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

options.fmm = fmm;  % fmm for single-layer potentials
options.fmmDLP = fmm; % fmm for double-layer potentials
options.matFreeWalls = false; % W2W interactions are done without a matrix
                              % this is slow
                              
options.saveWallMat  = false; % save exact factorization on disk
options.haveWallMats = false; % Do we have wall matrices computed already?

options.outOfCore    = false; % Compute inexact factorization out-of-core
                              % Also block applying exact inverse
% if exact factorization and out-of-core, then we block applying exact inverse
options.memsize      = 1;     % memory usage in gb for out-of-core

options.fastDirect   = false; % Use fast direct solver for preconditioner
options.HODLRforW2W  = false;  % Use HODLR to compress wall2wall interactions
prams.lev_max        = 3;     % Maximum level to go in HODLR
prams.etolFD         = 1e-6;  % tolerance in hodlr to decide l

%if 1
%  !mkdir /workspace/gokberk/DLD_Matrices/DLDN20/
%  options.wallMatFile = '/workspace/gokberk/DLD_Matrices/DLDN20/wallMat';
%else
%  !mkdir ./output/LargeMatDLD/
%  options.wallMatFile = './output/LargeMatDLD/wallMat';
%end
% directory where we save large wall matrices
options.wallMatFile = './output/wallMat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.antiAlias = true; % upsampling for anti-aliasing
options.reparameterization = true; % reparametrization of membrane
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA
options.repulsion = repulsion; % repulsion between vesicles
% Minimum distance (ratio of max. spacing) after which repulsion is active
prams.minDist      = 0.5; 

options.timeAdap = true; % time adaptivity
options.nsdc = 0;  % number of sdc corrections
prams.areaLenTol = rhoAL; 
% tolerance for errors in area-length for time-stepping

if ~options.streaming
  options.saveData = true;  
end
% Name of log file
options.logFile  = ['./output/fullDLDRuns/' prams.runName '.log'];
% Name of binary data file for storing vesicle information
options.dataFile = ['./output/fullDLDRuns/' prams.runName '_Data.bin'];

% Set options and parameters that the user doesn't
[options,prams] = initVes2D(options,prams);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE WALLS
oc = curve;

% Generate the full DLD
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
if ~isempty(Xpost)
% Find where the smallest gap size is, place vesicles there and scale the
% maximum velocity there, Dy is reached at that point
[xpost,ypost] = oc.getXY(Xpost);
% move shape such that rectangle around it is centered at (0,0)
xup = interpft(xpost,512); yup = interpft(ypost,512);
x2 = xpost-(max(xup)-prams.Dpostx/2); x2up = interpft(x2,512);
y2 = ypost-(max(yup)-prams.Dposty/2); y2up = interpft(y2,512);

% look at the points in the horizontal direction and find the largest size
% across the shape in y-direction
maxGsize = -1E+4;
sort_x2 = sort(x2up);
for ij = 1 : numel(x2up)
  ids = find(abs(x2up-sort_x2(ij))<=1e-2);
  ys = y2up(ids);
  gapSize = max(ys)-min(ys);
  
  if gapSize>=maxGsize 
    maxGsize = gapSize;
    [~,idx] = max(y2up(ids));
    maxDist2Left = x2up(ids(idx))+prams.Dpostx/2;  
  end
end
else % circular post
maxDist2Left = prams.Dpostx/2;
end

% Pillar gradient
pillarGrad = tan(theta);

% Zig-zagging ones at the top, displacing ones at the bottom
if isempty(whichCell)
  prams.nv = 2; % e.g. [D;ZZ]
  prams.viscCont = VCs;  % Viscosity contrasts of vesicles
else
  prams.nv = 1;
  prams.viscCont = VCs(whichCell);
end

centx = (-L/2+4*prams.Dx+4*prams.Dpostx)*ones(2,1);
centy = 3*pillarGrad*(prams.Dx+prams.Dpostx) + ...
    [(prams.Dy+prams.Dposty) (-prams.Dy-prams.Dposty)];
if prams.nv ~= 2
  centy = centy(whichCell);
end

% Put back after one column
prams.xfreeze = centx + prams.Dpostx+prams.Dx;
prams.yZigZag = centy - prams.Dy/2-prams.Dposty/2;
prams.nPutBackMax = ceil(prams.periodN*nRepeatPeriod);

angles = 0.05*ones(prams.nv,1);

load relaxed64.mat
Xves = [interpft(X(1:end/2),prams.N);interpft(X(end/2+1:end),prams.N)];
IA = oc.getIncAngle(X);
X = zeros(2*prams.N,prams.nv);

for k = 1 : prams.nv
  X(1:end/2,k) = cos(-IA+angles(k))*Xves(1:end/2)-sin(-IA+angles(k))*Xves(end/2+1:end);
  X(end/2+1:end,k) = sin(-IA+angles(k))*Xves(1:end/2)+cos(-IA+angles(k))*Xves(end/2+1:end);
  
  X(:,k) = [X(1:end/2,k)-mean(X(1:end/2,k))+centx(k);...
      X(end/2+1:end,k)-mean(X(end/2+1:end,k))+centy(k)];
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZE BEFORE RUNNING SIMULATION
if 0
  figure(1);clf;hold on;
  [intWallx,intWally] = oc.getXY(XwallsInt);
  plot([extWallx;extWallx(1)],[extWally;extWally(1)],'k')
  plot([intWallx;intWallx(1,:)],[intWally;intWally(1,:)],'k','linewidth',2)
  plot([X(1:end/2,:);X(1,:)],[X(end/2+1:end,:);X(end/2+1,:)],'r','linewidth',2)
  axis equal
  
  pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET THE TIME STEPPING DETAILS
% Get the length and velocity scale of simulation to decide on params.
[~,~,arcLen] = oc.geomProp(X); 
lenScale = max(arcLen); % length scale of the simulation

% Temporal Resolution
prams.T = 500;   % time horizon
% time scale for simulation
% # of time steps (if constant time stepping) or Dt for the first step
% (adaptive)
prams.m = ceil(prams.T/(0.01*lenScale/options.farFieldSpeed)); 
prams.dtMax = 0.1*lenScale/options.farFieldSpeed; % maximum allowable time-step size
prams.dtMin = 1e-4*lenScale/options.farFieldSpeed; % minimum allowable time-step size  

% first compute velocity at the middle and scale the maximum velocity
gapX = -L/2+[2;3;4;5;6;7;8]*(prams.Dpostx+prams.Dx);
gapY = [1;2;3;4;5;6;7]*pillarGrad*(prams.Dpostx+prams.Dx);
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
Xfinal = Ves2D(X,[],XwallsInt,XwallsExt,prams,options,tt,[],[]);

