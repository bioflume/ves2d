clear; clc;
restoredefaultpath;

fprintf('Single vesicle in a small DLD device.\n');
options.farField = 'DLD'; % background velocity
options.farFieldSpeed = 1.2; % scaling of background velocity
prams.nv = 1; % number of vesicles
options.confined = true; % confined or unbounded geometry
options.diffDiscWalls = true; % walls are discretized by different Nbd

% Descriptive run name
prams.runName = 'st5a1_rateP4'; 
% This is the name given to all related files 

options.usePlot = ~true; % Plot on-the-fly
options.track = false; % trackers on membrane

% Compute the velocity field in device without vesicles
computeVelField   = false;
% Solve for passive tracers only 
tracersNoVesicles = false;
% Keep seeding vesicles in the entrance
options.streaming = true;
% Freeze vesicles coming to the end of the device
options.freezing = true;

if options.streaming
  % 1/scaling for area to check before streaming
  prams.streamRate  = 0.4;
  % Number of vesicles to be streamed
  prams.nSeed = 2;
  % Number of vesicles will be streamed (including the initial ones)
  prams.totnv = 100;
  % Viscosity contrast of the streamed vesicles (2 different)
  % the majority of the vesicles haS VC of prams.vesViscCont
  % the minority has VC of prams.ves2ViscCont
  prams.ves2ViscCont = 1;
  % concentration of the other vesicle; this is the percentage of the
  % total number of vesicles with prams.ves2ViscCont streamed. 
  prams.ves2Concent = 0.9; 
  % VCsXaY means that we want to separate one with VC = X from the other
  % one (so VC=X is the minority)
  !mkdir ./output/streamSim_VCs5a1_rateP4/
  prams.folderName = './output/streamSim_VCs5a1_rateP4/';
  options.saveData = false; % do not save data in the usual way
  options.freezing = true; % we have to freeze if we keep streaming
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Spatial Resolution
prams.N      = 64; % # of points per vesicle

prams.NbdInt = 48; % # of points on posts
prams.NbdExt = 3072; % # of points on exterior wall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% GEOMETRY
scale       = 0.1;
% DLD geometry
% Period of a device
% periodN = floor(1/epsilon);
prams.periodN = 6; 

% number of obstacles in x direction  (number of columns)
prams.ncol = ceil(1.5*prams.periodN);
% number of obstacles in y direction (number of rows)
prams.nrow    = 6; 

% Row-shift fraction
prams.epsilon = 1/prams.periodN; 

% Sizes
prams.Dpost = 15*scale; % diameter of a post in micron * scale
prams.Dy    = 10*scale; % vertical separation in micron * scale
prams.Dx    = 10*scale; % horizontal separation in micron * scale

% Vesicle geometry
DeffVes     = 3; % short-axis thickness in microns
reducedArea = 0.65; % vesicle reduced area
scaleVes    = scale*DeffVes/2.23;
prams.scaleVes = scaleVes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Vesicle
prams.kappa    = 1e-2; % bending coefficient
prams.vesViscCont = 5; % Vesicle viscosity contrast

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for numerics
prams.gmresTol = 1e-7;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

options.fmm = true;  % fmm for single-layer potentials
options.fmmDLP = true; % fmm for double-layer potentials
options.matFreeWalls = false; % W2W interactions are done without a matrix

options.saveWallMat  = false; % save exact factorization on disk
options.haveWallMats = false; % Do we have wall matrices computed already?

options.outOfCore    = false; % Compute inexact factorization out-of-core
                              % Also block applying exact inverse
% if exact factorization and out-of-core, then we block applying exact inverse
options.memsize      = 1;     % memory usage in gb for out-of-core

options.fastDirect   = false;  % Use fast direct solver
options.HODLRforW2W  = true;   % Use HODLR to compress wall2wall interactions
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
prams.minDist      = 0.8; 

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
%load best10shapes.mat
%Xpost = Xposts(:,best10(1));
%prams.Dy = Dys(best10(1));

%Xpost = [interpft(Xpost(1:end/2),prams.NbdInt);...
%  interpft(Xpost(end/2+1:end),prams.NbdInt)];
Xpost = [];

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

% Single vesicle 
if ~options.streaming
  
  prams.nv = 1;
  % start at the vertical location parallel to the 2nd pillar
  centy = -prams.Dy*0.45; 
  centx = xmin+radiusRight+0.6*prams.Dy; 

  angles = 0*ones(prams.nv,1);
  
  irelaxedShape = true; % start with a relaxed vesicle shape

  if ~irelaxedShape
    X = oc.initConfig(prams.N,'nv',prams.nv,'angle',angles,...
       'scale',scaleVes,'center',[centx;centy],'reducedArea',reducedArea);
  else
    load relaxed64.mat
    X = [interpft(X(1:end/2),prams.N);interpft(X(end/2+1:end),prams.N)];
    IA = oc.getIncAngle(X);
    Xnew = zeros(size(X));
    for i = 1 : prams.N
      Xnew(i) = cos(-IA+angles)*X(i)-sin(-IA+angles)*X(i+prams.N);
      Xnew(i+prams.N) = sin(-IA+angles)*X(i)+cos(-IA+angles)*X(i+prams.N);
    end
    X = [Xnew(1:end/2)-mean(Xnew(1:end/2))+centx;Xnew(end/2+1:end)-...
      mean(Xnew(end/2+1:end))+centy];
  end

  % viscosity contrast
  prams.viscCont = prams.vesViscCont*ones(prams.nv,1);  
 
else % streaming 

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

% fprintf('Number of points on exterior wall %2d\n',prams.NbdExt);
% fprintf('Number of interior walls %2d\n',prams.nvbd-1);
% fprintf('Number of unknowns %5d\n',2*prams.NbdExt+2*prams.NbdInt*(prams.nvbd-1)+3*(prams.nvbd-1));
% 
% [reducedArea,area,length] = oc.geomProp(XwallsInt);
% intSpacing = length(1)/prams.NbdInt;
% fprintf('arc-length spacing for interior wall %1.4e\n',length(1)/prams.NbdInt);
% 
% [reducedArea,area,length] = oc.geomProp(XwallsExt);
% extSpacing = length/prams.NbdExt;
% fprintf('arc-length spacing for exterior wall %1.4e\n',length/prams.NbdExt);
% fprintf('Exterior wall should have %2d points \n',3*ceil(length/intSpacing));
% 
% pause

clear extWallx; clear extWally; clear intWallx; clear intWally;
% [~,VesicleArea,~] = oc.geomProp(X);
% Deff = sqrt(4*VesicleArea/pi); % if scale = 0.1, Deff = 0.4 for any RA
% pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


if ~computeVelField
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
  % Run vesicle code  
  Xfinal = Ves2D(X,[],XwallsInt,XwallsExt,prams,options,[],[],[]);
  
else
  prams.nv = 0;
  X = [];  
  
  ntra_x = 1000;
  ntra_y = 1000;
  
  xtra = linspace(-L/2,L/2,ntra_x);
  ytra = linspace(-H/2,H/2,ntra_y);
  [xxtra,yytra] = meshgrid(xtra,ytra);

  % Remove the grid points inside the obstacles
  [intWallx,intWally] = oc.getXY(XwallsInt);
  [extWallx,extWally] = oc.getXY(XwallsExt);
  
  velxTra = zeros(size(xxtra));
  velyTra = zeros(size(yytra));

  Xtra = [xxtra(:);yytra(:)]; 
  om = monitor(X,options,prams);
  tt = tstep(options,prams,om);
  % here, we will remove the grid points outside the exterior domain  
  [vel,allIdx,rmIdxOut] = tt.computeVelFieldNoVes(prams,[],XwallsInt,...
      XwallsExt,Xtra,1);
  
  velx = vel(1:end/2);
  vely = vel(end/2+1:end);

  velxTra(allIdx) = velx;
  velyTra(allIdx) = vely;
  if options.usePlot
    figure(1);clf;hold on;
    plot3([extWallx;extWallx(1)],[extWally;extWally(1)],...
        10*ones(size([extWallx;extWallx(1)])),'k','linewidth',2)
    plot3([intWallx;intWallx(1,:)],[intWally;intWally(1,:)],...
        10*ones(size([intWallx;intWallx(1,:)])),'k','linewidth',2)
    axis equal
  
    plot(xxtra,yytra,'k.','markersize',10);
    plot(xxtra(allIdx),yytra(allIdx),'r.','markersize',10)
    surf(xxtra,yytra,sqrt(velxTra.^2+velyTra.^2))
    
    figure(2);clf; hold on;
    plot([extWallx;extWallx(1)],[extWally;extWally(1)],'k','linewidth',2)
    plot([intWallx;intWallx(1,:)],[intWally;intWally(1,:)],'k','linewidth',2)
    quiver(xxtra,yytra,velxTra,velyTra)
  end
  
  save VelocityFieldN6Stronger rmIdxOut xxtra yytra velxTra velyTra extWallx extWally intWallx intWally allIdx
end
