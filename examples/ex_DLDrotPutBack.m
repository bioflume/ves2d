clear; clc;
restoredefaultpath;

% MAT-FREE setting periodic BCs?
setBCmatFree = false;

options.farField = 'rotDLD'; % background velocity
options.farFieldSpeed = 1.2; % scaling of background velocity
prams.nv = 1; % number of vesicles
options.confined = true; % confined or unbounded geometry
options.diffDiscWalls = true; % walls are discretized by different Nbd
options.putBackDLD = true; % put back into the 1st gap when it reaches xfreeze

% Descriptive run name
prams.runName = 'Fval_VCs5a10p5inCase1'; 
% This is the name given to all related files 

options.usePlot = ~true; % Plot on-the-fly
options.track = false; % trackers on membrane

% Compute the velocity field in device without vesicles
computeVelField   = false;
% Solve for passive tracers only 
tracersNoVesicles = false;
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
prams.periodN = 5; 

% number of obstacles in x direction  (number of columns)
prams.ncol    = 3; 

% number of obstacles in y direction (number of rows) 
prams.nrow    = 4; 
% So that away from top and bottom walls

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
% Spatial Resolution
prams.N      = 32; % # of points per vesicle

prams.NbdInt = 64; % # of points on posts
prams.NbdExt = 1280; % # of points on exterior wall 
% (1/4 arcLengthSpacing for pillars)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for numerics
prams.gmresTol = 1e-7;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

options.fmm = true;  % fmm for single-layer potentials
options.fmmDLP = false; % fmm for double-layer potentials
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
options.repulsion = true; % repulsion between vesicles
% Minimum distance (ratio of max. spacing) after which repulsion is active
prams.minDist      = 0.5; 

options.timeAdap = true; % time adaptivity
options.nsdc = 0;  % number of sdc corrections
prams.areaLenTol = 5e-4; 
% tolerance for errors in area-length for time-stepping

if ~options.streaming
  options.saveData = true;  
end
% Name of log file
options.logFile  = ['./output/' prams.runName '.log'];
% Name of binary data file for storing vesicle information
options.dataFile = ['./output/' prams.runName '_Data.bin'];

% Set options and parameters that the user doesn't
[options,prams] = initVes2D(options,prams);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE WALLS
oc = curve;

% Define the post shape, if empty, it is a circular post shape
if 1
%   load OptimPostsFree.mat
%   load triXpostSharper
  load OptVCs5a10PostsNew
  whichPost = 1;
  [Xpost,Dpostx,Dposty] = oc.samplePost(prams.NbdInt,...
      reshape(bsplinePoints(:,bestIds(whichPost)),8,2));

  %   Xpost = [interpft(Xpost(1:end/2),prams.NbdInt);...
%       interpft(Xpost(end/2+1:end),prams.NbdInt)];
%   Dpostx = max(Xpost(1:end/2))-min(Xpost(1:end/2));
%   Dposty = max(Xpost(end/2+1:end))-min(Xpost(end/2+1:end));

  prams.Dy = 2.5-Dposty;
  prams.Dx = 2.5-Dpostx;
  options.farFieldSpeed = 1.2/prams.Dy;
else
  Xpost = [];
end
% Generate the rotated DLD device (UNIT CELL MODEL)
[XwallsInt,XwallsExt,L,H,radiusRight,radiusUp,...
  prams.Dpostx,prams.Dposty] = oc.initConfigDLDRot(2,prams.NbdInt,...
  prams.NbdExt,Xpost,prams.Dpost,prams.epsilon,prams.Dx,prams.Dy,...
  prams.nrow,prams.ncol);

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

% Zig-zagging ones at the top, displacing ones at the bottom
prams.nv = 2; % e.g. [D;ZZ]
prams.viscCont = [5;10];  % Viscosity contrasts of vesicles

% Pillar gradient
pillarGrad = (prams.Dposty+prams.Dy)*prams.epsilon/(prams.Dpostx+prams.Dx);

ratioToDy = 0.5;
% Top to bottom, left to right
% The one supposed to displace is at the top

% For 4 rows x 3 columns
centy(1) = prams.Dy+prams.Dposty;
centy(2) = 0;

centx(1) = xmin+1/2*prams.Dx+maxDist2Left;
centx(2) = xmin+1/2*prams.Dx+maxDist2Left;

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
 
% Freeze the cells after they move one column
prams.xfreeze = centx(1) + prams.Dpostx+prams.Dx;
prams.yZigZag = [prams.Dy/2+prams.Dposty/2; -prams.Dy/2-prams.Dposty/2];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND THE BOUNDARY CONDITION USING WHOLE DEVICE
fileName = ['./output/' prams.runName '_settingBC.log'];
fid = fopen(fileName,'w');

% 1) Set options and prams for whole DLD device
optionsWhole.farField = 'DLD';
optionsWhole.farFieldSpeed = options.farFieldSpeed;
optionsWhole.confined = true;
optionsWhole.diffDiscWalls = true; % walls are discretized by different Nbd
optionsWhole.antiAlias = true;
pramsWhole.runName = [prams.runName 'Whole'];
optionsWhole.logFile  = ['output/' pramsWhole.runName '.log'];
% Name of binary data file for storing vesicle information
optionsWhole.dataFile = ['output/' pramsWhole.runName '_Data.bin'];
pramsWhole.periodN = prams.periodN;
pramsWhole.ncol = 9;
pramsWhole.nrow = 12;
pramsWhole.epsilon = 1/prams.periodN;
pramsWhole.Dpost = prams.Dpost;
pramsWhole.Dx = prams.Dx;
pramsWhole.Dy = prams.Dy;
pramsWhole.NbdInt = prams.NbdInt;
pramsWhole.NbdExt = 3584;
pramsWhole.gmresTol = prams.gmresTol;  % tolerance for gmres

optionsWhole.fastDirect = false;
optionsWhole.fmm = true;  % fmm for single-layer potentials
optionsWhole.fmmDLP = true; % fmm for double-layer potentials
optionsWhole.matFreeWalls = false; % W2W interactions are done without a matrix
optionsWhole.HODLRforW2W  = false;  % Use HODLR to compress wall2wall interactions
[optionsWhole,pramsWhole] = initVes2D(optionsWhole,pramsWhole);

% 2) Generate the whole DLD device
[XwallsInt_whole,XwallsExt_whole,Lwhole,Hwhole,~,~,~,...
    pramsWhole.Dpostx,pramsWhole.Dposty] = ...
    oc.initConfigDLD('wholeDLD',pramsWhole.NbdInt,pramsWhole.NbdExt,...
    Xpost,pramsWhole.Dpost,pramsWhole.epsilon,pramsWhole.periodN,...
    pramsWhole.Dx,pramsWhole.Dy,pramsWhole.nrow,pramsWhole.ncol);
pramsWhole.nvbdInt = size(XwallsInt_whole,2);
pramsWhole.nvbdExt = 1;
pramsWhole.nvbd = pramsWhole.nvbdInt + 1;

% 3) Move exterior wall for periodic device to the middle of the whole
% device (2 columns from left and right and 3 rows from top and 2 rows from
% bottom w.r.t. to the first column of the whole device)

% center of the first post in the periodic device
xcentper = -L/2+prams.Dx/2+prams.Dpostx/2; 
ycentper = -1.5*(prams.Dy+prams.Dposty);

% this center is supposed to move to
xcentwhole = -Lwhole/2+4*prams.Dx+4.5*prams.Dpostx-maxDist2Left;
ycentwhole = -2.5*(prams.Dposty+prams.Dy) + ...
    pillarGrad*3*(prams.Dpostx+prams.Dx);


% the amount we need to move is 
dispx = xcentwhole-xcentper;
dispy = ycentwhole-ycentper;

% target points (exterior wall of periodic model)
Xtra = [XwallsExt(1:end/2)+dispx;XwallsExt(end/2+1:end)+dispy];

% first compute velocity at the middle and scale the maximum velocity
gapX = -Lwhole/2+[2;3;4;5;6;7;8]*(pramsWhole.Dpostx+pramsWhole.Dx);
gapY = [1;2;3;4;5;6;7]*pillarGrad*(pramsWhole.Dpostx+pramsWhole.Dx);
gapY1 = gapY+prams.Dposty+prams.Dy;
gapY2 = gapY-prams.Dposty-prams.Dy;

% DEBUG
idebug = 0;
if idebug
  figure(1); clf; hold on;
  plot([XwallsExt_whole(1:end/2);XwallsExt_whole(1)],[XwallsExt_whole(end/2+1:end);XwallsExt_whole(end/2+1)],'k','linewidth',2)
  plot([XwallsInt_whole(1:end/2,:);XwallsInt_whole(1,:)],[XwallsInt_whole(end/2+1:end,:);XwallsInt_whole(end/2+1,:)],'k','linewidth',2)
  plot([Xtra(1:end/2);Xtra(1)],[Xtra(end/2+1:end);Xtra(end/2+1)],'r-o','markersize',1,'markerfacecolor','r','linewidth',2)
  plot([XwallsInt(1:end/2,:);XwallsInt(1,:)]+dispx,[XwallsInt(end/2+1:end,:);XwallsInt(end/2+1,:)]+dispy,'r','linewidth',2)
  plot(gapX,gapY,'bo','markersize',6,'markerfacecolor','b')
  plot(gapX,gapY2,'bo','markersize',6,'markerfacecolor','b')
  axis equal
  
  [~,~,len] = oc.geomProp(XwallsInt);
  intSpacing = len(1)/prams.NbdInt;
  fprintf('arc-length spacing for interior wall %1.4e\n',len(1)/prams.NbdInt);
  
  [~,~,len] = oc.geomProp(XwallsExt);
  extSpacing = len/prams.NbdExt;
  fprintf('arc-length spacing for small exterior wall %1.4e\n',len/prams.NbdExt);
  
  [~,~,len] = oc.geomProp(XwallsExt_whole);
  extSpacingWhole = len/pramsWhole.NbdExt;
  fprintf('arc-length spacing for whole exterior wall %1.4e\n',len/pramsWhole.NbdExt);
  
  pause
end

% 4) Construct tstep structure for both DLD models and get BCs
om = monitor(X,options,prams);

tWholeBuild = tic;
ttWhole = tstep(optionsWhole,pramsWhole,om);
tWholeBuild = toc(tWholeBuild);
message = ['Building tstep for whole DLD takes ' num2str(tWholeBuild,'%2.2e') ' seconds'];
fprintf(fid,'%s\n',message);
disp(message)

% find the scaling
tScale = tic;
if ~setBCmatFree
  [vel,~,~] = ttWhole.computeVelFieldNoVes(pramsWhole,[],XwallsInt_whole,...
      XwallsExt_whole,[gapX;gapX;gapX;gapY1;gapY;gapY2],0);
else
  [vel,~,~] = ttWhole.computeVelFieldMatFree(pramsWhole,XwallsInt_whole,...
      XwallsExt_whole,[gapX;gapX;gapX;gapY1;gapY;gapY2],[],0);  
end
tScale = toc(tScale);
message = ['Finding the scaling for farField takes ' num2str(tScale,'%2.2e') ' seconds'];
fprintf(fid,'%s\n',message);
disp(message)
message = ['X-velocity at the gap is ' num2str(mean(vel(1:end/2)),'%2.2e')];
fprintf(fid,'%s\n',message);
disp(message)
message = ['Y-velocity at the gap is ' num2str(mean(vel(end/2+1:end)),'%2.2e')];
fprintf(fid,'%s\n',message);
disp(message)

% Scale the maximum velocity input to match the max. velocity at the gap
scaling = optionsWhole.farFieldSpeed/mean(vel(1:end/2));
optionsWhole.farFieldSpeed = optionsWhole.farFieldSpeed*scaling;
message = ['Scaling is ' num2str(scaling,'%2.2e')];
fprintf(fid,'%s\n',message);
disp(message)

% Set the BC for the whole device again
ttWhole.farField = @(X,Xint) ttWhole.bgFlow(X,optionsWhole.farField,...
    'Speed',optionsWhole.farFieldSpeed,...
    'intWalls',XwallsInt_whole,'nrow',pramsWhole.nrow,...
    'ncol',pramsWhole.ncol,'Dpostx',pramsWhole.Dpostx,...
    'Dposty',pramsWhole.Dposty,'GapX',pramsWhole.Dx,'GapY',pramsWhole.Dy,...
    'epsilon',pramsWhole.epsilon);

% Now compute the velocity on the periodic device's exterior wall
% this is fast
if ~setBCmatFree
  [gPer,~,~] = ttWhole.computeVelFieldNoVes(pramsWhole,[],XwallsInt_whole,...
      XwallsExt_whole,Xtra,0);
else
  [gPer,~,~] = ttWhole.computeVelFieldMatFree(pramsWhole,XwallsInt_whole,...
      XwallsExt_whole,Xtra,scaling,0);  
end

% 5) Set the BCs for periodic device
tt = tstep(options,prams,om); % for unit cell DLD model

tt.farField = @(X,Xint) tt.bgFlow(X,options.farField,...
    'intWalls',XwallsInt,'velG',gPer,'Speed',1);
% Now the boundary conditions are ready.
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


if ~computeVelField
  
  % Run vesicle code  
  Xfinal = Ves2D(X,[],XwallsInt,XwallsExt,prams,options,tt,[],[]);
  
else
  % Compute velocity in domain and check how it looks like in the periodic device  
  % Turn on FMM
  tt.fmm = true;
  tt.fmmDLP = true;
  
  ntra_x = 400;
  ntra_y = 400;
  
  xtra = linspace(-L/2,L/2,ntra_x);
  ytra = linspace(min(extWally),max(extWally),ntra_y);
  [xxtra,yytra] = meshgrid(xtra,ytra);

  % Remove the grid points inside the obstacles
  [intWallx,intWally] = oc.getXY(XwallsInt);

  velxTra = zeros(size(xxtra));
  velyTra = zeros(size(yytra));
  
  Xtra = [xxtra(:);yytra(:)];     
  [vel,allIdx,rmIdxOut] = tt.computeVelFieldNoVes(prams,[],XwallsInt,...
      XwallsExt,Xtra,1);  
  velx = vel(1:end/2);
  vely = vel(end/2+1:end);

  velxTra(allIdx) = velx;
  velyTra(allIdx) = vely;  
  
  velxTra(rmIdxOut) = NaN;
  velyTra(rmIdxOut) = NaN;
  
  % save velocity field and other necessary variables for plotting
  fileName = [prams.runName '_Velocity.mat'];
  save(fileName,'extWallx','extWally','intWallx','intWally','xxtra','yytra',...
  'allIdx','velxTra','velyTra','rmIdxOut','allIdx')
  
end
