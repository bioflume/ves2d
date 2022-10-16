classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.  Handles both implicit and explicit vesicle-vesicle
% interactions, different inextensibility conditions, viscosity
% contrast, solid walls vs. unbounded flows.  This class also
% implements the adaptive time stepping strategy where the errors in
% length, area, and the residual are monitored to adaptively choose a
% time step size.


properties
runName       % descriptive name of the run
om            % Monitor class
Xcoeff        % Coefficients that are used in the discretization of the
              % time derivative when doing IMEX.  For example, if doing
              % IMEX Euler, Xcoeff = [1]
rhsCoeff      % Coefficients that are used in the discretization of the
              % explicit terms of IMEX.  In our case, this is
              % configuration of the arclength and layer potentials are
              % discretized.  For example, if doing IMEX Euler,
              % rhscoeff = [1]
beta          % Term multiplying implicit term in discretization of
              % derivative
order         % Time stepping order
dt            % Time step size
currentTime   % current time needed for adaptive time stepping
finalTime     % time horizon
solver        % method1, method2, method3, or method4.  Changes how
              % inextensiblity condition is treated
Galpert       % Single-layer stokes potential matrix using Alpert
Gbarnett      % layer potential in equation 4.12 of Alex and 
              % Shravan's paper
D             % Double-layer stokes potential matrix
lapDLP        % Double-layer laplace potential matrix
DLPnoCorr     % Double-layer stokes potential matrix without correction
SLPnoCorr     % Single-layer stokes potential matrix without correction
SLPnoCorrRegul% regularized with diag = sum of the entries in the row

wallRestrict  % Rate that the preconditioner is restricted for
              % the solid wall components

wallDLP       % Double-layer potential due to solid walls
wallN0        % Modificiation of double-layer potential to 
              % remove the rank deficiency on outer boundary
              
wallDLPint    % Double-layer potential due to interior solid walls
wallDLPext    % Double-layer potential due to exterior solid walls    

wallDLPnoCorr % DLP without corrections for the singular terms
wallDLPintNoCorr % DLP without corrections for the singular terms
wallDLPextNoCorr % DLP without corrections for the singular terms

farField      % Far field boundary condition
confined      % whether geometry is bounded or not
diffDiscWalls % whether there are two sets of walls disc. w/ diff Ns (e.g. DLD)

bdiagVes      % precomputed inverse of block-diagonal precondtioner
              % only for vesicle-vesicle interactions
bdiagTen
bdiagWall     % precomputed inverse of block-diagonal precondtioner
              % only for wall-wall interactions
              
bdiagWallFD   % precomputed inverse of b-d preconditioner using fast direct
              % solver for DLP and a Schur complement
vesves        % Discretization of vesicle-vesicle interactions.
              % Can be implicit or explicit
fmm           % with or without the FMM
fmmDLP        % use FMM for the double-layer potential
profile       % display profile times throughout simulation
bending       % bending is defined so that the vesicle tends to
              % a prescribed configuration
gmresTol      % GMRES tolerance
gmresMaxIter  % maximum number of gmres iterations
timeAdap      % Using adaptive time stepping
orderGL       % Order of Gauss-Lobatto points
GLpts         % Gauss-Lobatto points 
SDCcorrect    % timeStep changes depending on if we are doing
              % an sdc correction or forming a provisional
              % solution
areaLenTol    % maximum allowable error in area and length
nsdc          % number of sdc iterations to take
NearV2V       % near-singular integration for vesicle to vesicle
NearW2V       % near-singular integration for wall to vesicle
NearV2W       % near-singular integration for vesicle to wall 
NearW2W       % near-singular integration for wall to wall 

NearV2Wint    % near-singular structure for vesicle to inner walls
NearWint2V    % near-singular structure for inner walls to vesicle
NearV2Wext    % near-singular structure for vesicle to external wall
NearWext2V    % near-singular structure for external wall to vesicle
NearWint2Wext % near-singular structure for internal wall to external wall
NearWint2Wint % near-singular structure for internal wall to internal wall
NearWext2Wint % near-singular structure for external wall to internal wall

op            % class for evaluating potential so that we don't have
              % to keep building quadrature matricies
opWall        % class for walls              

opWallInt     % class for interior walls in DLD
opWallExt     % class for exterior wall in DLD

adhesion      % use adhesion in the model
periodic      % fudges in periodicity into the geometry
adStrength    % strength of adhesion
adRange       % range of adhesion
antiAlias     % use anti-aliasing

fmmPrecision  % precision of fmm

repulsion     % use repulsion in the model
repStrength   % repulsion strength
minDist       % minimum distance to activate repulsion

betaUp        % maximum amount time step is scaled up
betaDown      % maximum amount time step is scaled down
alpha         % safety factor scaling for new time step size

dtMax         % maximum time step size
dtMin         % minimum time step size
betaInc       % tolerance allowing us to use derivatives to increase dt
betaDec       % tolerance allowing us to use derivatives to decrease dt

matFreeWalls  % Compute wall-wall interactions matrix free
HODLRforW2W   % whether we use HODLR to compress wall2wall interactions (not fact. of preconditioner)
wallDLPandRSmat % wall2wall interaction matrix computed in initialConfined

% following three are needed when HODLR is used to compress M11
wallDLPandRSmatM12 
wallDLPandRSmatM21
wallDLPandRSmatM22

haveWallMats  % do we have wall matrices computed before?
fastDirect    % use fast direct solver to factorize preconditioner for walls
saveWallMat   % save or do not save the wall matrices

outOfCore     % flag for out-of-core matrix-vector multiplications

wallMatFile   % if we save the wall matrices (also hodlr matrices), fileName
maxbsize      % maximum block size for performing mat-vector multiplications
hodlrWalls    % hodlr for wall2wall interactions
hodlr         % hodlr class keeping all necessary matrices needed for preco.

alsoExplicit % if explicit solve is asked to be done
usePreco     % use block-diagonal preco?
useSpecPreco % use spectral preconditioner?

randVesicles % used for optimization only
Xshapes % random shapes
sigShapes % tensions of these cells
Xrand  % random vesicles in simulation

% These might be saved when we are scaling RHS for Wall2Wall interactions
% and hence do not compute again, just scale it
etaExt
etaInt
RS

% Inverse of the blocks of wall-2-wall interaction matrix
invM11Ext
invM11Int
invM22

saveVinf
saveVtotal
randFreqs
randAmps

verbose
end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(options,prams,om)
% o=tstep(options,prams,monitor): constructor.  Initialize class.  
% Take all elements of options and prams needed by the time stepper, and 
% monitor class.

o.runName = prams.runName; % descriptive name of the run (for files)

o.om = om; % monitor class
o.order = options.order; % Time stepping order
o.dt = prams.T/prams.m; % Time step size

% Get coefficients for doing time integration
[o.Xcoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);

% Method always starts at time 0
o.currentTime = 0;

% Need the time horizon for adaptive time stepping
o.finalTime = prams.T;

o.solver = options.inextens; % Discretization of inextensibility
o.vesves = options.vesves; % Vesicle-vesicle interactions

% fast-multipole method to compute interactions
o.fmm = options.fmm; % fmm for computing single-layer potentials
o.fmmDLP = options.fmmDLP; % fmm for computing double-layer potentials

% disp profile times
o.profile = options.profile;
% user-defined prefered configuration of vesicle
o.bending = options.bending;
% GMRES tolerance
o.gmresTol = prams.gmresTol;
% maximum number of GMRES iterations
o.gmresMaxIter = prams.gmresMaxIter;

% Far field boundary condition built as a function handle
o.farField = @(X,Xint) o.bgFlow(X,options.farField,...
    'Speed',options.farFieldSpeed,'vortexSize',prams.vortexSize,...
    'intWalls',Xint,'nrow',prams.nrow,...
    'ncol',prams.ncol,'Dpostx',prams.Dpostx,'Dposty',prams.Dposty,...
    'GapX',prams.Dx,'GapY',prams.Dy,'epsilon',...
    prams.epsilon);

% Confined or unbounded geometry
o.confined = options.confined;
% if there are walls discretized with different number of points
o.diffDiscWalls = options.diffDiscWalls;

% Time adaptivity flag
o.timeAdap = options.timeAdap;
% Gauss-Lobatto order
o.orderGL = options.orderGL;
% load Gauss-Lobatto points
o.GLpts = o.gaussLobatto(o.orderGL);

% if the shape is being corrected after each time step, then the
% areaLenTol is exactly the tolerance for the error
% accumulated at each time step.  
o.areaLenTol = prams.areaLenTol;

% number of sdc iterations to take
o.nsdc = options.nsdc;
if o.order > 1
  o.SDCcorrect = false;
end

% parameters for adaptive time stepping
o.betaUp = prams.betaUp; % constant factor to increase time step size
o.betaDown = prams.betaDown; % constant factor to decrease time step size
o.alpha = prams.alpha; % safety factor for increasing time step size
o.dtMax = prams.dtMax; % maximum allowable time step size
o.dtMin = prams.dtMin; % minimum allowable time step size
o.betaInc = prams.betaInc; % tolerance checking asymptotic assumption
o.betaDec = prams.betaDec; % tolerance checking asymptotic assumption
 
% These two are not used in the current code
o.periodic = options.periodic;
o.wallRestrict = 1;

% strength of adhesion (not used in the current code)
o.adhesion = options.adhesion;
o.adStrength = prams.adStrength;
% range of adhesion
o.adRange = prams.adRange;

% Repulsion between vesicles and vesicles-walls.
% if there is only one vesicle, turn off repulsion
o.repulsion = options.repulsion;
if prams.nv == 1 && ~options.confined
  o.repulsion = false;
end
% scaling the strength of the repulsion
o.repStrength = prams.repStrength;
% scaling the range of repulsion
o.minDist = prams.minDist;

% this is used for optimization (random vesicles are put around vesicle and
% forces are computed with respect to them)
o.randVesicles = options.randVesicles;
% if empty, no random vesicles to be added
% if not, then it is either addBelow or addAbove, so we add random vesicles
if ~isempty(o.randVesicles) 
  load vesShapesFromDense.mat
   % make sure we have same number of points
  o.Xshapes = [interpft(Xshapes(1:end/2,:),prams.N);...
      interpft(Xshapes(end/2+1:end,:),prams.N)];
  o.sigShapes = sigShapes;
  o.Xrand = [];
end
% anti-aliasing flag
o.antiAlias = options.antiAlias;

% fmm precision (see poten.m for more details)
o.fmmPrecision = options.fmmPrecision;
% build poten class for vesicles
o.op = poten(prams.N,o.fmmPrecision,o.antiAlias);

% explicit solve is also done?
o.alsoExplicit = options.alsoExplicit;

% use preco?
o.usePreco = options.usePreco;
o.useSpecPreco = false;
o.bdiagVes = [];
% out-of-core (wall2wall interactions and preco.)
o.outOfCore = options.outOfCore;

% for random velocities in shear flow
o.randFreqs = [];
o.randAmps = [];
o.saveVinf = options.saveVinf;
if o.saveVinf
  fid = fopen(['./output/' o.runName '_vInf.bin'],'w');
  fclose(fid);
end

o.saveVtotal = options.saveVtotal;
if o.saveVtotal 
  %fid = fopen(['./output/' o.runName '_vTotal.bin'],'w');
  %fclose(fid);
  %fid = fopen(['./output/' o.runName '_vTotalDecoup.bin'],'w');
  %fclose(fid);
  fid = fopen(['./output/freeCouetteShortRuns/vTotalFiles/' ...
      o.runName '_vTotalPreComp.bin'],'w');
  fclose(fid);
end
   
o.verbose = options.verbose;

% build poten classes for walls
if options.confined
  % initialize density and RS for computing w2w interactions only  
  o.etaExt = [];
  o.etaInt = [];
  o.RS = [];
  
  % flag for computing the W2W interactions with a precomp. matrix or not
  % if matFreeWalls = true, then we compute W2W interactions at every time
  % step either using FMM or not. If matFreeWalls = false, then we use the
  % precomputed matrix and apply it to a vector of density and RS to get
  % the wall2wall interactions. If fastDirect = true, then the wall2wall
  % interaction matrix is compressed.
  o.matFreeWalls = options.matFreeWalls;  
  
  % if there is one wall, it does not speed up the computation much
  if prams.nvbd == 1
    o.matFreeWalls = true;
  end
  
  % flag for using fast direct solver for walls
  o.fastDirect = options.fastDirect;  
   
  % Use HODLR to compress wall2wall interactions?
  o.HODLRforW2W = options.HODLRforW2W;
  
   % if wall2wall is handled matrix free, then there is nothing to compress
  if o.matFreeWalls
    o.HODLRforW2W = false;
  end
  
  % if wall2wall is compressed with HODLR, then matFreeWalls is false.
  if o.HODLRforW2W
    o.matFreeWalls = false;
  end
  
  % do we have wall matrices computed before
  o.haveWallMats = options.haveWallMats;  
  % save or do not save wall matrices
  o.saveWallMat = options.saveWallMat;
  if o.verbose 
  message = ['Building poten classes for walls...'];
  o.om.writeMessage(message,'%s\n');
  tPot0 = tic;
  end
  % if all walls are discretized with the same Nbd 
  if ~options.diffDiscWalls
    o.opWall = poten(prams.Nbd,o.fmmPrecision,o.antiAlias);
    o.opWallInt = [];
    o.opWallExt = [];
  else % if there are two sets of walls discretized separately
    o.opWallInt = poten(prams.NbdInt,o.fmmPrecision,o.antiAlias);
    o.opWallExt = poten(prams.NbdExt,o.fmmPrecision,o.antiAlias);
    o.opWall = [];
  end
  if o.verbose
  message = ['DONE, it took ' num2str(toc(tPot0),'%2.2e') ' seconds'];
  o.om.writeMessage(message,'%s\n');
  o.om.writeMessage(' ','%s\n');
  end
  % filename to save the matrices necessary for HODLR (or if the matrices
  % are already computed and saved, use them)
  o.wallMatFile = options.wallMatFile;

  % if walls are discretized with the same number of points, we do not have
  % hodlr implemented for that.
  if ~o.diffDiscWalls
    o.fastDirect = false;
    o.HODLRforW2W = false;
  end
  
  % if we are to use fastDirect, then construct hodlr class
  if o.fastDirect
    % maximum block size (maxbsize) is computed if we block the matrix 
    % operations
    memsize = options.memsize;
    o.maxbsize = floor(1e9*memsize/(8*(prams.NbdExt*2+...
      prams.NbdInt*2*(prams.nvbd-1))));
  
    % log-file to keep the details of HODLR
    hodlrLogFile = ['./output/' prams.runName '_HODLR.log'];
    
    % hodlr class for compressing and factorizing wall2wall interactions
    o.hodlr = hodlr(prams.lev_max,prams.etolFD,memsize,o.outOfCore,1,...
        hodlrLogFile);
    message = 'HODLR is in use for compressing and factorizing wall2wall ints.';
    o.om.writeMessage(message,'%s\n');
    o.om.writeMessage(' ','%s\n');  
  end
  
  % if we want to compress wall2wall interactions, then construct hodlr
  % class, this is different than fastDirect, it factorizes also.
  if o.HODLRforW2W 
    memsize = options.memsize;
    o.maxbsize = floor(1e9*memsize/(8*(prams.NbdExt*2+...
      prams.NbdInt*2*(prams.nvbd-1))));
  
    % log-file to keep the details of HODLR
    hodlrLogFile = ['./output/' prams.runName '_HODLR.log'];
    
    % hodlr class for compressing wall2wall interactions
    o.hodlrWalls = hodlr(prams.lev_max, prams.etolFD, memsize, ...
        o.outOfCore,1,hodlrLogFile);    
    message = 'HODLR is in use only for compressing wall2wall ints.';
    o.om.writeMessage(message,'%s\n');
    o.om.writeMessage(' ','%s\n');  
  end
  
else
  o.opWall = [];
  o.opWallInt = [];
  o.opWallExt = [];
  o.fastDirect = false;
end


end % tstep: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xcoeff,rhsCoeff,beta] = getCoeff(o,order)
% [Xcoeff,rhsCoeff,beta] = getCoeff(order) generates the coefficients
% required to discretize the derivative.  First-order time  derivatives
% are approximated by beta*x^{N+1} + rhscoeff.*[x^{N} x^{N-1} ...]
% Explicit terms (operators) are discretized at Xcoeff.*[x^{N} x^{N-1}
% ...] All rules are from Ascher, Ruuth, Wetton 1995.

if (order == 4) % fourth-order
  beta = 25/12;
  Xcoeff = [-1; 4; -6; 4]; 
  rhsCoeff = [-1/4; 4/3; -3; 4];
elseif (order == 3) % third-order
  beta = 11/6;
  Xcoeff = [1; -3; 3];
  rhsCoeff = [1/3; -3/2; 3];
elseif (order == 2) % second-order
  beta = 1.5;
  Xcoeff = [-1; 2];
  rhsCoeff = [-0.5; 2];
else % first-order
  beta = 1;
  Xcoeff = 1;
  rhsCoeff = 1;
end

end % getCoeff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [walls,wallsInt,wallsExt] = initialConfined(o,prams,Xwalls,...
        XwallsInt,XwallsExt)
% [walls,wallsInt,wallsExt] = 
% initialConfined(o,prams,Xwalls,XwallsInt,XwallsExt) builds 
% an object of class capsules for the solid walls.

% if the walls are discretized with the same Nbd
if ~o.diffDiscWalls
  
  potWall = o.opWall; % poten class for walls
  % velocity on solid walls coming from no-slip boundary condition
  [uwalls,~] = o.farField(Xwalls,[]);

  % build the walls
  walls = capsules(Xwalls,[],uwalls,zeros(prams.nvbd,1),...
      zeros(prams.nvbd,1),o.antiAlias);
  % set the upsampling rate
  if o.antiAlias
    walls.setUpRate(potWall);
  end
  % build the double-layer potential matrix for walls and save on memory
  o.wallDLP = potWall.stokesDLmatrix(walls);
  % build the double-layer matrix for walls without any correction if FMM
  % is used
  if o.fmmDLP && o.matFreeWalls
    o.wallDLPnoCorr = potWall.stokesDLmatrixNoCorr(walls);
  else
    o.wallDLPnoCorr = [];
  end
  
  % N0 to remove rank-1 deficiency
  o.wallN0 = potWall.stokesN0matrix(walls);
  
  % block diagonal preconditioner for solid walls
  o.bdiagWall = o.wallsPrecond(walls);

  wallsInt = [];
  wallsExt = [];
else % if there are two sets of walls discretized with different Nbds
  
  % poten classes for walls
  potWallInt = o.opWallInt;
  potWallExt = o.opWallExt;
  % velocity on solid walls coming from no slip boundary condition
  [uwallsExt,uwallsInt] = o.farField(XwallsExt,XwallsInt);
  
  % build the walls
  wallsInt = capsules(XwallsInt,[],uwallsInt,...
      zeros(prams.nvbdInt,1),zeros(prams.nvbdInt,1),o.antiAlias);
  wallsExt = capsules(XwallsExt,[],uwallsExt,0,0,o.antiAlias);
  
  % set the upsampling rate for walls
  if o.antiAlias
    wallsInt.setUpRate(potWallInt);
    wallsExt.setUpRate(potWallExt);
  end
    
  % build the double layer potential matrix for walls and save on memory
  if isempty(o.wallDLPint)
    o.wallDLPint = potWallInt.stokesDLmatrix(wallsInt);
  end
  if isempty(o.wallDLPext)
    o.wallDLPext = potWallExt.stokesDLmatrix(wallsExt);
  end
  
  % build the DLP for walls without any correction if FMM is on
  if o.fmmDLP && o.matFreeWalls
    if isempty(o.wallDLPintNoCorr)
      o.wallDLPintNoCorr = potWallInt.stokesDLmatrixNoCorr(wallsInt);
    end
    if isempty(o.wallDLPextNoCorr)
      o.wallDLPextNoCorr = potWallExt.stokesDLmatrixNoCorr(wallsExt);    
    end
  else
    o.wallDLPintNoCorr = [];
    o.wallDLPextNoCorr = [];
  end
  
  % N0 is computed only for the exterior wall
  if isempty(o.wallN0)
    o.wallN0 = potWallExt.stokesN0matrix(wallsExt);
  end
  
  % block diagonal preconditioner for solid walls
  if ~o.fastDirect % inverse
    if isempty(o.bdiagWall)  
      o.bdiagWall = o.wallsPrecondDiffDisc(wallsInt,wallsExt);
    end
    % if do not use fastDirect but want to compress wall2wall, then do it
    % here
    if isempty(o.wallDLPandRSmatM12) 
      if o.HODLRforW2W
        o.compressWall2Wall(wallsInt,wallsExt);
      end
    end
  else % using fast direct solver 
    if isempty(o.bdiagWallFD)
      [~] = o.wallsPrecondDiffDisc(wallsInt,wallsExt);
      o.bdiagWallFD = @(rhs) o.applyWallsPrecondFD(prams,rhs);
    end
    % if fastDirect, and also want to compress wall2wall, then do it inside
    % fastDirect because matrix is already built there.
  end
  
  walls = [];
end

end % initialConfined


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xstore,sigStore,uStore,etaStore,RSstore,Xtra] = ...
  firstSteps(o,options,prams,Xinit,walls,Xtra,pressTar)

% [Xstore,sigStore,uStore,etaStore,RSstore] = ...
%   firstSteps(options,prams,Xinit,sigInit,uInit,...
%   walls,om,pressTar)
% refines the first time step [0,dt] and uses a first-order, then
% second-order, then third-order, etc time stepping to find the vesicle
% position, tension, velocity, and density function eta defined on the
% solid walls and the rotlets and stokeslets at t = dt.  Returns
% Xstore, sigStore, uStore etaStore, so that they are immediatly ready
% to use for higher-order time stepping.  This routine also computes
% the tension and density functions of the initial condition.  This is
% needed in SDC since we require these functions at all Gauss-Lobatto
% points in order to compute the residual

% Higher order time stepping and nsdc>=1 are not implemented for 
% cases where the walls are discretized with different Nbds 
% (for which we use firstStepsDLD).

om = o.om; % load the monitor class

N = size(Xinit,1)/2; % points per vesicle
nv = size(Xinit,2); % number of vesicles
if o.confined 
  Xwalls = walls.X; % discretization points of solid walls
  Nbd = size(Xwalls,1)/2; % points per wall
  nvbd = size(Xwalls,2); % number of wall components
else
  Xwalls = [];
  Nbd = 0;
  nvbd = 0;
end

% time horizon is enough time steps so that we can
% continue with higher-order time integrator
prams.T = prams.T/prams.m*(o.order-1);

% For second-order time stepping, this keeps the local error from time
% step 0 to dt on the order of dt^3
if o.order ~=1
  mR = min(ceil(prams.m/32),100)*o.order^2;
  mR = mR * (o.order - 1);
else
  mR = 1;
end

% number of time steps to take in range [0,dt]
prams.m = mR*(o.order-1);


Xstore = zeros(2*N,nv,prams.m+1);
sigStore = zeros(N,nv,prams.m+1);
uStore = zeros(2*N,nv,prams.m+1);
etaStore = zeros(2*Nbd,nvbd,prams.m+1);
RSstore = zeros(3,nvbd,prams.m+1);

Xstore(:,:,1) = Xinit;
vesicle = capsules(Xinit,zeros(N,nv),[],...
    prams.kappa,prams.viscCont,o.antiAlias);
if o.antiAlias
  vesicle.setUpRate(o.op);
end

% need intial tension, density function, rotlets, and stokeslets so
% that we can do SDC updates, couling with advection-diffusion system,
% etc
[sigStore(:,:,1),etaStore(:,:,1),~,~,RSstore(:,:,1),u,iter] = ...
    vesicle.computeSigAndEta(o,walls,[],[]);

% delete previous versions of files and write some initial
% options/parameters to files and the console
om.initializeFiles(Xinit,sigStore(:,:,1),Xwalls,[],[],Xtra,pressTar);

message = ['GMRES took ',num2str(iter),...
    ' iterations to find intial tension and density function'];
om.writeMessage(message,'%s\n');
om.writeMessage(' ','%s\n');

for n = 1:prams.m
  time = n*prams.T/prams.m;

  % take the highester order possible
  options.order = min(n,o.order);
  tt = tstep(options,prams);
  
  % no SDC corrections
  tt.SDCcorrect = false;
  
  tt.wallDLP = o.wallDLP;
  tt.wallN0 = o.wallN0;
  
  % build inverse of the wall-to-wall intearctions includes diagonal and
  % off-diagonal terms.  That is, it could be used to solve the stokes
  % equations in a multiply-connected domain with no vesicles
  if ~o.fastDirect
    tt.bdiagWall = o.bdiagWall;
  else
   tt.bdiagWallFD = o.bdiagWallFD;
  end
 
  % take one time step
  updatePreco = true;
  [X,sigma,u,eta,RS,iter,iflag] = tt.timeStep(...
      Xstore(:,:,n-tt.order+1:n),...
      sigStore(:,:,n-tt.order+1:n),...
      uStore(:,:,n-tt.order+1:n),...
      etaStore(:,:,n-tt.order+1:n),...
      RSstore(:,:,n-tt.order+1:n),...
      [],[],[],[],[],[],[],...
      prams.kappa,prams.viscCont,walls,...
      updatePreco,vesicle);
  

  % compute the next position of tracers
  if numel(Xtra) > 1
    vel = o.tracersVel(X,sigma,u,...
        prams.kappa,prams.viscCont,walls,eta,RS,Xtra);
    Xtra = Xtra + tt.dt*vel;
  end

  % Required for adaptive time stepping which have not been
  % implemented yet for high-order time stepping
  accept = true;
  dtScale = 0;
  res = 0;
  
  % save and update position, tension, velocity, and
  % density function
  Xstore(:,:,n+1) = X;
  sigStore(:,:,n+1) = sigma;
  uStore(:,:,n+1) = u;
  etaStore(:,:,n+1) = eta;
  RSstore(:,:,n+1) = RS;
  

  % save data, write to log file, write to console as
  % requested
  terminate = om.outputInfo(X,sigma,u,eta,RS,...
      Xwalls,Xtra,time,iter,dtScale,res,1,[],iflag);
end 
% end of using small time steps to get the simulation far enough
% to use the desired order time integrator

% Only want every mR^{th} solution as these correspond to
% time steps that are multiples of dt
% Only need to do this if the time stepping order is not 1
if o.order > 1
  if options.pressure
    op = poten(N,o.fmmPrecision);
    % compute the pressure and stress due to the vesicles
    % and the solid walls
    [press,stress1,stress2] = op.pressAndStress(...
        X,sigma,u,prams.kappa,prams.viscCont,...
        walls,pressTar,eta,RS,options.confined,...
        options.fmm);
    
    if options.saveData
      om.writePressure(press);
      om.writeStress(stress1(1:end/2),'11');
      om.writeStress(stress1(1+end/2:end),'12');
      om.writeStress(stress2(1:end/2),'21');
      om.writeStress(stress2(1+end/2:end),'22');
    end
  end % ~options.pressure
  % write pressure and stress to output file at time dt

  Xstore = Xstore(:,:,1:mR:end);
  sigStore = sigStore(:,:,1:mR:end);
  uStore = uStore(:,:,1:mR:end);
  etaStore = etaStore(:,:,1:mR:end);
  RSstore = RSstore(:,:,1:mR:end);
end

end % firstSteps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xall,sigAll,uAll,etaInt,etaExt,RS,Xtra] = ...
        firstStepsDLD(o,prams,Xall,sigAll,uAll,...
        wallsInt,wallsExt,rem_ids,Xtra,istreaming,allVesViscConts)
% [Xall,sigAll,uAll,etaInt,etaExt,RS,Xtra] = ...
%   firstStepsDLD(o,prams,Xall,sigAll,uAll,...
%   wallsInt,wallsExt,rem_ids,Xtra)
% refines the first time step [0,dt] and uses only first-order
% time stepping to find the vesicle
% position, tension, velocity, and density function eta defined on the
% solid walls and the rotlets and stokeslets at t = dt.  Returns
% Xall, sigAll, uAll, etaInt, etaExt,, so that they are immediatly ready
% to use for higher-order time stepping. 

% Unlike firstSteps (only for one set of walls disc. w/ the same Nbds) 
% above, this routine works only for nsdc = 0 and first
% order time stepping. nsdc>= 1 is not implemented for the cases where the 
% walls are discretized separately.
om = o.om;

% If we are seeding vesicles, reach only the active ones in the first step
Xinit = Xall(:,rem_ids);

N = size(Xinit,1)/2; % points per vesicle
nv = size(Xinit,2); % number of vesicles

XwallsInt = wallsInt.X; % points on the interior walls
XwallsExt = wallsExt.X; % points on the exterior walls

% Compute intial tension, density function, rotlets, and stokeslets
vesicle = capsules(Xinit,zeros(N,nv),[],...
    prams.kappa,prams.viscCont,o.antiAlias);
if o.antiAlias
  vesicle.setUpRate(o.op);
end

[sigInit,~,etaInt,etaExt,RS,uInit,iter] = ...
    vesicle.computeSigAndEta(o,[],wallsInt,wallsExt);

uAll(:,rem_ids) = uInit;
sigAll(:,rem_ids) = sigInit;


% delete previous versions of files and write some initial
% options/parameters to files and the console
om.initializeFiles(Xall,sigAll,[],XwallsInt,XwallsExt,Xtra,[]);

if istreaming
% if we are streaming, we save at every accepted time step
% initially we save structures and viscosity contrasts of vesicles and
% other details
fileName = [prams.folderName prams.runName '_DataInst_0.mat'];
save(fileName,'XwallsInt','XwallsExt','Xinit','rem_ids','allVesViscConts',...
    'prams')
end

message = ['GMRES took ',num2str(iter),...
    ' iterations to find intial tension and density function'];
om.writeMessage(message,'%s\n');
om.writeMessage(' ','%s\n');


end % firstStepsDLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,sigma,u,eta,etaInt,etaExt,RS,iter,accept,dtScale,normRes,...
    iflag,collUprate] = timeStepGL(o,Xstore,sigStore,uStore,...
    etaStore,etaIntStore,etaExtStore,RSstore,kappa,viscCont,...
    walls,wallsInt,wallsExt)
% [X,sigma,u,eta,etaInt,etaExt,RS,iter,accept,dtScale,normRes,...
%     iflag,collUprate] = timeStepGL(o,Xstore,sigStore,uStore,...
%     etaStore,etaIntStore,etaExtStore,RSstore,kappa,viscCont,...
%     walls,wallsInt,wallsExt)
% takes the desired number of time steps at Gauss-Lobatto points to
% find solution at dt given solution at 0.  Calls o.timeStep which is
% what we use for the constant sized time step integration.  Returns a
% flag that says if the solution was accepted or not and the amount the
% time step was scaled.  errors is the norm of the residual and iflag
% is tells if any of the gmres runs failed.

om = o.om; % load the monitor class

N = size(Xstore,1)/2; % number of points per vesicle
nv = size(Xstore,2); % number of vesicles

% Get the discretization information about walls if confined
Nbd = size(etaStore,1)/2; % number of points per wall (if all disc. w/ Nbd)
NbdInt = size(etaIntStore,1)/2; % number of points per interior wall
NbdExt = size(etaExtStore,1)/2; % number of points per exterior wall
nvbdInt = size(etaIntStore,2);
nvbdExt = size(etaExtStore,2);
nvbd = size(etaStore,2)+nvbdInt+nvbdExt;

% initial areas and lengths of vesicles (these must be conserved)
a0 = om.area; 
l0 = om.length;
% build the curve class
oc = curve;

% Need to save the position, tension, velocity, density
% functions, rotlets/stokeslets at the Gauss-Lobatto
% time steps

% for vesicles
Xprov = zeros(2*N,nv,abs(o.orderGL));
sigmaProv = zeros(N,nv,abs(o.orderGL));
uProv = zeros(2*N,nv,abs(o.orderGL));
% for walls
etaProv = zeros(2*Nbd,size(etaStore,2),abs(o.orderGL));
etaIntProv = zeros(2*NbdInt,nvbdInt,abs(o.orderGL));
etaExtProv = zeros(2*NbdExt,nvbdExt,abs(o.orderGL));
RSprov = zeros(3,nvbd,abs(o.orderGL));

% Store the initial conditions in the provisional
% solution which will store the solution at all
% Gauss-Lobatto points
for k = 1:nv
  Xprov(:,k,1) = Xstore(:,k);
  sigmaProv(:,k,1) = sigStore(:,k);
  uProv(:,k,1) = uStore(:,k);
end

if ~o.diffDiscWalls
  for k = 1:nvbd
    etaProv(:,k,1) = etaStore(:,k);
    RSprov(:,k,1) = RSstore(:,k);
  end
else
  etaExtProv(:,1,1) = etaExtStore(:);
  RSprov(:,1,1) = RSstore(:,1);
  for k = 1 : nvbdInt
    etaIntProv(:,k,1) = etaIntStore(:,k);
    RSprov(:,k+1,1) = RSstore(:,k+1);
  end
end

% need the single-layer potential at all levels
% of the provisional solution
Galpert = zeros(2*N,2*N,nv,abs(o.orderGL));
Gbarnett = [];

% need the double-layer potential at all levels
% of the provisional solution
if any(viscCont ~= 1)
  D = zeros(2*N,2*N,nv,abs(o.orderGL));
else
  D = [];
end

% Errors that are solved for in each iteration.  The first column of
% deltaX will always be zero since we assume the inital condition is
% exact
deltaX = zeros(2*N,nv,abs(o.orderGL));
deltaSigma = zeros(N,nv,abs(o.orderGL));
deltaEta = zeros(2*Nbd,size(etaStore,2),abs(o.orderGL));
deltaEtaInt = zeros(2*NbdInt,nvbdInt,abs(o.orderGL));
deltaEtaExt = zeros(2*NbdExt,nvbdExt,abs(o.orderGL));
deltaRS = zeros(3,nvbd,abs(o.orderGL));


X = Xprov(:,:,1);
sigma = sigmaProv(:,:,1);
u = uProv(:,:,1);

% vesicle configuration at first Gauss-Lobatto point
vesicle(1) = capsules(X,sigma,u,kappa,viscCont,o.antiAlias);
% compute required upsampling rate
if o.antiAlias
  vesicle(1).setUpRate(o.op);
end

% compute intial reduced area, area, and length if using adaptive
% time stepping
if o.timeAdap
  Xinit_tadap = X;
  [~,aInit,lInit] = oc.geomProp(X);
end

% need to save the time step size
dt = o.dt;

% time step sizes use for Gauss-Lobatto points
dtimes = diff(o.GLpts)*dt/2;

% START OF FORMING PROVISIONAL SOLUTION
iflag = 0; % initialize gmres flag as everything okay
iter = 0; % initialize number of gmres iterations
updatePreco = true; % Form preconditioner at initial configuartion

for k = 1:abs(o.orderGL)-1
  o.dt = dtimes(k);
  o.SDCcorrect = false;
  
  % form provisional solution at next Gauss-Lobatto point
  [X,sigma,u,eta,etaInt,etaExt,RS,subIter,iflagTemp] = o.timeStep(...
      Xprov(:,:,k-o.order+1:k),sigmaProv(:,:,k-o.order+1:k),...
      uProv(:,:,k-o.order+1:k),etaProv(:,:,k-o.order+1:k),...
      etaIntProv(:,:,k-o.order+1:k),etaExtProv(:,:,k-o.order+1:k),...
      RSprov(:,:,k-o.order+1:k),[],[],[],[],[],[],[],[],[],[],[],...
      kappa,viscCont,walls,wallsInt,wallsExt,updatePreco,vesicle(k));
  
  % running counter for the total number of gmres iterations
  % per time step.  GMRES is used at each Gauss-lobatto time
  % substep which we need to add up for a fair comparison to
  % the traditional time stepping method
  iter = iter + subIter;
  
  % don't want to precompute the preconditioner at the next time
  % step to save CPU time.  Since the configuration doesn't 
  % change by much, the preconditioner shouldn't change
  % by much either and the old preconditioner should
  % work fine
  updatePreco = false;
  
  % if any of the gmres iterations fail, assign the failure to
  % flag to iflag for monitor to output
  if iflagTemp ~= 0
    iflag = iflagTemp;
  end
  
  if o.nsdc > 0
      
    % want to save single-layer potential for computing the residual.
    % This will not save the  one, but this is computed in
    % computeResidual
    Galpert(:,:,:,k) = o.Galpert;
    if any(viscCont ~= 1)
      D(:,:,:,k) = o.D;
    end
    
    % Near structure for vesicles and walls
    NearV2V{k} = o.NearV2V;
    NearW2V{k} = o.NearW2V;
    NearV2W{k} = o.NearV2W;
    NearW2W{k} = o.NearW2W;
    NearWint2V{k} = o.NearWint2V;
    NearWext2V{k} = o.NearWext2V;
    NearV2Wint{k} = o.NearV2Wint;
    NearV2Wext{k} = o.NearV2Wext;
    NearWext2Wint{k} = o.NearWext2Wint;
    NearWint2Wext{k} = o.NearWint2Wext;
    NearWint2Wint{k} = o.NearWint2Wint;

    % need to save if doing SDC corrections 
    vesicle(k+1) = capsules(X,sigma,u,kappa,viscCont,o.antiAlias);
    % compute required upsampling rate
    if o.antiAlias
      vesicle(k+1).setUpRate(o.op);
    end
  end

  % save the provisional solution at the intermediate 
  % Gauss-Lobatto point
  Xprov(:,:,k+1) = X;
  sigmaProv(:,:,k+1) = sigma;
  uProv(:,:,k+1) = u;
  etaProv(:,:,k+1) = eta;
  etaIntProv(:,:,k+1) = etaInt;
  etaExtProv(:,:,k+1) = etaExt;
  RSprov(:,:,k+1) = RS;
 
end % k


% need to save the near-singular integration structure and
% layer-potential matricies at the final Gauss-Lobatto point for future
% SDC sweeps and computing the residual, if there are any SDC sweeps
if o.nsdc > 0
  
  % Need near-singular integration strucutures for final state
  % for sdc updates and computing the residual
  if o.confined
    if ~o.diffDiscWalls % if walls discretized with Nbd
      % Need vesicle to vesicle and vesicle to wall interactions
      [NearV2V{abs(o.orderGL)},NearV2W{abs(o.orderGL)}] = ...
        vesicle(abs(o.orderGL)).getZone(walls,3);
      % Only need wall to vesicle interactions.  Wall to wall 
      % interactions should also use near-singular integration since
      % they may be close to one another
      if nvbd == 1 
        [NearW2W{abs(o.orderGL)},NearW2V{abs(o.orderGL)}] = ...
            walls.getZone(vesicle(abs(o.orderGL)),2);
      else
        [NearW2W{abs(o.orderGL)},NearW2V{abs(o.orderGL)}] = ...
            walls.getZone(vesicle(abs(o.orderGL)),3);
      end
      
      NearV2Wint{abs(o.orderGL)} = [];
      NearWint2V{abs(o.orderGL)} = [];
      NearWint2Wint{abs(o.orderGL)} = [];
      NearV2Wext{abs(o.orderGL)} = [];
      NearWext2V{abs(o.orderGL)} = [];
      NearWint2Wext{abs(o.orderGL)} = [];
      NearWext2Wint{abs(o.orderGL)} = [];
    else % if walls discretized with NbdInt and NbdExt
      [NearV2V{abs(o.orderGL)},NearV2Wint{abs(o.orderGL)}] = ...
        vesicle(abs(o.orderGL)).getZone(wallsInt,3);
      [NearWint2Wint{abs(o.orderGL)},NearWint2V{abs(o.orderGL)}] = ...
        wallsInt.getZone(vesicle(abs(o.orderGL)),3);
      [~,NearV2Wext{abs(o.orderGL)}] = ...
        vesicle(abs(o.orderGL)).getZone(wallsExt,2);
      [~,NearWext2V{abs(o.orderGL)}] = ...
        wallsExt.getZone(vesicle(abs(o.orderGL)),2);
      [~,NearWint2Wext{abs(o.orderGL)}] = ...
        wallsInt.getZone(wallsExt,2);
      [~,NearWext2Wint{abs(o.orderGL)}] = ...
        wallsExt.getZone(wallsInt,2);    
        
      NearV2W{abs(o.orderGL)} = [];
      NearW2V{abs(o.orderGL)} = [];
      NearW2W{abs(o.orderGL)} = [];
    end
  else
    NearV2V{abs(o.orderGL)} = vesicle(abs(o.orderGL)).getZone([],1);
    
    % no solid walls, so only need vesicle-vesicle intearactions
    NearV2W{abs(o.orderGL)} = [];
    NearW2V{abs(o.orderGL)} = [];
    NearW2W{abs(o.orderGL)} = [];
    NearV2Wint{abs(o.orderGL)} = [];
    NearWint2V{abs(o.orderGL)} = [];
    NearWint2Wint{abs(o.orderGL)} = [];
    NearV2Wext{abs(o.orderGL)} = [];
    NearWext2V{abs(o.orderGL)} = [];
    NearWint2Wext{abs(o.orderGL)} = [];
    NearWext2Wint{abs(o.orderGL)} = [];
  end
  
  % need single-layer potential at the final state of the 
  % provisional solution to form the residual
  Galpert(:,:,:,abs(o.orderGL)) = o.op.stokesSLmatrix(vesicle(o.orderGL));
  if any(viscCont ~= 1)
    D(:,:,:,abs(o.orderGL)) = o.op.stokesDLmatrix(vesicle(o.orderGL));
  end
  
end
% END OF FORMING PROVISIONAL SOLUTION AND THE SINGLE-LAYER
% POTENTIAL AT THESE INTERMEDIATE STATES

o.dt = dt;
% change back to original time step
 
%color = ['r' 'g' 'b' 'k' 'c' 'm' 'y'];
%figure(2); clf; hold on
%for k = 1:abs(o.orderGL)
%  plot(sigmaProv(:,1,k),color(k))
%end
%figure(3); clf; hold on
%for k = 1:abs(o.orderGL)
%  plot(squeeze(etaProv(1:end/2,1,k)),color(k))
%end
%figure(5);clf;hold on;
%for k = 1:abs(o.orderGL)
%  plot(etaProv(end/2+1:end,1,k),color(k))
%end
%pause
%% DEBUG: TO MAKE SURE THAT THE TENSION IS CONTINUOUS FROM BETWEEN
%% THE DIFFERENT GAUSS-LOBATTO POINTS
% save the very first provisional solution so that it can be taken
% as the solution if sdc corrections increase the error but this
% original provisional solution has the smallest error
Xinit = Xprov;
uInit = uProv;
sigmaInit = sigmaProv;
etaInit = etaProv;
etaIntInit = etaIntProv;
etaExtInit = etaExtProv;
RSinit = RSprov;

% compute the integrand that is in the residual and the residual 
% of the picard integral formulation of the dynamic equation
% If the residual is not desired, then errors is set to 0 and
% this will never refine the time step size
if o.nsdc > 0
  % form the residual as well as the velocities due to the different
  % components which are necessary to form the right-hand sides when
  % doing sdc corrections  
  [vesVel,divVesVel,wallVel,wallIntVel,wallExtVel,residual] = ...
      o.computeVelocities(vesicle,Galpert,Gbarnett,D,walls,wallsInt,...
      wallsExt,etaProv,etaIntProv,etaExtProv,RSprov,...
      NearV2V,NearV2W,NearW2V,NearW2W,NearV2Wint,NearWint2V,NearV2Wext,...
      NearWext2V,NearWint2Wext,NearWint2Wint,NearWext2Wint);

  % save the size of the residual at the final Gauss-Lobatto
  % point.  We may have more of these later if we do a full
  % deferred correction method.  Use the maximum L2 error where
  % the maximum is taken over all the vesicles
  normRes = max(abs(residual(:,:,end)));
  normRes = max(normRes);
else
  normRes = 0;
end

% print the provisional solution's residual and errors
[~,a,l] = oc.geomProp(Xprov(:,:,end));
% THIS PART HAS BEEN REMOVED FOR STREAMING, I ADDED THIS AND CAN IGNORE
% eaBefore = abs(a - a0)./abs(a0);
% elBefore = abs(l - l0)./abs(l0);
% message = ['sdcCount   ' num2str(0,'%2d') ...
%   ': Residual = ' num2str(normRes(end),'%5.2e') ...
%   ', eA = ' num2str(max(eaBefore),'%5.4e') ...
%   ', eL = ' num2str(max(elBefore),'%5.4e')];
% om.writeMessage(message,'%s\n');


% Start doing SDC corrections
for sdcCount = 1:o.nsdc
  % use preco that was formed at the provisional solution
  updatePreco = false; 
  for n = 1:abs(o.orderGL) - 1
    o.dt = dtimes(n);

    o.SDCcorrect = true;
    o.Galpert = Galpert(:,:,:,n+1);
    if any(viscCont ~= 1)
      o.D = D(:,:,:,n+1);
    end
    
    o.NearV2V = NearV2V{n+1};
    o.NearV2W = NearV2W{n+1};
    o.NearW2V = NearW2V{n+1};
    o.NearW2W = NearW2W{n+1};
    o.NearV2Wint = NearV2Wint{n+1};
    o.NearWint2V = NearWint2V{n+1};
    o.NearV2Wext = NearV2Wext{n+1};
    o.NearWext2V = NearWext2V{n+1};
    o.NearWint2Wext = NearWint2Wext{n+1};
    o.NearWint2Wint = NearWint2Wint{n+1};
    o.NearWext2Wint = NearWext2Wint{n+1};
 
  
    % Form the sdc update
    [X,sigma,u,eta,etaInt,etaExt,RS,subIter,iflagTemp] = o.timeStep(...
        Xprov(:,:,n+1),sigmaProv(:,:,n+1),...
        uStore,etaProv(:,:,n+1),etaIntProv(:,:,n+1),etaExtProv(:,:,n+1),...
        RSprov(:,:,n+1),...
        deltaX(:,:,n),deltaSigma(:,:,n),...
        deltaEta(:,:,n),deltaEtaInt(:,:,n),deltaEtaExt(:,:,n),...
        deltaRS(:,:,n),...
        residual(:,:,n+1) - residual(:,:,n),...
        vesVel(:,:,n+1),wallVel(:,:,n+1),wallIntVel(:,:,n+1),...
        wallExtVel(:,:,n+1),kappa,viscCont,walls,wallsInt,wallsExt,...
        updatePreco,vesicle(n+1),vesicle(1).sa,vesicle(1).IK);
 
    iter = iter + subIter;
    if iflagTemp ~= 0
      iflag = iflagTemp;
    end
    % turn correct off since it is not used when forming the
    % provisional solution
    o.SDCcorrect = false;
    
    % approximations of the error
    deltaX(:,:,n+1) = X;
    deltaSigma(:,:,n+1) = sigma;
    deltaEta(:,:,n+1) = eta;
    deltaEtaInt(:,:,n+1) = etaInt;
    deltaEtaExt(:,:,n+1) = etaExt;
    deltaRS(:,:,n+1) = RS;
    
  end
  % go back to original time step
  o.dt = dt;
  
  % update provision solution
  alpha = 1;
  Xprov = Xprov + alpha*deltaX;
  
  [~,a,l] = oc.geomProp(Xprov(:,:,end));
  ea = abs(a - a0)./abs(a0);
  el = abs(l - l0)./abs(l0);

%   if max(ea) > max(eaBefore) || max(el) > max(elBefore)
%     message = ['SDC increases the error, do not accept corrected solution'];
%     om.writeMessage(message,'%s\n');
%     % Do not accept SDC correction
%     Xprov = Xprov - alpha*deltaX;
%   else
    sigmaProv = sigmaProv + alpha*deltaSigma;
    etaProv = etaProv + alpha*deltaEta;
    etaIntProv = etaIntProv + alpha*deltaEtaInt;
    etaExtProv = etaExtProv + alpha*deltaEtaExt;
    RSprov = RSprov + alpha*deltaRS;  
    
%   end
  
  normRes = max(abs(residual(:,:,end)));
  normRes = max(normRes);

  message = ['sdcCount   ' num2str(sdcCount,'%2d') ...
    ': Residual = ' num2str(normRes(end),'%5.2e') ...
    ', eA = ' num2str(max(ea),'%5.4e') ...
    ', eL = ' num2str(max(el),'%5.4e')];
  om.writeMessage(message,'%s\n');
  
  if sdcCount < o.nsdc
    for n = 1:abs(o.orderGL)
      % update the vesicle objects and the single-layer potentials
      vesicle(n) = capsules(Xprov(:,:,n),sigmaProv(:,:,n),...
          [],kappa,viscCont,o.antiAlias);
      % find correct upsampling rate
      if o.antiAlias
        vesicle(n).setUpRate(o.op);
      end
      
      Galpert(:,:,:,n) = o.op.stokesSLmatrix(vesicle(n));
      if any(viscCont ~= 1)
        D(:,:,:,n) = o.op.stokesDLmatrix(vesicle(n));
      else
        D = [];
      end
    end
    
    [vesVel,divVesVel,wallVel,wallIntVel,wallExtVel,residual] = ...
      o.computeVelocities(vesicle,Galpert,Gbarnett,D,walls,wallsInt,...
      wallsExt,etaProv,etaIntProv,etaExtProv,RSprov,NearV2V,NearV2W,...
      NearW2V,NearW2W,NearV2Wint,NearWint2V,NearV2Wext,NearWext2V,...
      NearWint2Wext,NearWint2Wint,NearWext2Wint);
  end
  % if we are only recording the error in area and length, don't need
  % to compute the final residual
  
end
% End of doing SDC corrections
% update solution with an SDC iteration

if o.timeAdap 
  % if doing adaptive time stepping, get new time step
  % size and time step scaling
  [accept,dtScale,collUprate] = ...
      o.newTimeStepSizeDerivs(Xprov(:,:,end),a,l,...
      Xinit_tadap,aInit,lInit,walls,wallsInt,wallsExt);
else
  % if not doing adaptive time stepping, keep the same
  % time step size and always accept solution  
  accept = true;
  dtScale = 1;
  collUprate = 1;
end

if accept
  % take the new solution
  X = Xprov(:,:,end);
  sigma = sigmaProv(:,:,end);
  u = uProv(:,:,end);
  eta = etaProv(:,:,end);
  etaInt = etaIntProv(:,:,end);
  etaExt = etaExtProv(:,:,end);
  RS = RSprov(:,:,end);
else
  % revert to the old solution
  X = Xstore;
  sigma = sigStore;
  u = uStore;
  eta = etaStore;
  etaInt = etaIntProv(:,:,end);
  etaExt = etaExtProv(:,:,end);
  RS = RSstore;
end

if accept && o.periodic
  X = oc.addAndRemove(X,walls,o.near,o.fmm);
end

end % timeStepGL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [accept,dtScale] = newTimeStepSize(o,...
      aFinal,lFinal,aInit,lInit,accept,om)
  
% THIS IS BRYAN'S ADAPTIVE TIME STEPPING WHICH IS NOT USED IN THE CURRENT
% VERSION OF THE CODE. THIS IS NOT UP-TO-DATE.
  
% [accept,dtScale] = newTimeStepSize(...
%      aFinal,lFinal,aInit,lInit,accept,om)
% finds a new time step size based on the change in the area and length
% from the first to final Gauss- Lobatto points.  The output is a flag
% indicating acceptance or rejection, and the amount the time step is
% scaled.  To increase the likelihood that the a new time step size is
% accepted, the time step size is never increased if the previous time
% step was rejected

alpha = o.alpha;
% buffer so that we aren't trying to keep error exactly at 1.  This is a
% safeguard so that the next time step is accepted with higher
% probability
betaUp = o.betaUp;
% allowable upscaling change in time step size
betaDown = o.betaDown;
% allowable downscaling change in time step size

errArea = abs(aInit - aFinal);
errLength = abs(lInit - lFinal);
% absolute errors in area and length
tauArea = max(aInit,aFinal)*o.areaLenTol*o.dt;
tauLength = max(lInit,lFinal)*o.areaLenTol*o.dt;
% Tolerance for errArea and errLength
err = max(max(errArea./tauArea),max(errLength./tauLength));
% Maximum of relative errors in area and length
% Want this quantity to be as close to 1 as possible

actualOrder = o.expectedOrder;
% order of time stepping method
dtOPT = err^(-1/(actualOrder))*o.dt;
% optimal time step size

dtOld = o.dt;
if accept
  o.dt = alpha^(1/actualOrder) * ...
      min(betaUp*o.dt,max(dtOPT,betaDown*o.dt));
else
  o.dt = alpha^(1/(actualOrder)) * ...
      min(o.dt,max(dtOPT,betaDown*o.dt));
  % don't want to scale up this time step if it was previously rejected.
  % This hopefully gets rid of pattern where a solution alternates
  % between being accepted and rejected.
end
% safety factor added to the optimal time step size also, time step size
% is not scaled up or down too fast For safety factor, take 1/p root.
% In our formulation, this makes alpha the desired value for err
% regardless of the order of the method.
dtScale = o.dt/dtOld;
% time step scaling

if err > 1
  accept = false;
  % reject time step because the error is too large
  message = ['Time Step REJECTED with error ' ...
      num2str(err,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['Time step scaled by           ' ... 
      num2str(dtScale,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['New time step size is         ' ...
      num2str(o.dt,'%4.2e')];
  om.writeMessage(message,'%s\n')
  om.writeMessage(' ','%s\n')
else
  accept = true;
  % accept the solution because the error is small
  message = ['Time Step ACCEPTED with error ' ...
      num2str(err,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['Time step scaled by           ' ... 
      num2str(dtScale,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['New time step size is         ' ...
      num2str(o.dt,'%4.2e')];
  om.writeMessage(message,'%s\n')
end

end % newTimeStepSize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [accept,dtScale,collUprate] = ...
  newTimeStepSizeDerivs(o,X,aTdT,lTdT,Xold,aT,lT,walls,wallsInt,wallsExt)
% [accept,dtScale,collUprate] = ...
%   newTimeStepSizeDerivs(o,X,aTdT,lTdT,Xold,aT,lT,walls,wallsInt,wallsExt)

% finds a new time step size based on the change in the area and length
% after taking a time step and the derivatives dL/dt
%, dA/dt. This is adjusted for low-resolution simulations. If it decides 
% that the asymptotic assumption (Bryan's formulation above) can be used, 
% it uses that to find the new time step. This requires high 
% spatio-temporal resolution. If this is not the case, then it simply 
% increases or decreases the solution by a constant factor. 
% It also checks if there is any collision, if so the solution is rejected.
% The output is a flag
% indicating acceptance or rejection, and the amount the time step is
% scaled. 

om = o.om;
oc = curve;

nv = size(X,2);
% target error which we are allowed to commit at each time step
errTarget = o.areaLenTol; % given tolerance (\rho_AL)
errTargetMin = o.areaLenTol * 0.5; % creates a buffer zone (\rho_min)

% max and min allowable time step sizes
dtMax = o.dtMax;
dtMin = o.dtMin;

% buffer in increasing the time step size (fixed to 0.9)
alpha = o.alpha;

% constant factor by which the time step size is increased
betaUp = o.betaUp; 

% constant factor by which the time step size is reduced
betaDown = o.betaDown; 

% area and length of the provisional solution
aTdT = aTdT'; lTdT = lTdT';

% area and length at the previous time step
aT = aT'; lT = lT';

% compute the derivatives
[qA,qL] = oc.computeAreaLengthDerivs(X,Xold,o.dt);

qA = abs(qA)./aT; % dA/dt / A(t)
qL = abs(qL)./lT; % dL/dt / L(t)

% tolerances which allow us to use derivatives to compute new dt
betaInc = o.betaInc; % 1e-3, we are more careful in increasing dt
betaDec = o.betaDec; % 1e-1, we are more aggresive in decreasing dt

% COMPUTE dt BASED ON AREA

% (actual) error committed in this time step
errActArea = abs(aTdT-aT)./aT;

if max(errActArea) < errTarget && max(errActArea) >= errTargetMin
    acceptA = true;
    dtNewA = min(o.dt,dtMax);
    
elseif max(errActArea) < errTargetMin 
    acceptA = true;
    if max(abs(qA*o.dt-errActArea)./errActArea) < betaInc
        dtNewA = alpha*errTarget/max(abs(qA));
    else
        dtNewA = betaUp * o.dt;
    end
    % if the asymptotic assumption reduces the time step size, increase by
    % constant factor
    if dtNewA < o.dt
        dtNewA = betaUp * o.dt;
    end
    dtNewA = min(dtNewA,dtMax);
elseif max(errActArea) >= errTarget
    acceptA = false;
    if max(abs(qA*o.dt-errActArea)./errActArea) < betaDec
        dtNewA = alpha*errTarget/max(abs(qA));
    else
        dtNewA = betaDown * o.dt;
    end
    % if the asymptotic assumption increases or does not change 
    % the time step size, decrease by constant factor
    if 1.1*dtNewA >= o.dt
        dtNewA = betaDown * o.dt;
    end
    
end

% COMPUTE dt BASED ON LENGTH

% (actual) error committed in this time step
errActLength = abs(lTdT-lT)./lT;

if max(errActLength) < errTarget && max(errActLength) >= errTargetMin
    acceptL = true;
    dtNewL = min(o.dt,dtMax);
    
elseif max(errActLength) < errTargetMin
    acceptL = true;
    if max(abs(qL*o.dt-errActLength)./errActLength) < betaInc
        dtNewL = alpha*errTarget/max(abs(qL));
    else
        dtNewL = betaUp * o.dt;
    end
    % if the asymptotic assumption reduces the time step size, increase by
    % constant factor
    if dtNewL < o.dt
        dtNewL = betaUp * o.dt;
    end
    dtNewL = min(dtNewL,dtMax);
elseif max(errActLength) >= errTarget
    acceptL = false;
    if max(abs(qL*o.dt-errActLength)./errActLength) < betaDec
        dtNewL = alpha*errTarget/max(abs(qL));
    else
        dtNewL = betaDown * o.dt;
    end
    % if the asymptotic assumption increases or does not change 
    % the time step size, decrease by constant factor
    if 1.1*dtNewL >= o.dt
        dtNewL = betaDown * o.dt;
    end
    
end
% Time step size based on error in length and error
dtNewAL = min(dtNewA,dtNewL);
err = max(max(errActArea),max(errActLength));
acceptAL = (acceptA & acceptL);

% Collision detection
% See if we need to upsample, if so upsample
vesicleProv = capsules(X,[],[],[],[],o.antiAlias);
if o.antiAlias
  vesicleProv.setUpRate(o.op);
  collUprate = vesicleProv.uprate;
  Nup = size(X,1)/2 * collUprate;
  Xup = [interpft(X(1:end/2,:),Nup);interpft(X(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],[],[],o.antiAlias);
else
  vesicleUp = vesicleProv;
  collUprate = 1;
end
  
% Need near structure
if o.confined
  if ~o.diffDiscWalls  
    [NearV2V,NearV2W] = vesicleUp.getZone(walls,3);
  else
    [NearV2V,NearV2Wint] = vesicleUp.getZone(wallsInt,3);  
    [~,NearV2Wext] = vesicleUp.getZone(wallsExt,2);
  end
else
  NearV2V = vesicleUp.getZone([],1);  
  NearV2W = [];
  NearV2Wint = [];
  NearV2Wext = [];
end
  
if ~o.diffDiscWalls
  [icollisionVes,icollisionWall] = ...
    vesicleUp.collision(walls,NearV2V,NearV2W,o.fmm,o.op);
else
  % in some DLD cases wallsExt has too many points, this would slow down 
  % the code to check collision. Usually vesicle stays away from the 
  % exterior wall in DLD examples.
  if wallsExt.N <= 2048  
    [~,icollisionWallExt] = vesicleUp.collision(wallsExt,...
      NearV2V,NearV2Wext,o.fmm,o.op);
  else
    icollisionWallExt = 0;
  end
  [icollisionVes,icollisionWallInt] = vesicleUp.collision(wallsInt,...
    NearV2V,NearV2Wint,o.fmm,o.op);
  icollisionWall = icollisionWallInt || icollisionWallExt;
end
% Check for collisions 
acceptColl = true;
dtNewColl = dtNewAL;

if icollisionVes

  message = ['VESICLES HAVE CROSSED'];
  om.writeMessage(message,'%s\n')
  
  global ncollVes
  ncollVes = ncollVes + 1;
  dtNewColl = o.dt * betaDown;
  acceptColl = false;
end
if icollisionWall
  
  message = ['VESICLES HAVE CROSSED SOLID WALL'];
  om.writeMessage(message,'%s\n')
  
  global ncollWal
  ncollWal = ncollWal + 1;
  dtNewColl = o.dt * betaDown;
  acceptColl = false;
end
% check for collisions


% choose the minimum time step size
dtOld = o.dt;
o.dt = min(dtNewColl,dtNewAL);
accept = acceptAL & acceptColl;
% if vesicles cross, reject time step size anyways and choose the smallest
% dt

% if dt is required to be smaller then dtMin then do not accept dtMin and
% go down
if ~acceptColl 
  if o.dt <= dtMin
    message = '!! Time step size is less than the minimum value b/c of collision';
    om.writeMessage(message,'%s\n')
  end
else
  if o.dt <= dtMin
    o.dt = dtMin;
    accept = true;
    message = 'Error in area-length requires smaller time step size than the minimum';
    om.writeMessage(message,'%s\n')
    message = 'So, go on with the minimum time step size, comprimise the tolerance';
    om.writeMessage(message,'%s\n')
  end
end

dtScale = o.dt/dtOld;
% time step scaling

if ~accept
  % reject time step because the error is too large or collision occurs
  message = ['Time Step REJECTED with error ' ...
      num2str(err,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['Time step scaled by           ' ... 
      num2str(dtScale,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['New time step size is         ' ...
      num2str(o.dt,'%4.2e')];
  om.writeMessage(message,'%s\n')
  om.writeMessage(' ','%s\n')
else
  % accept the solution because the error is small
  message = ['Time Step ACCEPTED with error ' ...
      num2str(err,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['Time step scaled by           ' ... 
      num2str(dtScale,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['New time step size is         ' ...
      num2str(o.dt,'%4.2e')];
  om.writeMessage(message,'%s\n')
end

end % newTimeStepSizeDerivs

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,sigma,u,eta,etaInt,etaExt,RS,iter,iflag] = timeStep(o,...
    Xstore,sigStore,uStore,etaStore,etaIntStore,etaExtStore,RSstore,...
    deltaX,deltaSig,deltaEta,deltaEtaInt,deltaEtaExt,deltaRS,...
    diffResidual,vesVel,wallVel,wallIntVel,wallExtVel,kappa,viscCont,...
    walls,wallsInt,wallsExt,updatePreco,vesicle,sa,IK)

% [X,sigma,u,eta,etaInt,etaExt,RS,iter,iflag] = timeStep(o,...
%     Xstore,sigStore,uStore,etaStore,etaIntStore,etaExtStore,RSstore,...
%     deltaX,deltaSig,deltaEta,deltaEtaInt,deltaEtaExt,deltaRS,...
%     diffResidual,vesVel,wallVel,wallIntVel,wallExtVel,kappa,viscCont,...
%     walls,wallsInt,wallsExt,updatePreco,vesicle,sa,IK)
% uses implicit vesicle-vesicle interactions and
% discretizes the inextensibility condition in three different ways
% (method1, method2, or method 3).  Must pass in the vesicle positions,
% tension, and velocity from enough previous time steps (depends on
% o.order).  Returns a new positions, tension, velocity, density
% function defined on the solid walls and the number of required GMRES
% iterations if o.SDCcorrect=true, then it uses deltaX, deltaSig, etc
% to compute the right-hand sides needed for sdc updates updatePreco is
% a flog that decides if the block-diagonal preconditioner should be
% updated or not
% NOTE THAT A LOT OF THE FEATURES ARE NOT IMPLEMENTED WITH EXPLICIT
% VESICLE-VESICLE INTERACTIONS.

N = size(Xstore,1)/2; % Number of points per vesicle
nv = size(Xstore,2); % Number of vesicles
if o.confined
  if ~o.diffDiscWalls
    Xwalls = walls.X; % discretization points of solid walls
    XwallsInt = [];
    XwallsExt = [];
  else
    XwallsInt = wallsInt.X; % interior walls
    XwallsExt = wallsExt.X;
    Xwalls = [];
  end
else
  Xwalls = [];
  XwallsInt = [];
  XwallsExt = [];
end
Nbd = size(Xwalls,1)/2; % Number of points on the solid walls
NbdInt = size(XwallsInt,1)/2; % Number of points on the interior walls
NbdExt = size(XwallsExt,1)/2; % Number of points on the exterior walls
% number of solid wall components
nvbdInt = size(XwallsInt,2); % # of interior walls
nvbdExt = size(XwallsExt,2); % # of exterior walls
nvbdSme = size(Xwalls,2);    % # of walls of the same discretization
nvbd = nvbdSme + nvbdInt + nvbdExt;
% constant that appears in front of time derivative in
% vesicle dynamical equations
alpha = (1 + viscCont)/2; 

% Form linear combinations of previous time steps needed for Ascher,
% Ruuth, and Wetton IMEX methods
Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
uM = zeros(2*N,nv);
Xo = zeros(2*N,nv);
etaM = zeros(2*Nbd,nvbdSme);
etaMint = zeros(2*NbdInt,nvbdInt);
etaMext = zeros(2*NbdExt,nvbdExt);
RSm = zeros(3,nvbd);
for k = 1:o.order
  Xm = Xm + Xstore(:,:,k)*o.Xcoeff(k);
  sigmaM = sigmaM + sigStore(:,:,k)*o.Xcoeff(k);
  uM = uM + uStore(:,:,k)*o.Xcoeff(k);
  if o.confined
    if ~o.diffDiscWalls
      etaM = etaM + etaStore(:,:,k)*o.Xcoeff(k);
    else
      etaMint = etaMint + etaIntStore(:,:,k)*o.Xcoeff(k);
      etaMext = etaMext + etaExtStore(:,:,k)*o.Xcoeff(k);
    end
    RSm = RSm + RSstore(:,:,k)*o.Xcoeff(k);
  end
  Xo = Xo + Xstore(:,:,k)*o.rhsCoeff(k);
end

% build an object vesicle that contains tangent vector, jacobian, etc.
if o.order ~= 1
  vesicle = capsules(Xm,sigmaM,uM,kappa,viscCont,o.antiAlias);
  if o.antiAlias
    vesicle.setUpRate(o.op);
  end
end

% Build single layer potential matrix and put it in current object
% If we are doing an sdc update, this is already precomputed and 
% stored from when we formed the provisional solution
op = o.op;
if ~o.SDCcorrect
  o.Galpert = op.stokesSLmatrix(vesicle);
%  o.Gbarnett = op.laplaceSLcomplexMatrix(vesicle);
end


% Compute double-layer potential matrix due to each vesicle
% independent of the others.  Matrix is zero if there is no
% viscosity contrast
if ~o.SDCcorrect
  if any(viscCont ~= 1)
    o.D = op.stokesDLmatrix(vesicle);
    if o.fmmDLP
      Nup = op.LPuprate * vesicle.N;
      Xup = [interpft(vesicle.X(1:end/2,:),Nup);...
        interpft(vesicle.X(end/2+1:end,:),Nup)];
      vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0);
      o.DLPnoCorr = op.stokesDLmatrixNoCorr(vesicleUp);
     
    else
      o.DLPnoCorr = [];
    end
  else
    o.D = [];
    o.DLPnoCorr = [];
  end
  
  if o.fmm
    Nup = op.LPuprate * vesicle.N;
    Xup = [interpft(vesicle.X(1:end/2,:),Nup);...
      interpft(vesicle.X(end/2+1:end,:),Nup)];
    vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0); 
    o.SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);
  else
    o.SLPnoCorr = [];    
  end
end

% Only form near-singular integration structure if not doing an SDC
% update.  Otherwise this was formed and saved when forming the
% provisional solution
if ~o.SDCcorrect
  % Structures for deciding who
  % is close, how close it is, who is closest, etc., needed in nearSingInt
  if o.confined
    if~o.diffDiscWalls
      % Need vesicle to vesicle and vesicle to wall interactions
      [o.NearV2V,o.NearV2W] = vesicle.getZone(walls,3);
      
      % Only need wall to vesicle interactions.  Wall to wall
      % interactions should also use near-singular integration since
      % they may be close to one another
      if nvbd == 1
        [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,2);
      else
        if isempty(o.NearW2W)
          [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,3);
        else
          % there is no need to compute W2W again, since they do not move  
          [~,o.NearW2V] = walls.getZone(vesicle,2);
        end
      end
      o.NearV2Wint = [];
      o.NearWint2Wint = [];
      o.NearWint2V = [];
      o.NearV2Wext = [];
      o.NearWext2V = [];
      o.NearWint2Wext = [];
      o.NearWext2Wint = [];
    else
      % When walls are discretized separately, we cannot call getZone 
      % for all the walls. That function is not adjusted for that. Instead
      % this if statement is added and near structures are formed
      % separately.
      [o.NearV2V,o.NearV2Wint] = vesicle.getZone(wallsInt,3);
      
      if isempty(o.NearWint2Wint)
        [o.NearWint2Wint,o.NearWint2V] = wallsInt.getZone(vesicle,3);
      else
        % there is no need to compute W2W again, since they do not move
        [~,o.NearWint2V] = wallsInt.getZone(vesicle,2);
      end
      
      [~,o.NearV2Wext] = vesicle.getZone(wallsExt,2);
      [~,o.NearWext2V] = wallsExt.getZone(vesicle,2);
      
      if isempty(o.NearWint2Wext)
        % there is no need to compute W2W again, since they do not move
        [~,o.NearWint2Wext] = wallsInt.getZone(wallsExt,2);
      end
      
      if isempty(o.NearWext2Wint)
        % there is no need to compute W2W again, since they do not move
        [~,o.NearWext2Wint] = wallsExt.getZone(wallsInt,2);    
      end

      o.NearV2W = [];
      o.NearW2W = [];
      o.NearW2V = [];
    end
  else
    % no solid walls, so only need vesicle-vesicle intearactions
    o.NearV2V = vesicle.getZone([],1);
    
    o.NearV2W = [];
    o.NearW2V = [];
    o.NearW2W = [];
    o.NearV2Wint = [];
    o.NearWint2Wint = [];
    o.NearWint2V = [];
    o.NearV2Wext = [];
    o.NearWext2V = [];
    o.NearWint2Wext = [];
    o.NearWext2Wint = [];
  end
  
end

% Parts of rhs from previous solution.  The right-hand-side depends on
% whether we are doing an SDC correction or forming the provisional
% solution.
if ~o.SDCcorrect
  rhs1 = Xo;
  rhs2 = zeros(N,nv);
  if o.confined
    if ~o.diffDiscWalls  
      rhs3 = walls.u;
      rhs3Ext = [];
      rhs3Int = [];
    else
      rhs3Ext = wallsExt.u;
      rhs3Int = wallsInt.u;
      rhs3 = [];
    end
  else
    rhs3 = [];
    rhs3Ext = [];
    rhs3Int = [];
  end
else
  rhs1 = deltaX + diffResidual;
  if any(vesicle.viscCont ~= 1)
    z = zeros(2*N*nv,1);
    for k = 1:nv
      z(2*(k-1)*N+1:2*k*N) = rhs1(:,k);
    end
    z = o.IminusD(z,vesicle);
    for k = 1:nv
      rhs1(:,k) = z(2*(k-1)*N+1:2*k*N)/alpha(k);
    end
  end
  if strcmp(o.solver,'method1')
    rhs2 = ones(N,nv);
  else
    rhs2 = -vesicle.surfaceDiv(vesVel);
  end
  if o.confined
    if~o.diffDiscWalls
      rhs3 = -wallVel + walls.u;
      rhs3Ext = [];
      rhs3Int = [];
    else
      rhs3Ext = -wallExtVel + wallsExt.u;
      rhs3Int = -wallIntVel + wallsInt.u;
      rhs3 = [];
    end
  else
    rhs3 = [];
    rhs3Ext = [];
    rhs3Int = [];
  end
end


% START TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP
% vesicle-vesicle and vesicle-wall interactions are handled
% implicitly in TimeMatVec
% END TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VISCOSITY CONTRAST
if any(vesicle.viscCont ~= 1)
  
  % Need to add jump to double-layer potential if using near-singular
  % integration so that we can compute boundary values for the
  % near-singular integration algorithm
  jump = 1/2*(1-vesicle.viscCont);
  DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);

  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
    kernelDirect = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLnewfmm;
    kernelDirect = @op.exactStokesDL;
  end
  
  if ~o.SDCcorrect
    density = Xo;
  else
    density = deltaX;
  end

  % Use near-singular integration to compute double-layer
  % potential from previous solution 
  Fdlp = op.nearSingInt(vesicle,density,DLP,o.DLPnoCorr,...
       o.NearV2V,kernel,kernelDirect,vesicle,true,false);
      
  if o.confined
    if ~o.diffDiscWalls  
      FDLPwall = op.nearSingInt(vesicle,density,DLP,[],...
        o.NearV2W,kernel,kernelDirect,walls,false,false);
      FDLPwallInt = [];
      FDLPwallExt = [];
    else
      FDLPwallInt = op.nearSingInt(vesicle,density,DLP,[],o.NearV2Wint,...
          kernel,kernelDirect,wallsInt,false,false);
      FDLPwallExt = op.nearSingInt(vesicle,density,DLP,[],o.NearV2Wext,...
          kernel,kernelDirect,wallsExt,false,false);
      FDLPwall = [];
    end
  else
    FDLPwall = [];
    FDLPwallInt = [];
    FDLPwallExt = [];
  end
    
else
  % If no viscosity contrast, there is no velocity induced due to a
  % viscosity contrast  
  Fdlp = zeros(2*N,nv);
  FDLPwall = zeros(2*Nbd,nvbdSme);
  FDLPwallInt = zeros(2*NbdInt,nvbdInt);
  FDLPwallExt = zeros(2*NbdExt,nvbdExt);
end

% add in viscosity contrast term due to each vesicle independent of the
% others (o.D * Xo) from the previous solution followed by the term due
% to all other vesicles (Fdlp)
if (any(viscCont ~= 1) && ~o.SDCcorrect)
  DXo = op.exactStokesDLdiag(vesicle,o.D,Xo);
  rhs1 = rhs1 - (Fdlp + DXo) * diag(1./alpha);
end

% compute the double-layer potential due to all other vesicles from the
% appropriate linear combination of previous time steps.  Depends on
% time stepping order and vesicle-vesicle discretization
rhs3 = rhs3 + FDLPwall/o.dt;
rhs3Ext = rhs3Ext + FDLPwallExt/o.dt;
rhs3Int = rhs3Int + FDLPwallInt/o.dt;


% START COMPUTING SINGLE-LAYER POTENTIAL FOR REPULSION
if o.repulsion 
  % Repulsion is handled explicitly between vesicles, vesicles-walls.  
  if ~o.fmm
      kernel = @op.exactStokesSL;
      kernelDirect = @op.exactStokesSL;
  else
      kernel = @op.exactStokesSLfmm;
      kernelDirect = @op.exactStokesSL;
  end
  
  if ~o.SDCcorrect
    Xrep = Xo;
  else
    Xrep = deltaX;
  end
  
  if ~o.confined
    repulsion = vesicle.repulsionScheme(Xrep,o.repStrength,o.minDist,...
        [],[],[]);
  else
    if ~o.diffDiscWalls  
      repulsion = vesicle.repulsionScheme(Xrep,o.repStrength,o.minDist,...
          walls,[],[]);
    else
      repulsion = vesicle.repulsionScheme(Xrep,o.repStrength,o.minDist,...
          [],wallsInt,wallsExt);
    end
  end

  Frepulsion = op.exactStokesSLdiag(vesicle,o.Galpert,repulsion);
  % diagonal term of repulsion


  SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
    
  % Use near-singular integration to compute single-layer potential
  % due to all other vesicles.  Need to pass function
  % op.exactStokesSL so that the layer potential can be computed at
  % far points and Lagrange interpolation points
  Frepulsion = Frepulsion + ...
    op.nearSingInt(vesicle,repulsion,SLP,o.SLPnoCorr,...
        o.NearV2V,kernel,kernelDirect,vesicle,true,false);
    
  % Evaluate the velocity on the walls due to the vesicles
  if o.confined
      if ~o.diffDiscWalls
        FREPwall = op.nearSingInt(vesicle,repulsion,SLP,[],...
            o.NearV2W,kernel,kernelDirect,walls,false,false);
        FREPwallInt = [];
        FREPwallExt = [];
      else
        FREPwallInt = op.nearSingInt(vesicle,repulsion,SLP,[],...
            o.NearV2Wint,kernel,kernelDirect,wallsInt,false,false);    
        FREPwallExt = op.nearSingInt(vesicle,repulsion,SLP,[],...
            o.NearV2Wext,kernel,kernelDirect,wallsExt,false,false);    
        FREPwall = [];
      end
  else
      FREPwall = [];
      FREPwallInt = [];
      FREPwallExt = [];
  end
    
  % keep the number of timesteps when repulsion is nonzero
  if norm(Frepulsion) ~= 0 
    global repuls
    repuls = repuls + 1;
  end
    
  rhs1 = rhs1 + o.dt*Frepulsion*diag(1./alpha);
  rhs3 = rhs3 - FREPwall;
  rhs3Ext = rhs3Ext - FREPwallExt;
  rhs3Int = rhs3Int - FREPwallInt;
end
% END COMPUTING SINGLE-LAYER POTENTIALS FOR REPULSION 



% START TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS
% This is done implicitly in TimeMatVec
if ~o.confined
  if ~o.SDCcorrect
    % Add in far-field condition (extensional, shear, etc.)
    vInf = o.farField(Xm,[]);
    rhs1 = rhs1 + o.dt*vInf*diag(1./alpha);
    if o.saveVinf
      file = ['./output/' o.runName '_vInf.bin'];
      fid = fopen(file,'a');
      fwrite(fid,vInf(:),'double');
      fclose(fid);
    end
  end
end
% END TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS


% START TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION
if (strcmp(o.solver,'method1'))
  % If using method1, vesicle-vesicle interactions and the presence of a
  % viscosity contrast is irrelavent
  
  % rhs2 is the right-hand side for the inextensibility condition
  if ~o.SDCcorrect
    rhs2 = rhs2 + vesicle.surfaceDiv(Xo); 
  else
    divf = zeros(N,nv);
    for k = 1:nv
      divf(:,k) = curve.arcDeriv(Xm(1:N,k),1,1./sa(:,k),...
          IK(:,k)).^2 + ...
                  curve.arcDeriv(Xm(N+1:2*N,k),1,1./sa(:,k),...
          IK(:,k)).^2;
    end
    rhs2 = 1/2*(rhs2 - divf);
  end
else 
  % If using method2, method3, or method4, vesicle-vesicle interaction
  % affects the right-hand side of the inextensibility condition

  if ~o.confined && ~o.SDCcorrect
    if any(viscCont ~= 1)
      if o.profile
        tic
      end
      rhs2 = rhs2 - vesicle.surfaceDiv(...
          o.solveIminusD(o.farField(Xm,[]),vesicle));
      if o.profile
        fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
      end
    else
      rhs2 = rhs2 - vesicle.surfaceDiv(o.farField(Xm,[]));
    end
  end
  % add in term from farfield
end
% END TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION

% The next makes sure that the rhs is all order one rather than have rhs3
% being order 1/o.dt and other two parts (rhs1 and rhs2) being order 1.
% This of course needs to be compensated in the TimeMatVec routine
if (any(vesicle.viscCont ~= 1) && ...
      strcmp(o.vesves,'implicit') && o.confined)
  rhs3 = rhs3 * o.dt;
  rhs3Int = rhs3Int * o.dt;
  rhs3Ext = rhs3Ext * o.dt;
end


% START TO COMPUTE RIGHT-HAND SIDE DUE TO NEW BENDING TERM
% This term is always treated fully explicitly.  The highest derivative
% it contains is second-order and it appears through a curvature
if o.bending
  f = vesicle.newBending;
  % New part of traction jump if we have variable bending
  rhs1 = rhs1 + o.dt * op.exactStokesSLdiag(vesicle,o.Galpert,f);
end
% END TO COMPUTE RIGHT-HAND SIDE DUE TO NEW BENDING TERM.  THIS WILL
% ONLY SUPPORT A SINGLE VESICLE IN AN UNBOUNDED FLOW

% Stack the right-hand sides in an alternating with respect to the
% vesicle fashion
rhs = [rhs1; rhs2];
rhs = rhs(:);
if ~o.diffDiscWalls
  rhs = [rhs; rhs3(:)];
else
  rhs = [rhs; rhs3Ext(:); rhs3Int(:)];
end
% Add on the no-slip boundary conditions on the solid walls
% Rotlet and Stokeslet equations
rhs = [rhs; zeros(3*(nvbd-1),1)];

% Use a preconditioner (block-diagonal preconditioner is implemented)
usePreco = o.usePreco;
useSpecPreco = o.useSpecPreco;
    
% START BUILDING BLOCK-DIAGONAL PRECONDITIONER
if usePreco && ~useSpecPreco
  % only build the preconditioner if updatePreco == true  
  if updatePreco

    % Build differential operators. 
    % Compute bending, tension, and surface divergence of current
    % vesicle configuration
    [Ben,Ten,Div] = vesicle.computeDerivs;

    if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
      bdiagVes.L = zeros(3*N,3*N,nv);
      bdiagVes.U = zeros(3*N,3*N,nv);
    elseif strcmp(o.solver,'method3')
      bdiagVes.GT = zeros(2*N,N,nv); % SLP of tension
      bdiagVes.DGT = zeros(N,N,nv); % divergence of SLP of tension
      bdiagVes.DGB = zeros(N,2*N,nv); % divergence of SLP of bending
      % schur complement of the lower right N by N block
      bdiagVes.schur = zeros(2*N,2*N,nv); 
    elseif strcmp(o.solver,'method4')
      bdiagVes.GT = zeros(2*N,N,nv); % SLP of tension
      % inverse of identity plus SLP of Bending
      bdiagVes.IpBen = zeros(2*N,2*N,nv); 
      bdiagVes.DGB = zeros(N,2*N,nv); % divergence of SLP of bending
      % schur complement of the upper left 2*N by 2*N block
      bdiagVes.schur = zeros(N,N,nv); 
    end
    
    % Build block-diagonal preconditioner of self-vesicle 
    % intearctions in matrix form
    for k=1:nv
      if strcmp(o.solver,'method1')
        if any(vesicle.viscCont ~= 1)
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
                o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            o.beta*Div(:,:,k) zeros(N)]);
        else
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*eye(2*N) + ...
                o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            o.beta*Div(:,:,k) zeros(N)]);
        end
      elseif strcmp(o.solver,'method2')
        if any(vesicle.viscCont ~= 1)
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu( ...
            [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
                o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            -vesicle.kappa*Div(:,:,k)*(alpha(k)*eye(2*N) - o.D(:,:,k))\...
                (o.Galpert(:,:,k)*Ben(:,:,k)) ...
            Div(:,:,k)*(alpha(k)*eye(2*N) - o.D(:,:,k))\...
                (o.Galpert(:,:,k)*Ten(:,:,k))]);
        else
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*eye(2*N) + ...
                o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt*o.Galpert(:,:,k)*Ten(:,:,k); ...
            -vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k) ...
            Div(:,:,k)*o.Galpert(:,:,k)*Ten(:,:,k)]);
        end
      elseif strcmp(o.solver,'method3')
        % schur complement of lower right block of method2
        bdiagVes.GT(:,:,k) = o.Galpert(:,:,k)*Ten(:,:,k);
        bdiagVes.DGT(:,:,k) = (Div(:,:,k)*o.Galpert(:,:,k)*...
            Ten(:,:,k))\eye(N);
        bdiagVes.DGB(:,:,k) = ...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k);
        bdiagVes.schur(:,:,k) = ...
          inv((o.beta*eye(2*N) + vesicle.kappa*o.dt*...
          o.Galpert(:,:,k)*Ben(:,:,k)) - ...
          o.dt*o.Galpert(:,:,k)*Ten(:,:,k)*bdiagVes.DGT(:,:,k)*...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k));

      elseif strcmp(o.solver,'method4')
        % schur complement of upper left block of method2
        bdiagVes.GT(:,:,k) = o.Galpert(:,:,k)*Ten(:,:,k);
        bdiagVes.IpBen(:,:,k) = inv(o.beta*eye(2*N) + ...
            vesicle.kappa*o.dt*o.Galpert(:,:,k)*Ben(:,:,k));
        bdiagVes.DGB(:,:,k) = ...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k);
        bdiagVes.schur(:,:,k) = ...
          inv(Div(:,:,k)*(-vesicle.kappa*o.dt*o.Galpert(:,:,k)*...
          Ben(:,:,k)*bdiagVes.IpBen(:,:,k) + eye(2*N))*...
          bdiagVes.GT(:,:,k));
      end
    end
    o.bdiagVes = bdiagVes;
  end % updatePreco
end % usePreco

% SOLVING THE SYSTEM USING GMRES
% Start tic-toc for GMRES iterations
if o.verbose
message = ['Using GMRES to solve...'];
o.om.writeMessage(message,'%s\n');
tGMRES = tic;
end

warning off
% any warning is printed to the terminal and the log file so
% don't need the native matlab version
initGMRES = [Xm;sigmaM];
initGMRES = initGMRES(:);
if o.confined 
  RS = RSm(:,2:end);
  if ~o.diffDiscWalls
    initGMRES = [initGMRES;etaM(:);RS(:)];
  else
    initGMRES = [initGMRES;etaMext(:);etaMint(:);RS(:)];
  end
end

% Use GMRES to solve for new positions, tension, density
% function defined on the solid walls, and rotlets/stokeslets
if usePreco && ~useSpecPreco
  [Xn,iflag,~,I,~] = gmres(@(X) o.TimeMatVec(X,vesicle,walls,...
      wallsInt,wallsExt),rhs,[],o.gmresTol,o.gmresMaxIter,...
      @o.preconditionerBD,[],initGMRES);
  iter = I(2);    

elseif useSpecPreco
  
  [Xn,iflag,~,I,~] = gmres(@(X) o.TimeMatVec(X,vesicle,walls,...
      wallsInt,wallsExt),rhs,[],o.gmresTol,o.gmresMaxIter,...
      @(z) o.preconditionerSpectral(vesicle,[],[],z),[],initGMRES);
  iter = I(2);      
else
  [Xn,iflag,~,I,~] = gmres(@(X) o.TimeMatVec(X,vesicle,walls,...
      wallsInt,wallsExt),rhs,[],o.gmresTol,o.gmresMaxIter);
  iter = I(2);
end
warning on

if o.verbose
message = ['DONE, it took ' num2str(toc(tGMRES),'%2.2e') ' seconds'];
o.om.writeMessage(message,'%s\n');
end
% END OF SOLVING THE SYSTEM USING GMRES

% allocate space for positions, tension, and density function
X = zeros(2*N,nv);
sigma = zeros(N,nv);
eta = zeros(2*Nbd,nvbdSme);
etaInt = zeros(2*NbdInt,nvbdInt);
etaExt = zeros(NbdExt,nvbdExt);
RS = zeros(3,nvbd);

% unstack the positions and tensions
for k=1:nv
  X(:,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigma(:,k) = Xn((3*k-1)*N+1:3*k*N);
end

% unstack the density function
Xn = Xn(3*nv*N+1:end);
if ~o.diffDiscWalls
  for k = 1:nvbd
    eta(:,k) = Xn((k-1)*2*Nbd+1:2*k*Nbd);
  end
else
  etaExt = Xn(1:2*NbdExt);  % assuming nvbdExt = 1
  for k = 1:nvbdInt
    etaInt(:,k) = Xn(2*NbdExt+(k-1)*2*NbdInt+1:2*NbdExt+2*k*NbdInt);
  end  
end

% unstack the rotlets and stokeslets
if ~o.diffDiscWalls
  otlets = Xn(2*nvbd*Nbd+1:end);
else
  otlets = Xn(2*NbdExt+2*nvbdInt*NbdInt+1:end);  
end

for k = 2:nvbd
  istart = (k-2)*3+1;
  iend = 3*(k-1);
  RS(:,k) = otlets(istart:iend);
end

% Compute the velocity using the differencing stencil
u = (o.beta*X - Xo)/o.dt;

% if we want to save the explicit solution, too
if o.alsoExplicit
  [XnExp,~,~,~,~] = gmres(@(X) o.TimeMatVec(X,vesicle,walls,...
      wallsInt,wallsExt),rhs,[],o.gmresTol,1,...
      @o.preconditionerBD,[],initGMRES);
  XExp = zeros(2*N,nv);
  sigmaExp = zeros(N,nv);
  
  % unstack the positions and tensions
  for k = 1 : nv
    XExp(:,k) = XnExp((3*k-3)*N+1:(3*k-1)*N);
    sigmaExp(:,k) = XnExp((3*k-1)*N+1:3*k*N);
  end
  uExp = (o.beta*XExp-Xo)/o.dt;
  
  expFile = ['/workspace/gokberk/relaxationRuns/' o.runName ...
      'ExpData_t' num2str(o.currentTime) '.mat'];
  save(expFile,'Xm','sigmaM','uM','XExp','sigmaExp','uExp','X','sigma','u')
end

end % timeStep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = TimeMatVec(o,Xn,vesicle,walls,wallsInt,wallsExt)
% val = TimeMatVec(o,Xn,vesicle,walls,wallsInt,wallsExt) 
% MATVEC for GMRES in the IMEX scheme.
% Evaluations vesicle-vesicle and vesicle-boundary interaction formulas
% to the function Xn which contains both the position, tension, and
% density
% 
% - Xn : all state variables in the following order
%   ves1:x,y,sigma,  ves2:x,y,sigma, ... vesNv:x,y,sigma, (outer_solidwall1:fx,fy,
%   inner_solidwall_1:fx,fy; inner_solid_wall2:fx,fy; ...; inner_solid_walln:fx,fy;
%   stokeslet_rotlet_innerwall1, stokeslet_rolet_innerwall2....
%
% - vesicle: class capsules used to evaluate the operators for the GMRES 
% - walls: same thing as "vesicle" but for the confined walls geometry
% - wallsInt: inner solid walls for the cases when outer and inner walls
%   are discretized with different Nbds (e.g. DLD examples)
% - wallsExt: outer solid wall for the cases mentioned above
% - either we have walls or wallsInt and wallExt together. Three of them do
%   not exist at the same time

% counter for the number of matrix-vector multiplications
% that are required for the entire simulation
global matvecs  
matvecs = matvecs + 1;

op = o.op; % poten class
N = vesicle.N; % Number of points
nv = vesicle.nv; % Number of vesicles
if o.confined
  if ~o.diffDiscWalls
    Nbd = walls.N; % Number of points on walls
    nvbdSme = walls.nv; % Number of components to walls
    NbdInt = 0;
    NbdExt = 0;
    nvbdInt = 0;
    nvbdExt = 0;
  else
    NbdInt = wallsInt.N;
    NbdExt = wallsExt.N;
    nvbdInt = wallsInt.nv;
    nvbdExt = wallsExt.nv;
    Nbd = 0;
    nvbdSme = 0;
  end
else
  Nbd = 0;
  NbdInt = 0;
  NbdExt = 0;
  nvbdSme = 0;
  nvbdInt = 0;
  nvbdExt = 0;
end
% total number of solid walls
nvbd = nvbdSme + nvbdInt + nvbdExt;

% right-hand side that corresponds to position equation
valPos = zeros(2*N,nv);
% right-hand side that corresponds to inextensibilty equation
valTen = zeros(N,nv);
% right-hand side that corresponds to solid wall equation
if o.confined 
  if ~o.diffDiscWalls
    valWalls = zeros(2*Nbd,nvbdSme);
    valWallsInt = [];
    valWallsExt = [];
  else
    valWallsInt = zeros(2*NbdInt,nvbdInt);
    valWallsExt = zeros(2*NbdExt,nvbdExt);
    valWalls = [];
  end
  % right-hand side corresponding to the rotlets and stokeslets
  valLets = zeros(3*(nvbd-1),1);
end

% Unstack the position and tension from the input
Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
for k=1:nv
  Xm(1:2*N,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigmaM(:,k) = Xn((3*k-1)*N+1:3*k*N);
end

% Unstack the density function from the input
if o.confined
    
  eta = Xn(3*nv*N+1:end);
  % x for the system DLPRSmat*x = b (W2W handled by the pre-computed matrix)
  % much faster than FMM for walls since geometries do not change. It
  % requires memory for walls with large Nbds though.
  etaAll = eta; % stacked eta  
  if ~o.diffDiscWalls
    etaM = zeros(2*Nbd,nvbd);
    for k = 1:nvbd
      etaM(:,k) = eta((k-1)*2*Nbd+1:2*k*Nbd);
    end
    otlets = Xn(3*nv*N+2*nvbd*Nbd+1:end);
    etaMint = [];
    etaMext = [];
  else
    etaMext = eta(1:2*NbdExt);
    eta = eta(2*NbdExt+1:end);
    etaMint = zeros(2*NbdInt,nvbdInt);
    for k = 1:nvbdInt
      etaMint(:,k) = eta((k-1)*2*NbdInt+1:2*k*NbdInt);
    end
    otlets = Xn(3*nv*N+2*NbdExt+2*nvbdInt*NbdInt+1:end);  
    etaM = [];
  end
else
  etaM = [];
  etaMext = [];
  etaMint = [];
  otlets = [];
end

% otlets keeps stokeslets and rotlets of each wall. Ordered as
% [stokeslet1(component 1);stokeslet1(component 2);rotlet1;...
%  stokeslet2(component 1);stokeslet2(component 2);rotlet2;...];

% f is the traction jump stored as a 2N x nv matrix
f = vesicle.tracJump(Xm,sigmaM);

% constant that multiplies the time derivative in the 
% vesicle position equation
alpha = (1+vesicle.viscCont)/2; 

% Gf is the single-layer potential applied to the traction jump. 
Gf = op.exactStokesSLdiag(vesicle,o.Galpert,f);

% DXm is the double-layer potential applied to the position
if any(vesicle.viscCont ~= 1)
  DXm = op.exactStokesDLdiag(vesicle,o.D,Xm);
else
  DXm = zeros(2*N,nv);
end


% START COMPUTING REQUIRED SINGLE-LAYER POTENTIALS
% Evaluate single-layer potential due to all vesicles except itself and
% the single-layer potential due to all vesicles evaluated on the solid
% walls.  
if ~o.fmm
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
  kernelDirect = @op.exactStokesSL;
end

% Evaulate single-layer potential due to all other vesicles
% WITH near-singular integration.  FMM is optional
SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
Fslp = op.nearSingInt(vesicle,f,SLP,o.SLPnoCorr,...
    o.NearV2V,kernel,kernelDirect,vesicle,true,false);

if o.confined
  % Evaluate single-layer potential due to all vesicles on
  % the solid walls WITH near-singular integration
  if ~o.diffDiscWalls
    FSLPwall = op.nearSingInt(vesicle,f,SLP,[],...
      o.NearV2W,kernel,kernelDirect,walls,false,false);
    FSLPwallInt = [];
    FSLPwallExt = [];
  else
    FSLPwallInt = op.nearSingInt(vesicle,f,SLP,[],...
      o.NearV2Wint,kernel,kernelDirect,wallsInt,false,false);
    FSLPwallExt = op.nearSingInt(vesicle,f,SLP,[],...
      o.NearV2Wext,kernel,kernelDirect,wallsExt,false,false);
    FSLPwall = [];
  end
else
  FSLPwall = [];
  FSLPwallInt = [];
  FSLPwallExt = [];
end
% END COMPUTING REQUIRED SINGLE-LAYER POTENTIALS


% START COMPUTING REQUIRED DOUBLE-LAYER POTENTIALS FOR VISCOSITY
% CONTRAST
if any(vesicle.viscCont ~= 1)
    
  jump = 1/2*(1-vesicle.viscCont);
  DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
  
  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
    kernelDirect = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLnewfmm;
    kernelDirect = @op.exactStokesDL;
  end
    
  % Use near-singular integration to compute double-layer
  % potential due to V2V interactions. FMM optional.
  Fdlp = op.nearSingInt(vesicle,Xm,DLP,o.DLPnoCorr,...
    o.NearV2V,kernel,kernelDirect,vesicle,true,false);
  
  if o.confined
    if ~o.diffDiscWalls    
      FDLPwall = op.nearSingInt(vesicle,Xm,DLP,[],...
        o.NearV2W,kernel,kernelDirect,walls,false,false);
      FDLPwallInt = [];
      FDLPwallExt = [];
    else
      FDLPwallInt = op.nearSingInt(vesicle,Xm,DLP,[],...
        o.NearV2Wint,kernel,kernelDirect,wallsInt,false,false);
      FDLPwallExt = op.nearSingInt(vesicle,Xm,DLP,[],...
        o.NearV2Wext,kernel,kernelDirect,wallsExt,false,false);     
      FDLPwall = [];
    end
  else
    FDLPwall = [];
    FDLPwallInt = [];
    FDLPwallExt = [];
  end
  
else
  Fdlp = [];
  FDLPwall = [];
  FDLPwallInt = [];
  FDLPwallExt = [];
end
% END COMPUTING REQUIRED DOUBLE-LAYER POTENTIALS FOR VISCOSITY CONTRAST

% START OF EVALUATING DOUBLE-LAYER POTENTIALS DUE TO SOLID WALLS ON
% VESICLES
if o.confined
  jump = -1/2;  
  if ~o.diffDiscWalls  
    potWall = o.opWall;
    if ~o.fmmDLP
      kernel = @potWall.exactStokesDL;
      kernelDirect = @potWall.exactStokesDL;
    else
      kernel = @potWall.exactStokesDLnewfmm;
      kernelDirect = @potWall.exactStokesDL;
    end
    
    DLP = @(X) jump*X + potWall.exactStokesDLdiag(walls,o.wallDLP,X);
    Fwall2Ves = potWall.nearSingInt(walls,etaM,DLP,[],...
        o.NearW2V,kernel,kernelDirect,vesicle,false,false);
  else
    potWallInt = o.opWallInt; 
    if ~o.fmmDLP
      kernel = @potWallInt.exactStokesDL;
      kernelDirect = @potWallInt.exactStokesDL;
    else
      kernel = @potWallInt.exactStokesDLnewfmm;
      kernelDirect = @potWallInt.exactStokesDL;
    end
  
    DLP = @(X) jump*X + potWallInt.exactStokesDLdiag(wallsInt,...
        o.wallDLPint,X);
    FwallInt2Ves = potWallInt.nearSingInt(wallsInt,etaMint,DLP,[],...
        o.NearWint2V,kernel,kernelDirect,vesicle,false,false);
  
    potWallExt = o.opWallExt;
    if ~o.fmmDLP
      kernel = @potWallExt.exactStokesDL;
      kernelDirect = @potWallExt.exactStokesDL;
    else
      %kernel = @potWallExt.exactStokesDL;   
      % FMM for exterior wall on vesicle makes it slow
      kernel = @potWallExt.exactStokesDLnewfmm;
      kernelDirect = @potWallExt.exactStokesDL;
    end    
    
    DLP = @(X) jump*X + potWallExt.exactStokesDLdiag(wallsExt,...
        o.wallDLPext,X);
    FwallExt2Ves = potWallExt.nearSingInt(wallsExt,etaMext,DLP,[],...
        o.NearWext2V,kernel,kernelDirect,vesicle,false,false);
  
    Fwall2Ves = FwallInt2Ves + FwallExt2Ves;
  end
else
  Fwall2Ves = zeros(2*N,nv);
end
% END OF EVALUATING DOUBLE-LAYER POTENTIALS DUE TO SOLID WALLS

% START OF EVALUATING WALL TO WALL INTERACTIONS
if o.confined
  
  % only need to do wall to wall interactions if the domain is multiply
  % connected
  if (~o.diffDiscWalls && nvbd > 1)
    if o.matFreeWalls
      potWall = o.opWall;
      if ~o.fmmDLP
        kernel = @potWall.exactStokesDL;
        FDLPwall2wall = kernel(walls,etaM,[]);
      else
        % Since the double-layer potential is still expensive even with the
        % FMM, this eliminates the need to do one big FMM followed by a
        % bunch of small ones to subtract off the self-interaction term
        % which is calculated using the precomputed matrix  
        kernel = @potWall.exactStokesDLnewfmm;
        FDLPwall2wall = kernel(walls,etaM,o.wallDLPnoCorr);
      end % o.fmmDLP
      
    else %o.matFreeWalls
      % !!in-core without fast-direct solver 
      % (out-core and with FD will be implemented) 
      
      if ~o.fastDirect % if wallDLPandRSmat is not compressed
        wallAllRHS = o.wallDLPandRSmat*etaAll;  
      else % if wallDLPandRSmat is compressed
        wallAllRHS = o.applyWall2WallFD(NbdExt,NbdInt,nvbd,etaAll);     
      end %~o.fastDirect
      FDLPwall2wall = wallAllRHS(1:2*Nbd*nvbd);
      valLets = wallAllRHS(2*Nbd*nvbd+1:end);
    end % o.matFreeWalls
    
  elseif (o.diffDiscWalls && nvbd > 1)
  % if diffDiscWalls, then nvbd > 1
    if o.matFreeWalls
      % if wall-wall interactions are computed matrix-free
      potWallInt = o.opWallInt;
      potWallExt = o.opWallExt;
      if ~o.fmmDLP
        kernel = @potWallExt.exactStokesDL;
        [~,FDLPwallExt2wallInt] = kernel(wallsExt,etaMext,[],wallsInt.X,1);
      
        kernel = @potWallInt.exactStokesDL;
        FDLPwallInt2wallInt = kernel(wallsInt,etaMint,[]);
        [~,FDLPwallInt2wallExt] = kernel(wallsInt,etaMint,[],...
            wallsExt.X,1:wallsInt.nv);
      else
        kernel = @potWallInt.exactStokesDLnewfmm;
        FDLPwallInt2wallInt = kernel(wallsInt,etaMint,o.wallDLPintNoCorr);

        [~,FDLPwallExt2wallInt] = kernel(wallsExt,etaMext,[],wallsInt.X,1);
        [~,FDLPwallInt2wallExt] = kernel(wallsInt,...
            etaMint,[],wallsExt.X,1:wallsInt.nv);
   
      end % o.fmmDLP
      
    else %matFreeWalls
    % if we want to use pre-computed wall2wall DLP and RS
      if ~o.outOfCore % do not block the operations
        if~o.HODLRforW2W % if wallDLPandRSmat is not compressed
          wallAllRHS = o.wallDLPandRSmat*etaAll;
        else % if wallDLPandRSmat is compressed with HODLR
          wallAllRHS = o.applyWall2WallFD(NbdExt,NbdInt,nvbd,etaAll);
        end
        FDLPwall2wall = wallAllRHS(1:2*NbdExt+2*NbdInt*nvbdInt);
        valLets = wallAllRHS(2*NbdExt+2*NbdInt*nvbdInt+1:end);
      else
        % block the operation
        numRows = 2*NbdExt+2*NbdInt*nvbdInt+3*nvbdInt;
        memmapWallDLPandRS = o.wallDLPandRSmat;
        wallAllRHS = zeros(numRows,1);
        bsize = min(numRows,o.maxbsize);
        for i = 1:floor(numRows/bsize)
          istart = 1+(i-1)*bsize;
          iend   = i*bsize;
          row = memmapWallDLPandRS.Data.M(istart:iend,:);
          wallAllRHS(istart:iend,1) = row*etaAll;
        end
        istart = 1+i*bsize;
        row = memmapWallDLPandRS.Data.M(istart:end,:);
        wallAllRHS(istart:end,1) = row*etaAll;
        clear row;
        FDLPwall2wall = wallAllRHS(1:2*NbdExt+2*NbdInt*nvbdInt);
        valLets = wallAllRHS(2*NbdExt+2*NbdInt*nvbdInt+1:end);
      end % ~o.outOfCore

    end % o.matFreeWalls
  elseif nvbd == 1
    if ~o.matFreeWalls
      wallAllRHS = o.wallDLPandRSmat*etaAll;  
      valWalls = wallAllRHS(1:2*Nbd*nvbd);
    end
  end % ~o.diffDiscWalls && nvbd>1 
end % o.confined
% END OF EVALUATING WALL TO WALL INTERACTIONS


% START OF EVALUATING POTENTIAL DUE TO STOKESLETS AND ROTLETS
if nvbd > 1
  if ~o.diffDiscWalls  
    LetsWalls = zeros(2*Nbd,nvbd);
    LetsVes = zeros(2*N,nv);
    for k = 2:nvbd
      stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
      rotlet = otlets(3*(k-1));
      % compute velocity due to rotlets and stokeslets on the vesicles
      LetsVes = LetsVes + o.RSlets(vesicle.X,walls.center(:,k),...
          stokeslet,rotlet);
      % compute velocity due to rotlets and stokeslets on the solid walls
      if o.matFreeWalls    
        LetsWalls = LetsWalls + o.RSlets(walls.X,walls.center(:,k),...
            stokeslet,rotlet);
      end
      % if ~matFreeWalls, these are already computed above
    end
    % Integral constraints on the density function eta related
    % to the weights of the stokeslets and rotlets
    if o.matFreeWalls
      valLets = o.letsIntegrals(otlets,etaM,etaMint,walls,wallsInt);
    end
    % if ~matFreeWalls, these are already computed above
  else
      
    LetsWallsInt = zeros(2*NbdInt,nvbdInt);
    LetsWallsExt = zeros(2*NbdExt,nvbdExt);
    LetsVes = zeros(2*N,nv);
    for k = 2:nvbd
      stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
      rotlet = otlets(3*(k-1));
      % compute velocity due to rotlets and stokeslets on the vesicles
      LetsVes = LetsVes + o.RSlets(vesicle.X,wallsInt.center(:,k-1),...
          stokeslet,rotlet);

      % compute velocity due to rotlets and stokeslets on the solid walls
      if o.matFreeWalls
        LetsWallsInt = LetsWallsInt + o.RSlets(wallsInt.X,...
            wallsInt.center(:,k-1),stokeslet,rotlet);
        LetsWallsExt = LetsWallsExt + o.RSlets(wallsExt.X,...
            wallsInt.center(:,k-1),stokeslet,rotlet);
      end
    end
    % Integral constraints on the density function eta related
    % to the weights of the stokeslets and rotlets
    if o.matFreeWalls
      valLets = o.letsIntegrals(otlets,etaM,etaMint,walls,wallsInt);
    end      
  end
else
  LetsVes = [];
  LetsWalls = [];
  FDLPwall2wall = [];
  LetsWallsInt = [];
  LetsWallsExt = [];
  FDLPwallInt2wallInt = [];
  FDLPwallExt2wallInt = [];
  FDLPwallInt2wallExt = [];
end
% END OF EVALUATING POTENTIAL DUE TO STOKESLETS AND ROTLETS

% START OF EVALUATING VELOCITY ON VESICLES

if ~isempty(Gf)
  % self-bending and self-tension terms
  valPos = valPos - o.dt*Gf*diag(1./alpha);
end

if ~isempty(DXm)
  % self-viscosity contrast term
  valPos = valPos - o.beta*DXm*diag(1./alpha);
end
if ~isempty(Fslp)
  % single-layer potential due to all other vesicles
  valPos = valPos - o.dt*Fslp*diag(1./alpha);
end
if ~isempty(Fdlp)
  % double-layer potential due to all other vesicles
  valPos = valPos - o.beta*Fdlp*diag(1./alpha);
end

if o.confined
  if ~isempty(Fwall2Ves)    
    % velocity due to solid walls evaluated on vesicles 
    valPos = valPos - o.dt*Fwall2Ves*diag(1./alpha);
  end
  if ~isempty(LetsVes)
    % velocity on vesicles due to the rotlets and stokeslets
    valPos = valPos - o.dt*LetsVes*diag(1./alpha);
  end
end
% END OF EVALUATING VELOCITY ON VESICLES

% START OF EVALUATING VELOCITY ON WALLS
% evaluate velocity on solid walls due to the density function.
% self solid wall interaction
if o.confined
  if ~o.diffDiscWalls
    if o.matFreeWalls
      potWall = o.opWall;
      valWalls = valWalls - 1/2*etaM + ...
        potWall.exactStokesDLdiag(walls,o.wallDLP,etaM);
      valWalls(:,1) = valWalls(:,1) + ...
        potWall.exactStokesN0diag(walls,o.wallN0,etaM(:,1));
    end
  else
    if o.matFreeWalls
      potWallInt = o.opWallInt;
      valWallsInt = valWallsInt - 1/2*etaMint + ...
        potWallInt.exactStokesDLdiag(wallsInt,o.wallDLPint,etaMint);
  
      potWallExt = o.opWallExt;
      valWallsExt = valWallsExt - 1/2*etaMext + ...
        potWallExt.exactStokesDLdiag(wallsExt,o.wallDLPext,etaMext);
      valWallsExt(:,1) = valWallsExt(:,1) + ...
        potWallExt.exactStokesN0diag(wallsExt,o.wallN0,etaMext);    
    end
  end
end

if o.confined
  if ~o.diffDiscWalls
    if ~isempty(FSLPwall)  
      % velocity on walls due to the vesicle traction jump
      valWalls = valWalls + FSLPwall;
    end
    if ~isempty(FDLPwall)
      % velocity on walls due to the vesicle viscosity jump
      valWalls = valWalls + o.beta*FDLPwall/o.dt;
    end
    if o.matFreeWalls
      if ~isempty(FDLPwall2wall)  
        % velocity on walls due to all other walls
        valWalls = valWalls + FDLPwall2wall;    
      end
      if ~isempty(LetsWalls)
        % velocity on walls due to the rotlets and stokeslets
        valWalls = valWalls + LetsWalls;  
      end
    else
      if ~isempty(FDLPwall2wall) 
        valWalls = valWalls(:) + FDLPwall2wall;  
      end
    end

  else
    if ~isempty(FSLPwallInt)  
      % velocity on walls due to the vesicle traction jump  
      valWallsInt = valWallsInt + FSLPwallInt;
    end

    if ~isempty(FDLPwallInt)  
      % velocity on walls due to the vesicle viscosity jump
      valWallsInt = valWallsInt + o.beta*FDLPwallInt/o.dt;
    end
    
    if ~isempty(FSLPwallExt)  
      % velocity on walls due to the vesicle traction jump
      valWallsExt = valWallsExt + FSLPwallExt;
    end
    
    if ~isempty(FDLPwallExt)  
      % velocity on walls due to the vesicle viscosity jump  
      valWallsExt = valWallsExt + o.beta*FDLPwallExt/o.dt;
    end
    if o.matFreeWalls
      if ~isempty(FDLPwallExt2wallInt)  
        % velocity on walls due to all other walls
        valWallsInt = valWallsInt + FDLPwallExt2wallInt;
      end
      if ~isempty(FDLPwallInt2wallInt)  
        % velocity on walls due to all other walls
        valWallsInt = valWallsInt + FDLPwallInt2wallInt;
      end
      if ~isempty(LetsWallsInt) 
        % velocity on walls due to the rotlets and stokeslets  
        valWallsInt = valWallsInt + LetsWallsInt;
      end
      if ~isempty(FDLPwallInt2wallExt) 
        % velocity on walls due to all other walls
        valWallsExt = valWallsExt + FDLPwallInt2wallExt;
      end
      if ~isempty(LetsWallsExt) 
        % velocity on walls due to the rotlets and stokeslets
        valWallsExt = valWallsExt + LetsWallsExt;
      end
      % pad the valWalls
      valWalls = [valWallsExt(:);valWallsInt(:)];
    else
      valWalls = [valWallsExt(:);valWallsInt(:)];  
      if ~isempty(FDLPwall2wall)
        valWalls = valWalls + FDLPwall2wall;  
      end
    end % if o.matFreeWalls
  end %if ~o.diffDiscWalls
end % if o.confined
% END OF EVALUATING VELOCITY ON WALLS

% START OF EVALUATING INEXTENSIBILITY CONDITION
% Two possible discretizations of the inextensibility condition
if (strcmp(o.solver,'method1'))
  % compute surface divergence of the current GMRES iterate
  % method1 sets this equal to the surface divergence of
  % the previous time step
  valTen = o.beta * vesicle.surfaceDiv(Xm);

else
  if any(vesicle.viscCont ~= 1)
    if o.confined
      valTen = vesicle.surfaceDiv(...
        o.solveIminusD(Gf+Fslp+Fwall2Ves+LetsVes,vesicle));
    else
      valTen = vesicle.surfaceDiv(...
        o.solveIminusD(Gf+Fslp,vesicle));  
    end
    
  else
    valTen = -1/o.dt*vesicle.surfaceDiv(valPos);
  end
  % method2, method3, and method4 sets the surface divergence of the sum
  % of single-layer potentials due to bending and tension plus the
  % farField to zero.  The only difference between the two methods is
  % how the preconditioner is used.  method2 uses the possibly
  % ill-conditioned full matrix where as method3 and method 4 use the
  % two schur complements.  Eventually will phase out method2 and it
  % will be fully replaced by method3 and method4
end
% END OF EVALUATING INEXTENSIBILITY CONDITION

% beta times solution coming from time derivative
valPos = valPos + o.beta*Xm;

% Initialize output from vesicle and inextensibility equations to zero
val = zeros(3*N*nv,1);

% Stack val as [x-coordinate;ycoordinate;tension] repeated
% nv times for each vesicle
for k=1:nv
  val((k-1)*3*N+1:3*k*N) = [valPos(:,k);valTen(:,k)];
end
if (any(vesicle.viscCont ~= 1) && o.confined)
  % This combination of options causes problems with
  % the scaling of the preconditioner.  Need to
  % get rid of the potentially small value o.dt
  valWalls = valWalls * o.dt;
end

% Stack velocity along the solid walls in same manner as above
% Stack the stokeslets and rotlet componenets at the end
if o.confined
  val = [val;valWalls(:);valLets];
end
  
end % TimeMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = TimeMatVecLowAcc(o,Xn,vesicle)
% val = TimeMatVecLowAcc(o,Xn,vesicle) 
% A version of TimeMatVec that does not use a special quadrature for self
% interactions and near singular integration for near interactions
% 
% - Xn : all state variables in the following order
%   ves1:x,y,sigma,  ves2:x,y,sigma, ... vesNv:x,y,sigma, 
%
% - vesicle: class capsules used to evaluate the operators for the GMRES 

% counter for the number of matrix-vector multiplications
% that are required for the entire simulation
global matvecs  
matvecs = matvecs + 1;

op = o.op; % poten class
N = vesicle.N; % Number of points
nv = vesicle.nv; % Number of vesicles

% right-hand side that corresponds to position equation
valPos = zeros(2*N,nv);

% Unstack the position and tension from the input
Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
for k=1:nv
  Xm(1:2*N,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigmaM(:,k) = Xn((3*k-1)*N+1:3*k*N);
end

% f is the traction jump stored as a 2N x nv matrix
f = vesicle.tracJump(Xm,sigmaM);

% constant that multiplies the time derivative in the 
% vesicle position equation
alpha = (1+vesicle.viscCont)/2; 

% Gf is the single-layer potential applied to the traction jump. 
Gf = op.exactStokesSLdiag(vesicle,o.SLPnoCorr,f);


% DXm is the double-layer potential applied to the position
if any(vesicle.viscCont ~= 1)
  DXm = op.exactStokesDLdiag(vesicle,o.D,Xm);
else
  DXm = [];
end

% START COMPUTING REQUIRED SINGLE-LAYER POTENTIALS
% Evaluate single-layer potential due to all vesicles except itself and
% the single-layer potential due to all vesicles evaluated on the solid
% walls.  

% Evaulate single-layer potential due to all other vesicles
if ~o.fmm
  Fslp = op.exactStokesSL(vesicle,f,o.SLPnoCorr);
else
  Fslp = op.exactStokesSLfmm(vesicle,f,o.SLPnoCorr);
end


% START COMPUTING REQUIRED DOUBLE-LAYER POTENTIALS FOR VISCOSITY
% CONTRAST
if any(vesicle.viscCont ~= 1)
    
  % compute double-layer potential due to V2V interactions.
  if ~o.fmmDLP
    Fdlp = op.exactStokesDL(vesicle,Xm,o.DLPnoCorr);
  else
    Fdlp = op.exactStokesDLnewfmm(vesicle,Xm,o.DLPnoCorr);  
  end
  
else
  Fdlp = [];
end
% END COMPUTING REQUIRED DOUBLE-LAYER POTENTIALS FOR VISCOSITY CONTRAST

% START OF EVALUATING VELOCITY ON VESICLES

% self-bending and self-tension terms
valPos = valPos - o.dt*Gf*diag(1./alpha);

if ~isempty(DXm)
  % self-viscosity contrast term
  valPos = valPos - o.beta*DXm*diag(1./alpha);
end

% single-layer potential due to all other vesicles
valPos = valPos - o.dt*Fslp*diag(1./alpha);

if ~isempty(Fdlp)
  % double-layer potential due to all other vesicles
  valPos = valPos - o.beta*Fdlp*diag(1./alpha);
end
% END OF EVALUATING VELOCITY ON VESICLES

% START OF EVALUATING INEXTENSIBILITY CONDITION
% Two possible discretizations of the inextensibility condition
if (strcmp(o.solver,'method1'))
  % compute surface divergence of the current GMRES iterate
  % method1 sets this equal to the surface divergence of
  % the previous time step
  valTen = o.beta * vesicle.surfaceDiv(Xm);
else
  if any(vesicle.viscCont ~= 1)
    valTen = vesicle.surfaceDiv(o.solveIminusD(Gf+Fslp,vesicle));  
  else
    valTen = -1/o.dt*vesicle.surfaceDiv(valPos);
  end
  % method2, method3, and method4 sets the surface divergence of the sum
  % of single-layer potentials due to bending and tension plus the
  % farField to zero.  The only difference between the two methods is
  % how the preconditioner is used.  method2 uses the possibly
  % ill-conditioned full matrix where as method3 and method 4 use the
  % two schur complements.  Eventually will phase out method2 and it
  % will be fully replaced by method3 and method4
end
% END OF EVALUATING INEXTENSIBILITY CONDITION

% beta times solution coming from time derivative
valPos = valPos + o.beta*Xm;

% Initialize output from vesicle and inextensibility equations to zero
val = zeros(3*N*nv,1);

% Stack val as [x-coordinate;ycoordinate;tension] repeated
% nv times for each vesicle
for k=1:nv
  val((k-1)*3*N+1:3*k*N) = [valPos(:,k);valTen(:,k)];
end
  
end % TimeMatVecLowAcc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vesVel,divVesVel,wallVel,wallIntVel,wallExtVel,residual] = ...
      computeVelocities(o,vesicle,Galpert,Gbarnett,D,walls,wallsInt,...
      wallsExt,etaProv,etaIntProv,etaExtProv,RSprov,NearV2V,NearV2W,...
      NearW2V,NearW2W,NearV2Wint,NearWint2V,NearV2Wext,NearWext2V,...
      NearWint2Wext,NearWint2Wint,NearWext2Wint)
% [vesVel,divVesVel,wallVel,wallIntVel,wallExtVel,residual] = ...
%       computeVelocities(o,vesicle,Galpert,Gbarnett,D,walls,wallsInt,...
%       wallsExt,etaProv,etaIntProv,etaExtProv,RSprov,NearV2V,NearV2W,...
%       NearW2V,NearW2W,NearV2Wint,NearWint2V,NearV2Wext,NearWext2V,...
%       NearWint2Wext,NearWint2Wint,NearWext2Wint)      
% computes the velocity induced by the provisional vesicle position,
% tension, and density function on the vesicles and the walls.  Also
% returns the vesicle divergence of the velocity field and the residual
% of the picard integral formulation of the vesicle velocity.  These
% quantities are needed to form the modificiations of the right-hand
% side when doing SDC updates

N = vesicle.N;
nv = vesicle.nv;
if o.confined
  if ~o.diffDiscWalls  
    Nbd = walls.N;
    nvbd = walls.nv;
    NbdInt = 0;
    NbdExt = 0;
    nvbdInt = 0;
    nvbdExt = 0;
  else
    NbdInt = wallsInt.N;
    NbdExt = wallsExt.N;
    nvbdInt = wallsInt.nv;
    nvbdExt = wallsExt.nv;
    nvbd = nvbdInt + nvbdExt;
  end
else
  Nbd = 0;
  nvbd = 0;
  NbdInt = 0;
  NbdExt = 0;
  nvbdInt = 0;
  nvbdExt = 0;
end

% velocity on the vesicles due to the provisional solution 
vesVel = zeros(2*N,nv,abs(o.orderGL));
% vesicle divergence of the vesicles due to the provisional
% solution
divVesVel = zeros(N,nv,abs(o.orderGL));
% velocity on the solid walls due to the provisional solution
if o.confined
  if ~o.diffDiscWalls
    wallVel = zeros(2*Nbd,nvbd,abs(o.orderGL));
    wallIntVel = [];
    wallExtVel = [];
  else
    wallIntVel = zeros(2*NbdInt,nvbdInt,abs(o.orderGL));
    wallExtVel = zeros(2*NbdExt,nvbdExt,abs(o.orderGL));
    wallVel = [];
  end
else
  wallVel = [];
  wallIntVel = [];
  wallExtVel = [];
end

% residual of the Picard integral coming from the time
% derivative term of the vesicle position
residual = zeros(2*N,nv,abs(o.orderGL));


% need to use implicit so that other vesicles are used to compute
% the integrand z
vesves = o.vesves;
o.vesves = 'implicit';

% need to save the time stepping order
order = o.order;
o.order = 1;

[o.Xcoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);

% need to save current time step
dt = o.dt;

% avoid introducing numerical error by multiplying and dividing
% by a potentially small number
o.dt = 1;

if ~o.confined
  z = zeros(3*N*nv,1);
else
  if ~o.diffDiscWalls  
    z = zeros(3*N*nv + 2*Nbd*nvbd + 3*(nvbd-1),1);
  else
    z = zeros(3*N*nv + 2*NbdExt*nvbdExt + 2*NbdInt*nvbdInt + 3*(nvbd-1),1);  
  end
end

op = o.op;
for n = 1:abs(o.orderGL)
  % loop over the configurations
  o.Galpert = Galpert(:,:,:,n);

  % Already computed the single-layer potential when forming the
  % provisional solution Store the single- and double-layer potentials
  % in the object o
  if any(vesicle(n).viscCont ~= 1)
    o.D = D(:,:,:,n);
  end
  
  % retrieve near-singular integration strucutre at Gauss-Lobatto time
  % step n
  o.NearV2V = NearV2V{n};
  if o.confined
    o.NearV2W = NearV2W{n};
    o.NearW2V = NearW2V{n};
    o.NearV2Wint = NearV2Wint{n};
    o.NearWint2V = NearWint2V{n};
    o.NearV2Wext = NearV2Wext{n};
    o.NearWext2V = NearWext2V{n};
    
    if isempty(o.NearW2W)
      o.NearW2W = NearW2W{n};
    end
    if isempty(o.NearWint2Wext)
      o.NearWint2Wext = NearWint2Wext{n};
    end
    if isempty(o.NearWint2Wint)
      o.NearWint2Wint = NearWint2Wint{n};
    end
    if isempty(o.NearWext2Wint)
      o.NearWext2Wint = NearWext2Wint{n};
    end
  end
  
  % put positions and tension into z the way that TimeMatVec requires
  for k = 1:nv
    z(3*(k-1)*N+1:3*k*N) = [vesicle(n).X(:,k);vesicle(n).sig(:,k)];
  end
  
  % put in density function on solid walls into z the way that
  % TimeMatVec requires them
  if o.confined
    if ~o.diffDiscWalls
      for k = 1:nvbd
        istart = 3*nv*N + 2*(k-1)*Nbd + 1;
        iend = istart + 2*Nbd - 1;
        z(istart:iend) = etaProv(:,k,n);
      end
      
      % put the rotlet and stokeslet coefficients into z
      for k = 2:nvbd
        istart = 3*nv*N+2*nvbd*Nbd+3*(k-2)+1;
        iend = istart + 2;
        z(istart:iend) = RSprov(:,k,n);
      end
    else
      % assuming only one exterior wall
      istart = 3*nv*N+1;
      iend = istart + 2*NbdExt - 1;
      z(istart:iend) = etaExtProv(:,1,n);
      
      for k = 1:nvbdInt
        istart = 3*nv*N + 2*NbdExt + 2*(k-1)*NbdInt + 1;
        iend = istart + 2*NbdInt - 1;
        z(istart:iend) = etaIntProv(:,k,n);
      end
      
      % put the rotlet and stokeslet coefficients into z
      for k = 2:nvbdInt
        istart = 3*nv*N+2*NbdExt+2*nvbdInt*NbdInt+3*(k-2)+1;
        iend = istart + 2;
        z(istart:iend) = RSprov(:,k,n);
      end
    end
  end

  % don't want the double-layer contribution when computing the velocity
  % of the vesicle
  viscCont = vesicle(n).viscCont;
  vesicle(n).viscCont = ones(1,nv);
  
  % use TimeMatVec to find velocity on vesicles and solid
  % walls due to the provisional solution
  z = o.TimeMatVec(z,vesicle(n),walls);
  
  % set the viscosity contrast back to its original value
  vesicle(n).viscCont = viscCont;
  

  % form the velocity on the vesicle due to the current provisional
  % solution
  rhs = zeros(2*N,nv);
  for k = 1:nv
    istart = 3*(k-1)*N+1;
    iend = istart + 2*N - 1;
    rhs(:,k) = -z(istart:iend);
  end
  if ~o.confined
    rhs = rhs + o.farField(vesicle(n).X,[]);
  end
  vesVel(:,:,n) = vesicle(n).X + rhs;
  

  % need to apply inv(alpha*I - DLP) if there is a viscosity contrast to
  % obtain the velocity of the vesicles.
  if any(vesicle(n).viscCont ~= 1)
    vesVel(:,:,n) = o.solveIminusD(vesVel(:,:,n),vesicle(n));
  end 
  
  % surface divergence of the velocity of the vesicle due
  % to the provisional solution
  divVesVel(:,:,n) = vesicle(n).surfaceDiv(vesVel(:,:,n));

  if o.confined
    if ~o.diffDiscWalls
      for k = 1:nvbd
        istart = 3*nv*N + (k-1)*2*Nbd + 1;
        iend = istart + 2*Nbd - 1;
        wallVel(:,k,n) = z(istart:iend);
      end
    else
      istart = 3*nv*N + 1;
      iend = istart + 2*NbdExt - 1;
      wallExtVel(:,1,n) = z(istart:iend);
      for k = 1:nvbdInt
        istart = 3*nv*N + 2*NbdExt + (k-1)*2*NbdInt + 1;
        iend = istart + 2*NbdInt - 1;
        wallIntVel(:,k,n) = z(istart:iend);
      end
    end % ~o.diffDiscWalls
  end % ~o.confined

  % velocity on the solid walls due to the provisional solution
  if any(vesicle(n).viscCont ~= 1) && o.confined
    jump = 1/2*(1-vesicle(n).viscCont);
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle(n),o.D,X);
    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
      kernelDirect = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLnewfmm;
      kernelDirect = @op.exactStokesDL;
    end

    den = vesVel(:,:,n);
    % Evaulate the velocity due to the viscosity contrast
    % on the walls WITH near-singulation integration
    if ~o.diffDiscWalls
      FDLPwall = op.nearSingInt(vesicle(n),den,DLP,[],...
          o.NearV2W,kernel,kernelDirect,walls,false,false);
      wallVel(:,:,n) = wallVel(:,:,n) + FDLPwall;
    else
      FDLPwallInt = op.nearSingInt(vesicle(n),den,DLP,[],...
          o.NearV2Wint,kernel,kernelDirect,wallsInt,false,false);
      wallIntVel(:,:,n) = wallIntVel(:,:,n) + FDLPwallInt;
      
      FDLPwallExt = op.nearSingInt(vesicle(n),den,DLP,[],...
          o.NearV2Wext,kernel,kernelDirect,wallsExt,false,false);
      wallExtVel(:,:,n) = wallExtVel(:,:,n) + FDLPwallExt;
    end
  end
end % for o.orderGL

% change back to original time step and vesicle-vesicle interaction
o.dt = dt;
o.order = order;
[o.Xcoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);
o.vesves = vesves;

% integrate the vesicles velocity using quadrature rules 
% that are exact for polynomials defined at the 
% Gauss-Lobatto points
IvesVel = o.lobattoInt(vesVel);

% compute residual by adding the initial vesicle configuartion and
% subtracting the current vesicle configuartion
for n = 1:abs(o.orderGL)
  residual(:,:,n) = vesicle(1).X - vesicle(n).X + ...
      o.dt/2 * IvesVel(:,:,n);
end


end % computeVelocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,sigma,eta,RS,iter,iflag] = timeStepSimple(o,...
    Xstore,sigStore,etaStore,RSstore,viscCont,walls,vesicle)

N = size(Xstore,1)/2; % Number of points per vesicle
nv = size(Xstore,2); % Number of vesicles
if o.confined  
Xwalls = walls.X; % discretization points of solid walls
else
  Xwalls = [];
end

Nbd = size(Xwalls,1)/2; % Number of points on the solid walls

% number of solid wall components
nvbd = size(Xwalls,2);    % # of walls of the same discretization

% constant that appears in front of time derivative in
% vesicle dynamical equations
alpha = (1 + viscCont)/2; 

% Form linear combinations of previous time steps needed for Ascher,
% Ruuth, and Wetton IMEX methods
etaM = zeros(2*Nbd,nvbd);
RSm = zeros(3,nvbd);
Xm = Xstore;
sigmaM = sigStore;
if o.confined
  etaM = etaStore;
  RSm = RSstore;
end
Xo = Xstore;


% Build single layer potential matrix and put it in current object
% If we are doing an sdc update, this is already precomputed and 
% stored from when we formed the provisional solution
op = o.op;
o.Galpert = op.stokesSLmatrix(vesicle);


% Compute double-layer potential matrix due to each vesicle
% independent of the others.  Matrix is zero if there is no
% viscosity contrast
if any(viscCont ~= 1)
  o.D = op.stokesDLmatrix(vesicle);
  if o.fmmDLP
    Nup = op.LPuprate * vesicle.N;
    Xup = [interpft(vesicle.X(1:end/2,:),Nup);...
      interpft(vesicle.X(end/2+1:end,:),Nup)];
    vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0);
    o.DLPnoCorr = op.stokesDLmatrixNoCorr(vesicleUp);
  else
    o.DLPnoCorr = [];
  end
else
  o.D = [];
  o.DLPnoCorr = [];
end

if o.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(vesicle.X(1:end/2,:),Nup);...
    interpft(vesicle.X(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0); 
  o.SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);
else
  o.SLPnoCorr = [];    
end

% Only form near-singular integration structure 

% Structures for deciding who
% is close, how close it is, who is closest, etc., needed in nearSingInt
if o.confined
  % Need vesicle to vesicle and vesicle to wall interactions
  [o.NearV2V,o.NearV2W] = vesicle.getZone(walls,3);

  % Only need wall to vesicle interactions.  Wall to wall
  % interactions should also use near-singular integration since
  % they may be close to one another
  if nvbd == 1
    [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,2);
  else
    if isempty(o.NearW2W)
      [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,3);
    else
      % there is no need to compute W2W again, since they do not move  
      [~,o.NearW2V] = walls.getZone(vesicle,2);
    end
  end
else
% no solid walls, so only need vesicle-vesicle intearactions
if nv > 1
o.NearV2V = vesicle.getZone([],1);
else
o.NearV2V = [];
end

o.NearV2W = [];
o.NearW2V = [];
o.NearW2W = [];
end



% Parts of rhs from previous solution.  The right-hand-side depends on
% whether we are doing an SDC correction or forming the provisional
% solution.

rhs1 = Xo;
rhs2 = zeros(N,nv);
if o.confined
  rhs3 = walls.u;
else
  rhs3 = [];
end

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP
% vesicle-vesicle and vesicle-wall interactions are handled
% implicitly in TimeMatVec
% END TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VISCOSITY CONTRAST
if any(vesicle.viscCont ~= 1)
  
  % Need to add jump to double-layer potential if using near-singular
  % integration so that we can compute boundary values for the
  % near-singular integration algorithm
  jump = 1/2*(1-vesicle.viscCont);
  DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);

  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
    kernelDirect = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLnewfmm;
    kernelDirect = @op.exactStokesDL;
  end
  
  
  density = Xo;
  % Use near-singular integration to compute double-layer
  % potential from previous solution 
  Fdlp = op.nearSingInt(vesicle,density,DLP,o.DLPnoCorr,...
       o.NearV2V,kernel,kernelDirect,vesicle,true,false);
      
  if o.confined
    FDLPwall = op.nearSingInt(vesicle,density,DLP,[],...
      o.NearV2W,kernel,kernelDirect,walls,false,false);
  else
    FDLPwall = [];
  end
    
else
  % If no viscosity contrast, there is no velocity induced due to a
  % viscosity contrast  
  Fdlp = zeros(2*N,nv);
  FDLPwall = zeros(2*Nbd,nvbd);
end

% add in viscosity contrast term due to each vesicle independent of the
% others (o.D * Xo) from the previous solution followed by the term due
% to all other vesicles (Fdlp)
if (any(viscCont ~= 1) && ~o.SDCcorrect)
  DXo = op.exactStokesDLdiag(vesicle,o.D,Xo);
  rhs1 = rhs1 - (Fdlp + DXo) * diag(1./alpha);
end

% compute the double-layer potential due to all other vesicles from the
% appropriate linear combination of previous time steps.  Depends on
% time stepping order and vesicle-vesicle discretization
rhs3 = rhs3 + FDLPwall/o.dt;


% START COMPUTING SINGLE-LAYER POTENTIAL FOR REPULSION
if o.repulsion 
  % Repulsion is handled explicitly between vesicles, vesicles-walls.  
  if ~o.fmm
      kernel = @op.exactStokesSL;
      kernelDirect = @op.exactStokesSL;
  else
      kernel = @op.exactStokesSLfmm;
      kernelDirect = @op.exactStokesSL;
  end
 
  Xrep = Xo;
  if ~o.confined
    repulsion = vesicle.repulsionScheme(Xrep,o.repStrength,o.minDist,...
        [],[],[]);
  else
    repulsion = vesicle.repulsionScheme(Xrep,o.repStrength,o.minDist,...
        walls,[],[]);
  end

  Frepulsion = op.exactStokesSLdiag(vesicle,o.Galpert,repulsion);
  % diagonal term of repulsion


  SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
    
  % Use near-singular integration to compute single-layer potential
  % due to all other vesicles.  Need to pass function
  % op.exactStokesSL so that the layer potential can be computed at
  % far points and Lagrange interpolation points
  Frepulsion = Frepulsion + ...
    op.nearSingInt(vesicle,repulsion,SLP,o.SLPnoCorr,...
        o.NearV2V,kernel,kernelDirect,vesicle,true,false);
    
  % Evaluate the velocity on the walls due to the vesicles
  if o.confined
    FREPwall = op.nearSingInt(vesicle,repulsion,SLP,[],...
        o.NearV2W,kernel,kernelDirect,walls,false,false);      
  else
    FREPwall = [];
  end
    
  % keep the number of timesteps when repulsion is nonzero
  if norm(Frepulsion) ~= 0 
    global repuls
    repuls = repuls + 1;
  end
    
  rhs1 = rhs1 + o.dt*Frepulsion*diag(1./alpha);
  rhs3 = rhs3 - FREPwall;
end
% END COMPUTING SINGLE-LAYER POTENTIALS FOR REPULSION 



% START TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS
% This is done implicitly in TimeMatVec
if ~o.confined
  % Add in far-field condition (extensional, shear, etc.)
  vInf = o.farField(Xm);
  rhs1 = rhs1 + o.dt*vInf*diag(1./alpha);
end
% END TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS


% START TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION
% Vesicle-vesicle interactions and the presence of a
% viscosity contrast is irrelavent  
% rhs2 is the right-hand side for the inextensibility condition
rhs2 = rhs2 + vesicle.surfaceDiv(Xo); 
% END TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION

% The next makes sure that the rhs is all order one rather than have rhs3
% being order 1/o.dt and other two parts (rhs1 and rhs2) being order 1.
% This of course needs to be compensated in the TimeMatVec routine
if (any(vesicle.viscCont ~= 1) && ...
      strcmp(o.vesves,'implicit') && o.confined)
  rhs3 = rhs3 * o.dt;
end

% Stack the right-hand sides in an alternating with respect to the
% vesicle fashion
rhs = [rhs1; rhs2];
rhs = rhs(:);
rhs = [rhs; rhs3(:)];
% Add on the no-slip boundary conditions on the solid walls
% Rotlet and Stokeslet equations
rhs = [rhs; zeros(3*(nvbd-1),1)];

% Use a preconditioner (block-diagonal preconditioner is implemented)
usePreco = o.usePreco;
useSpecPreco = o.useSpecPreco;
    
% START BUILDING BLOCK-DIAGONAL PRECONDITIONER
if usePreco && ~useSpecPreco
  
  % Build differential operators. 
  % Compute bending, tension, and surface divergence of current
  % vesicle configuration
  [Ben,Ten,Div] = vesicle.computeDerivs;
  bdiagVes.L = zeros(3*N,3*N,nv);
  bdiagVes.U = zeros(3*N,3*N,nv);

  % Build block-diagonal preconditioner of self-vesicle 
  % intearctions in matrix form
  for k=1:nv
    if any(vesicle.viscCont ~= 1)
      [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
        [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
            o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
        -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
        o.beta*Div(:,:,k) zeros(N)]);
    else
      [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
        [o.beta*eye(2*N) + ...
            o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
        -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
        o.beta*Div(:,:,k) zeros(N)]);
    end
  end
  o.bdiagVes = bdiagVes;
end % usePreco

% SOLVING THE SYSTEM USING GMRES
warning off
% any warning is printed to the terminal and the log file so
% don't need the native matlab version
initGMRES = [Xm;sigmaM];
initGMRES = initGMRES(:);
if o.confined 
  RS = RSm(:,2:end);
  initGMRES = [initGMRES;etaM(:);RS(:)];
end

% Use GMRES to solve for new positions, tension, density
% function defined on the solid walls, and rotlets/stokeslets
if usePreco && ~useSpecPreco
  [Xn,iflag,~,I,~] = gmres(@(X) o.TimeMatVec(X,vesicle,walls,...
      [],[]),rhs,[],o.gmresTol,o.gmresMaxIter,...
      @o.preconditionerBD,[],initGMRES);
  iter = I(2);    

elseif useSpecPreco
  
  [Xn,iflag,~,I,~] = gmres(@(X) o.TimeMatVec(X,vesicle,walls,...
      [],[]),rhs,[],o.gmresTol,o.gmresMaxIter,...
      @(z) o.preconditionerSpectral(vesicle,[],[],z),[],initGMRES);
  iter = I(2);      
else
  [Xn,iflag,~,I,~] = gmres(@(X) o.TimeMatVec(X,vesicle,walls,...
      [],[]),rhs,[],o.gmresTol,o.gmresMaxIter);
  iter = I(2);
end
warning on
% END OF SOLVING THE SYSTEM USING GMRES

% allocate space for positions, tension, and density function
X = zeros(2*N,nv);
sigma = zeros(N,nv);
eta = zeros(2*Nbd,nvbd);
RS = zeros(3,nvbd);

% unstack the positions and tensions
for k=1:nv
  X(:,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigma(:,k) = Xn((3*k-1)*N+1:3*k*N);
end

% unstack the density function
Xn = Xn(3*nv*N+1:end);
for k = 1:nvbd
  eta(:,k) = Xn((k-1)*2*Nbd+1:2*k*Nbd);
end


% unstack the rotlets and stokeslets
otlets = Xn(2*nvbd*Nbd+1:end);
for k = 2:nvbd
  istart = (k-2)*3+1;
  iend = 3*(k-1);
  RS(:,k) = otlets(istart:iend);
end

end % timeStepSimple

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,sigma,etaInt,etaExt,RS,iter,iflag] = timeStepSimpleDiffDisc(o,...
    Xstore,sigStore,etaIntStore,etaExtStore,RSstore,viscCont,...
    wallsInt,wallsExt,vesicle)

N = size(Xstore,1)/2; % Number of points per vesicle
nv = size(Xstore,2); % Number of vesicles
XwallsInt = wallsInt.X; % interior walls
XwallsExt = wallsExt.X;
Xwalls = []; % discretization points of solid walls

Nbd = size(Xwalls,1)/2; % Number of points on the solid walls
NbdInt = size(XwallsInt,1)/2; % Number of points on the interior walls
NbdExt = size(XwallsExt,1)/2; % Number of points on the exterior walls

% number of solid wall components
nvbdInt = size(XwallsInt,2); % # of interior walls
nvbdExt = size(XwallsExt,2); % # of exterior walls
nvbdSme = size(Xwalls,2);    % # of walls of the same discretization
nvbd = nvbdSme + nvbdInt + nvbdExt;

% constant that appears in front of time derivative in
% vesicle dynamical equations
alpha = (1 + viscCont)/2; 

% Form linear combinations of previous time steps needed for Ascher,
% Ruuth, and Wetton IMEX methods
etaMint = zeros(2*NbdInt,nvbdInt);
etaMext = zeros(2*NbdExt,nvbdExt);
RSm = zeros(3,nvbd);
Xm = Xstore;
sigmaM = sigStore;
if o.confined
  etaMint = etaIntStore;
  etaMext = etaExtStore;
  RSm = RSstore;
end
Xo = Xstore;


% Build single layer potential matrix and put it in current object
% If we are doing an sdc update, this is already precomputed and 
% stored from when we formed the provisional solution
op = o.op;
o.Galpert = op.stokesSLmatrix(vesicle);


% Compute double-layer potential matrix due to each vesicle
% independent of the others.  Matrix is zero if there is no
% viscosity contrast
if any(viscCont ~= 1)
  o.D = op.stokesDLmatrix(vesicle);
  if o.fmmDLP
    Nup = op.LPuprate * vesicle.N;
    Xup = [interpft(vesicle.X(1:end/2,:),Nup);...
      interpft(vesicle.X(end/2+1:end,:),Nup)];
    vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0);
    o.DLPnoCorr = op.stokesDLmatrixNoCorr(vesicleUp);
  else
    o.DLPnoCorr = [];
  end
else
  o.D = [];
  o.DLPnoCorr = [];
end

if o.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(vesicle.X(1:end/2,:),Nup);...
    interpft(vesicle.X(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0); 
  o.SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);
else
  o.SLPnoCorr = [];    
end

% Only form near-singular integration structure 

% Structures for deciding who
% is close, how close it is, who is closest, etc., needed in nearSingInt
[o.NearV2V,o.NearV2Wint] = vesicle.getZone(wallsInt,3);
      
if isempty(o.NearWint2Wint)
[o.NearWint2Wint,o.NearWint2V] = wallsInt.getZone(vesicle,3);
else
% there is no need to compute W2W again, since they do not move
[~,o.NearWint2V] = wallsInt.getZone(vesicle,2);
end

[~,o.NearV2Wext] = vesicle.getZone(wallsExt,2);
[~,o.NearWext2V] = wallsExt.getZone(vesicle,2);

if isempty(o.NearWint2Wext)
% there is no need to compute W2W again, since they do not move
[~,o.NearWint2Wext] = wallsInt.getZone(wallsExt,2);
end

if isempty(o.NearWext2Wint)
% there is no need to compute W2W again, since they do not move
[~,o.NearWext2Wint] = wallsExt.getZone(wallsInt,2);    
end
      

% Parts of rhs from previous solution.  The right-hand-side depends on
% whether we are doing an SDC correction or forming the provisional
% solution.

rhs1 = Xo;
rhs2 = zeros(N,nv);
rhs3Ext = wallsExt.u;
rhs3Int = wallsInt.u;

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP
% vesicle-vesicle and vesicle-wall interactions are handled
% implicitly in TimeMatVec
% END TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VISCOSITY CONTRAST
if any(vesicle.viscCont ~= 1)
  
  % Need to add jump to double-layer potential if using near-singular
  % integration so that we can compute boundary values for the
  % near-singular integration algorithm
  jump = 1/2*(1-vesicle.viscCont);
  DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);

  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
    kernelDirect = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLnewfmm;
    kernelDirect = @op.exactStokesDL;
  end
  
  
  density = Xo;
  % Use near-singular integration to compute double-layer
  % potential from previous solution 
  Fdlp = op.nearSingInt(vesicle,density,DLP,o.DLPnoCorr,...
       o.NearV2V,kernel,kernelDirect,vesicle,true,false);
  FDLPwallInt = op.nearSingInt(vesicle,density,DLP,[],o.NearV2Wint,...
        kernel,kernelDirect,wallsInt,false,false);
  FDLPwallExt = op.nearSingInt(vesicle,density,DLP,[],o.NearV2Wext,...
        kernel,kernelDirect,wallsExt,false,false);  
else
  % If no viscosity contrast, there is no velocity induced due to a
  % viscosity contrast  
  Fdlp = zeros(2*N,nv);
  FDLPwallInt = zeros(2*NbdInt,nvbdInt);
  FDLPwallExt = zeros(2*NbdExt,nvbdExt);
end

% add in viscosity contrast term due to each vesicle independent of the
% others (o.D * Xo) from the previous solution followed by the term due
% to all other vesicles (Fdlp)
if (any(vesicle.viscCont ~= 1) && ~o.SDCcorrect)
  DXo = op.exactStokesDLdiag(vesicle,o.D,Xo);
  rhs1 = rhs1 - (Fdlp + DXo) * diag(1./alpha);
end

% compute the double-layer potential due to all other vesicles from the
% appropriate linear combination of previous time steps.  Depends on
% time stepping order and vesicle-vesicle discretization
rhs3Ext = rhs3Ext + FDLPwallExt/o.dt;
rhs3Int = rhs3Int + FDLPwallInt/o.dt;


% START COMPUTING SINGLE-LAYER POTENTIAL FOR REPULSION
if o.repulsion 
  % Repulsion is handled explicitly between vesicles, vesicles-walls.  
  if ~o.fmm
      kernel = @op.exactStokesSL;
      kernelDirect = @op.exactStokesSL;
  else
      kernel = @op.exactStokesSLfmm;
      kernelDirect = @op.exactStokesSL;
  end
 
  Xrep = Xo;
  repulsion = vesicle.repulsionSchemeSimple(Xrep,o.repStrength,o.minDist,...
        [],wallsInt,wallsExt);

  Frepulsion = op.exactStokesSLdiag(vesicle,o.Galpert,repulsion);
  % diagonal term of repulsion


  SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
    
  % Use near-singular integration to compute single-layer potential
  % due to all other vesicles.  Need to pass function
  % op.exactStokesSL so that the layer potential can be computed at
  % far points and Lagrange interpolation points
  Frepulsion = Frepulsion + ...
    op.nearSingInt(vesicle,repulsion,SLP,o.SLPnoCorr,...
        o.NearV2V,kernel,kernelDirect,vesicle,true,false);
    
  % Evaluate the velocity on the walls due to the vesicles
  FREPwallInt = op.nearSingInt(vesicle,repulsion,SLP,[],...
        o.NearV2Wint,kernel,kernelDirect,wallsInt,false,false);    
  FREPwallExt = op.nearSingInt(vesicle,repulsion,SLP,[],...
        o.NearV2Wext,kernel,kernelDirect,wallsExt,false,false);    
  
  % keep the number of timesteps when repulsion is nonzero
  if norm(Frepulsion) ~= 0 
    global repuls
    repuls = repuls + 1;
  end
    
  rhs1 = rhs1 + o.dt*Frepulsion*diag(1./alpha);
  rhs3Ext = rhs3Ext - FREPwallExt;
  rhs3Int = rhs3Int - FREPwallInt;
end
% END COMPUTING SINGLE-LAYER POTENTIALS FOR REPULSION 



% START TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS
% This is done implicitly in TimeMatVec
% END TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS


% START TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION
% Vesicle-vesicle interactions and the presence of a
% viscosity contrast is irrelavent  
% rhs2 is the right-hand side for the inextensibility condition
rhs2 = rhs2 + vesicle.surfaceDiv(Xo); 
% END TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION

% The next makes sure that the rhs is all order one rather than have rhs3
% being order 1/o.dt and other two parts (rhs1 and rhs2) being order 1.
% This of course needs to be compensated in the TimeMatVec routine
if (any(vesicle.viscCont ~= 1) && ...
      strcmp(o.vesves,'implicit') && o.confined)
  rhs3Int = rhs3Int * o.dt;
  rhs3Ext = rhs3Ext * o.dt;
end

% Stack the right-hand sides in an alternating with respect to the
% vesicle fashion
rhs = [rhs1; rhs2];
rhs = rhs(:);
rhs = [rhs; rhs3Ext(:); rhs3Int(:)];
% Add on the no-slip boundary conditions on the solid walls
% Rotlet and Stokeslet equations
rhs = [rhs; zeros(3*(nvbd-1),1)];

% Use a preconditioner (block-diagonal preconditioner is implemented)
usePreco = o.usePreco;
useSpecPreco = o.useSpecPreco;
    
% START BUILDING BLOCK-DIAGONAL PRECONDITIONER
if usePreco && ~useSpecPreco
  
  % Build differential operators. 
  % Compute bending, tension, and surface divergence of current
  % vesicle configuration
  [Ben,Ten,Div] = vesicle.computeDerivs;
  bdiagVes.L = zeros(3*N,3*N,nv);
  bdiagVes.U = zeros(3*N,3*N,nv);

  % Build block-diagonal preconditioner of self-vesicle 
  % intearctions in matrix form
  for k=1:nv
    if any(vesicle.viscCont ~= 1)
      [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
        [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
            o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
        -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
        o.beta*Div(:,:,k) zeros(N)]);
    else
      [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
        [o.beta*eye(2*N) + ...
            o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
        -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
        o.beta*Div(:,:,k) zeros(N)]);
    end
  end
  o.bdiagVes = bdiagVes;
end % usePreco

% SOLVING THE SYSTEM USING GMRES
warning off
% any warning is printed to the terminal and the log file so
% don't need the native matlab version
initGMRES = [Xm;sigmaM];
initGMRES = initGMRES(:);
RS = RSm(:,2:end);
initGMRES = [initGMRES;etaMext(:);etaMint(:);RS(:)];

% Use GMRES to solve for new positions, tension, density
% function defined on the solid walls, and rotlets/stokeslets
if usePreco && ~useSpecPreco
  [Xn,iflag,~,I,~] = gmres(@(X) o.TimeMatVec(X,vesicle,[],...
      wallsInt,wallsExt),rhs,[],o.gmresTol,o.gmresMaxIter,...
      @o.preconditionerBD,[],initGMRES);
  iter = I(2);    

elseif useSpecPreco
  
  [Xn,iflag,~,I,~] = gmres(@(X) o.TimeMatVec(X,vesicle,[],...
      wallsInt,wallsExt),rhs,[],o.gmresTol,o.gmresMaxIter,...
      @(z) o.preconditionerSpectral(vesicle,[],[],z),[],initGMRES);
  iter = I(2);      
else
  [Xn,iflag,~,I,~] = gmres(@(X) o.TimeMatVec(X,vesicle,[],...
      wallsInt,wallsExt),rhs,[],o.gmresTol,o.gmresMaxIter);
  iter = I(2);
end
warning on
% END OF SOLVING THE SYSTEM USING GMRES

% allocate space for positions, tension, and density function
X = zeros(2*N,nv);
sigma = zeros(N,nv);
etaInt = zeros(2*NbdInt,nvbdInt);
RS = zeros(3,nvbd);

% unstack the positions and tensions
for k=1:nv
  X(:,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigma(:,k) = Xn((3*k-1)*N+1:3*k*N);
end

% unstack the density function
Xn = Xn(3*nv*N+1:end);
etaExt = Xn(1:2*NbdExt);  % assuming nvbdExt = 1
for k = 1:nvbdInt
  etaInt(:,k) = Xn(2*NbdExt+(k-1)*2*NbdInt+1:2*NbdExt+2*k*NbdInt);
end  


% unstack the rotlets and stokeslets
otlets = Xn(2*NbdExt+2*nvbdInt*NbdInt+1:end);  
for k = 2:nvbd
  istart = (k-2)*3+1;
  iend = 3*(k-1);
  RS(:,k) = otlets(istart:iend);
end

end % timeStepSimpleDiffDisc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zOut = solveIminusD(o,zIn,vesicle)
% zOut = o.solveIminusD(zIn,vesicle) inverts the linear system \alpha*I
% - DLP where alpha=0.5*(1+viscCont) and DLP is the double-layer
% potential.  This is the opeartor that needs to be inverted to find the
% velocity when there is a viscosity contrast involved.  If there is no
% viscosity contrast, alpha = 1 and DLP = 0 so that zOut = zIn


% solve with block-diagonal preconditioned GMRES.  Integral equation is
% of the form identity + compact.  Need a bit more accuracy in this
% gmres solve as it is an inner iteration within an outer GMRES
% iteration.  This hides the fact that the GMRES solver is not linear
warning off
[zIn,flag,relres,iter] = gmres(@(X) o.IminusD(X,vesicle),zIn(:),...
    [],1e-2*o.gmresTol,min(2*vesicle.N*vesicle.nv,o.gmresMaxIter),...
    @(z) o.precoIminusD(z,vesicle));
warning on

% Sort output in appropriate format
zOut = zeros(2*vesicle.N,vesicle.nv);
for k = 1:vesicle.nv
  zOut(:,k) = zIn((k-1)*2*vesicle.N+1:2*k*vesicle.N);
end


end % solveIminusD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = IminusD(o,X,vesicle)
% val = IminusD(X,vesicle) does a matrix-vector multiply with the matrix
% (alpha*I - D) where alpha = 0.5*(1+viscCont) and D is the double-layer
% potential.  This matrix needs to be inverted when solving for the
% velocity field of a given vesicle configuration and tension

op = o.op;
N = vesicle.N;
nv = vesicle.nv;
alpha = 0.5*(1+vesicle.viscCont);

Xm = zeros(2*N,nv);
for k = 1:nv
  Xm(:,k) = X((k-1)*2*N+1:k*2*N);
end

% "jump" term since we are computing alpha * I - DLP
val = Xm*diag(alpha);

% self-interaction term
val = val - op.exactStokesDLdiag(vesicle,o.D,Xm);

% compute term for evaluating limiting value along boundary of domain
jump = 0.5*(1-vesicle.viscCont);
DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);

% kernel for evaulating the double-layer potential
if ~o.fmmDLP
  kernel = @op.exactStokesDL;
  kernelDirect = @op.exactStokesDL;
else
  kernel = @op.exactStokesDLnewfmm;
  kernelDirect = @op.exactStokesDL;
end

% potential due to other vesicles
Fdlp = op.nearSingInt(vesicle,Xm,DLP,o.DLPnoCorr,...
  o.NearV2V,kernel,kernelDirect,vesicle,true,false);

% add potential due to self and other vesicles.  Jump is already
% included
val = val - Fdlp;
val = val(:);

end % IminusD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = sigDenMatVec(o,sigma,vesicle,walls,wallsInt,wallsExt)
% val = sigDenMatVec(o,sigma,vesicle,walls,wallsInt,wallsExt) does the
% matvec multiply required to find the tension and density function of a
% given vesicle configuration and farfield or solid wall velocity field
N = vesicle.N; % Number of points
nv = vesicle.nv; % Number of vesicles
if o.confined
  if ~o.diffDiscWalls  
    Nbd = walls.N;
    nvbd = walls.nv;
    NbdInt = 0;
    NbdExt = 0;
    nvbdInt = 0;
    nvbdExt = 0;
  else
    NbdInt = wallsInt.N;
    NbdExt = wallsExt.N;
    nvbdInt = wallsInt.nv;
    nvbdExt = wallsExt.nv;
    nvbd = nvbdInt + nvbdExt;
    Nbd = 0;
  end
else
  Nbd = 0;
  nvbd = 0;
  NbdInt = 0;
  NbdExt = 0;
  nvbdInt = 0;
  nvbdExt = 0;
end

% need to do some manipulations to the TimeMatVec matrix applied only
% to the first component, which is the velocity
u = zeros(2*N,nv);

for k = 1:nv
  istart = (k-1)*3*N + 1;
  iend = istart + 2*N - 1;
  u(:,k) = sigma(istart:iend);
end

% want interactions to be implicit so that the most accurate tension and
% density functions are found

inextens = o.solver;
o.solver = 'method1';
dt = o.dt;
o.dt = 1;

% Do a matvec but let the incoming postitions be zero since they are
% handled by the initial condition
% Can overwrite sigma as we don't need it from this point onwards
sigma = o.TimeMatVec(sigma,vesicle,walls,wallsInt,wallsExt);

% change back to old vesicle-vesicle and vesicle-boundary interactions
o.solver = inextens;
o.dt = dt;

% part that corresponds to the velocity and tension
valVel = zeros(2*N,nv);
valTen = zeros(N,nv);
valDen = zeros(2*Nbd,nvbd);
valDenInt = zeros(2*NbdInt,nvbdInt);
valRS = zeros(3,nvbd-1);
for k = 1:nv
  valVel(:,k) = sigma((k-1)*3*N+1:(3*k-1)*N);
  valTen(:,k) = sigma((3*k-1)*N+1:3*k*N);
end

% part that corresponds to the density function
if ~o.diffDiscWalls
  for k = 1:nvbd
    valDen(:,k) = sigma(3*nv*N+(k-1)*2*Nbd+1:3*nv*N+k*2*Nbd);
  end
else
  valDenExt = sigma(3*nv*N+1:3*nv*N+2*NbdExt);
  for k = 1:nvbdInt
    valDenInt(:,k) = sigma(3*nv*N+2*NbdExt+(k-1)*2*NbdInt+1:...
      3*nv*N+2*NbdExt+k*2*NbdInt);      
  end
end

% part that corresponds to the Stokeslets and Rotlets
if ~o.diffDiscWalls
  for k = 2:nvbd
    istart = 3*nv*N + 2*nvbd*Nbd + 3*(k-2) + 1;
    iend = istart + 2;  
    valRS(:,k-1) = sigma(istart:iend);
  end
else
  for k = 2:nvbd
    istart = 3*nv*N + 2*NbdExt + 2*nvbdInt*NbdInt + 3*(k-2) + 1;
    iend = istart + 2;  
    valRS(:,k-1) = sigma(istart:iend);   
  end
end

% kernel for single-layer potential.  Only difference is if the FMM is
% used or not
op = o.op;
if ~o.fmm
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
  kernelDirect = @op.exactStokesSL;
end

% bending due to the velocity
f = vesicle.tracJump(u,zeros(N,nv));


SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
Fslp = op.nearSingInt(vesicle,f,SLP,o.SLPnoCorr,...
    o.NearV2V,kernel,kernelDirect,vesicle,true,false);
      
% Evaluate single-layer potential due to all vesicles on
% the solid walls WITH near-singular integration
if o.confined
  if ~o.diffDiscWalls
    FSLPwall = op.nearSingInt(vesicle,f,SLP,[],...
      o.NearV2W,kernel,kernelDirect,walls,false,false);
  else
    FSLPwallInt = op.nearSingInt(vesicle,f,SLP,[],...
        o.NearV2Wint,kernel,kernelDirect,wallsInt,false,false);
    FSLPwallExt = op.nearSingInt(vesicle,f,SLP,[],...
        o.NearV2Wext,kernel,kernelDirect,wallsExt,false,false);
  end
else
  FSLPwall = [];
  FSLPwallInt = [];
  FSLPwallExt = [];
end

% multiply top row of matrix by alpha
alpha = (1+vesicle.viscCont)/2; 
valVel = valVel * diag(alpha);

valVel = valVel + op.exactStokesSLdiag(vesicle,o.Galpert,f) + Fslp;

val = zeros(3*N*nv,1);
for k=1:nv
  val((k-1)*3*N+1:3*k*N) = [valVel(:,k);valTen(:,k)];
end

% subtract off terms that TimeMatVec introduces but we do not have in
% this linear system
% also,
% stack the different components coming from the inextensibility, solid
% walls, and rotlets/stokeslets
if o.confined
  if ~o.diffDiscWalls
    valDen = valDen - FSLPwall;
    val = [val(:);valDen(:);valRS(:)];
  else
    valDenInt = valDenInt - FSLPwallInt;
    valDenExt = valDenExt - FSLPwallExt;
    val = [val(:);valDenExt(:);valDenInt(:);valRS(:)];
  end
end

end % sigDenMatVec


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = sigDenMatVecLowAcc(o,sigma,vesicle)
% val = sigDenMatVecLowAcc(o,sigma,vesicle) does the
% matvec multiply required to find the tension and density function of a
% given vesicle configuration and farfield. It does not use a special
% quadrature for self interactions and near singular integration
N = vesicle.N; % Number of points
nv = vesicle.nv; % Number of vesicles

op = o.op;

% need to do some manipulations to the TimeMatVec matrix applied only
% to the first component, which is the velocity
u = zeros(2*N,nv);

for k = 1:nv
  istart = (k-1)*3*N + 1;
  iend = istart + 2*N - 1;
  u(:,k) = sigma(istart:iend);
end

% want interactions to be implicit so that the most accurate tension and
% density functions are found

inextens = o.solver;
o.solver = 'method1';
dt = o.dt;
o.dt = 1;

% Do a matvec but let the incoming postitions be zero since they are
% handled by the initial condition
% Can overwrite sigma as we don't need it from this point onwards
sigma = o.TimeMatVecLowAcc(sigma,vesicle);

% change back to old vesicle-vesicle and vesicle-boundary interactions
o.solver = inextens;
o.dt = dt;

% part that corresponds to the velocity and tension
valVel = zeros(2*N,nv);
valTen = zeros(N,nv);
for k = 1:nv
  valVel(:,k) = sigma((k-1)*3*N+1:(3*k-1)*N);
  valTen(:,k) = sigma((3*k-1)*N+1:3*k*N);
end

% bending due to the velocity
f = vesicle.tracJump(u,zeros(N,nv));
if ~o.fmm
  Fslp = op.exactStokesSL(vesicle,f,o.SLPnoCorr);  
else
  Fslp = op.exactStokesSLfmm(vesicle,f,o.SLPnoCorr);  
end
     
% subtract off terms that TimeMatVec introduces but we do not have in
% this linear system

% multiply top row of matrix by alpha
alpha = (1+vesicle.viscCont)/2; 
valVel = valVel * diag(alpha);

valVel = valVel + op.exactStokesSLdiag(vesicle,o.SLPnoCorr,f) + Fslp;

val = zeros(3*N*nv,1);
for k=1:nv
  val((k-1)*3*N+1:3*k*N) = [valVel(:,k);valTen(:,k)];
end

end % sigDenMatVecLowAcc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = letsIntegrals(o,otlets,etaM,etaMint,walls,wallsInt)
% z = letsIntegrals(o,otlets,etaM,etaMint,walls,wallsInt) integrates 
% the density function to enforce constraints on stokeslets and rotlets
% if there are walls with different discretizations (i.e. one large outer
% wall and multiple small inner walls), then takes only the inner walls)

if ~o.diffDiscWalls
  Nbd = walls.N;
  nvbd = walls.nv;
else
  NbdInt = wallsInt.N;
  nvbd = wallsInt.nv + 1; % assuming one outer wall
end

z = zeros(3*(nvbd-1),1);

if ~o.diffDiscWalls
  for k = 2:nvbd
    % two stokeslet terms per inner boundary  
    stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
    % one rotlet term per inner boundary
    rotlet = otlets(3*(k-1));
    % integral of density function dotted with [1;0]
    % is one stokeslet
    ind = 3*(k-2)+1;
    z(ind) = -2*pi*stokeslet(1) + ...
      sum(etaM(1:Nbd,k).*walls.sa(:,k))*2*pi/Nbd;
    % integral of density fuction dotted with [0;1]
    % is the other stokeslet
    z(ind+1) = -2*pi*stokeslet(2) + ...
      sum(etaM(Nbd+1:2*Nbd,k).*walls.sa(:,k))*2*pi/Nbd;
    % integral of density function dotted with (-y,x)
    % is the rotlet
    z(ind+2) = -2*pi*rotlet + sum(...
      ((walls.X(Nbd+1:2*Nbd,k)).*etaM(1:Nbd,k) - ...
      (walls.X(1:Nbd,k)).*etaM(Nbd+1:2*Nbd,k)).*...
      walls.sa(:,k))*2*pi/Nbd;
  end % k
else
  for k = 2:nvbd
    % two stokeslet terms per inner boundary  
    stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
    % one rotlet term per inner boundary
    rotlet = otlets(3*(k-1));
    % integral of density function dotted with [1;0]
    % is one stokeslet
    ind = 3*(k-2)+1;
    z(ind) = -2*pi*stokeslet(1) + ...
      sum(etaMint(1:NbdInt,k-1).*wallsInt.sa(:,k-1))*2*pi/NbdInt;
    % integral of density fuction dotted with [0;1]
    % is the other stokeslet
    z(ind+1) = -2*pi*stokeslet(2) + ...
      sum(etaMint(NbdInt+1:2*NbdInt,k-1).*wallsInt.sa(:,k-1))*2*pi/NbdInt;
    % integral of density function dotted with (-y,x)
    % is the rotlet
    z(ind+2) = -2*pi*rotlet + sum(...
      ((wallsInt.X(NbdInt+1:2*NbdInt,k-1)).*etaMint(1:NbdInt,k-1) - ...
      (wallsInt.X(1:NbdInt,k-1)).*etaMint(NbdInt+1:2*NbdInt,k-1)).*...
      wallsInt.sa(:,k-1))*2*pi/NbdInt;
  end % k
end % diffDiscWalls

end % letsIntegrals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = RSlets(o,X,center,stokeslet,rotlet)
% vel = RSlets(o,X,center,stokeslet,rotlet) evaluates the velocity due
% to the stokeslet and rotlet terms.  Center of the rotlet and
% stokeslet is contained in center

oc = curve;
% set of points where we are evaluating the velocity
[x,y] = oc.getXY(X);
% the center of the rotlet/stokeslet terms
[cx,cy] = oc.getXY(center);

% distance squared
rho2 = (x-cx).^2 + (y-cy).^2;

% x component of velocity due to the stokeslet and rotlet
LogTerm = -0.5*log(rho2)*stokeslet(1);
rorTerm = 1./rho2.*((x-cx).*(x-cx)*stokeslet(1) + ...
    (x-cx).*(y-cy)*stokeslet(2));
RotTerm = (y-cy)./rho2*rotlet;
velx = 1/4/pi*(LogTerm + rorTerm) + RotTerm;

% y component of velocity due to the stokeslet and rotlet
LogTerm = -0.5*log(rho2)*stokeslet(2);
rorTerm = 1./rho2.*((y-cy).*(x-cx)*stokeslet(1) + ...
    (y-cy).*(y-cy)*stokeslet(2));
RotTerm = -(x-cx)./rho2*rotlet;
vely = 1/4/pi*(LogTerm + rorTerm) + RotTerm;

% velocity
vel = [velx;vely];


end % RSlets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF DIFFERENT PRECONDITIONERS INCLUDING BLOCK-DIAGONAL, ONE FOR
% THE SYSTEM THAT SOLVES FOR THE TENSION AND DENSITY GIVEN A POSITION,
% MULTIGRID IDEAS, SCHUR COMPLEMENTS, AND ANALYTIC (BASED ON A
% CIRCLE).  THE BLOCK-DIAGONAL PRECONDITIONER IS THE MOST ROBUST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerBD(o,z)
% val = preconditionBD(z) applies the block diagonal preconditioner
% required by preconditioned-GMRES to the vector z

if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
  nv = size(o.bdiagVes.L,3); % number of vesicles
  N = size(o.bdiagVes.L,1)/3; % number of points
elseif strcmp(o.solver,'method3')
  nv = size(o.bdiagVes.schur,3); % number of vesicles
  N = size(o.bdiagVes.schur,1)/2; % number of vesicles
elseif strcmp(o.solver,'method4')
  nv = size(o.bdiagVes.schur,3); % number of vesicles
  N = size(o.bdiagVes.schur,1); % number of vesicles
end

% extract the position and tension part.  Solid walls is
% handled in the next section of this routine
zves = z(1:3*N*nv);

valVes = zeros(3*N*nv,1);
% precondition with the block diagonal preconditioner for the
  % vesicle position and tension
if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
  for k=1:nv
    valVes((k-1)*3*N+1:3*k*N) = o.bdiagVes.U(:,:,k)\...
      (o.bdiagVes.L(:,:,k)\zves((k-1)*3*N+1:3*k*N));
  end % k
  
elseif strcmp(o.solver,'method3')
  for k = 1:nv
    % seperate position and tension terms  
    istart = (k-1)*3*N+1;
    iend = istart + 2*N - 1;
    bx = zves(istart:iend); 
    istart = iend + 1;
    iend = istart + N - 1;
    bsig = zves(istart:iend);
    
    % use schur decomposition to operate preconditioner
    rhs = bx + o.dt*o.bdiagVes.GT(:,:,k)*o.bdiagVes.DGT(:,:,k)*bsig;
    pos = o.bdiagVes.schur(:,:,k)*rhs;
    rhs = bsig + o.bdiagVes.DGB(:,:,k)*pos;
    sig = o.bdiagVes.DGT(:,:,k)*rhs;
    valVes((k-1)*3*N+1:3*k*N) = [pos;sig];
    
  end % k
elseif strcmp(o.solver,'method4')
  for k = 1:nv
    % seperate position and tension terms  
    istart = (k-1)*3*N+1;
    iend = istart + 2*N - 1;
    bx = zves(istart:iend); 
    istart = iend + 1;
    iend = istart + N - 1;
    bsig = zves(istart:iend);
    
    % use schur decomposition to operate preconditioner
    rhs = bsig + o.bdiagVes.DGB(:,:,k)*o.bdiagVes.IpBen(:,:,k)*bx;
    sig = o.bdiagVes.schur(:,:,k)*rhs;
    rhs = bx + o.dt*o.bdiagVes.GT(:,:,k)*sig;
    pos = o.bdiagVes.IpBen(:,:,k)*rhs;
    valVes((k-1)*3*N+1:3*k*N) = [pos;sig];
    
  end % k

end % o.solver

% part of z from solid walls
zwalls = z(3*N*nv+1:end);
if ~o.fastDirect
  % exact factorization of (I+wallsDLP) 
  % this matrix is well-coditioned since it is in the form I + DLP
  if o.outOfCore
    % then we use memory map to reach entries of the matrix  
    memmapInv = o.bdiagWall;  
    
    % Block matrix-vector multiplication b/c we cannot extract this large
    % matrix on memory
    numRows = numel(memmapInv.Data.M(:,1));
    valWalls = zeros(numRows,1);
      
    bsize = min(numRows,o.maxbsize);
    for i = 1:floor(numRows/bsize)
      istart = 1+(i-1)*bsize;
      iend = i*bsize;
      row = memmapInv.Data.M(istart:iend,:);
      valWalls(istart:iend,1) = row*zwalls;
    end
    istart = 1+i*bsize;
    row = memmapInv.Data.M(istart:end,:);
    valWalls(istart:end,1) = row*zwalls;
    clear row;     
      
  else
    % if not out-of-core, then the inverse is on the memory      
    valWalls = o.bdiagWall * zwalls;
  end
  
else
  % inexact factorization of (I + wallsDLP) using HODLR  
  valWalls = o.bdiagWallFD(zwalls);  
end

% stack the two componenets of the preconditioner
val = [valVes;valWalls];


end % preconditionerBD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mat = wallsPrecond(o,walls)
% wallsPrecond(walls) computes the matrix which is the 
% exact inverse of
% the double-layer potential for stokes flow in a bounded domain.  Used
% in the preconditioner for vesicle simulations and capsules.m/computeEta
% which computes eta and RS when there is no vesicle.


Nbd = walls.N;
nvbd = walls.nv;
oc = curve;
[x,y] = oc.getXY(walls.X);
[nory,norx] = oc.getXY(walls.xt);
nory = -nory;
sa = walls.sa;
[cx,cy] = oc.getXY(walls.center);

% Allocate space for blocks of matrix that carries the double- layer
% potential, rotlets, and stokeslets to the velocity and the conditions
% in (A4) and (A5) in Rahimian et al.
M11 = zeros(2*Nbd*nvbd,2*Nbd*nvbd);
M12 = zeros(2*Nbd*nvbd,3*(nvbd-1));
M21 = zeros(3*(nvbd-1),2*Nbd*nvbd);


% Self interaction terms with the jump coming from the double layer
% potential
M11(1:2*Nbd,1:2*Nbd) = M11(1:2*Nbd,1:2*Nbd) + o.wallN0(:,:,1);
jump = - 1/2; 
for k = 1:nvbd
istart = (k-1)*2*Nbd+1;
iend = 2*k*Nbd;
M11(istart:iend,istart:iend) = M11(istart:iend,istart:iend) + ...
    jump*eye(2*Nbd) + o.wallDLP(:,:,k);
end


for ktar = 1:nvbd % loop over targets
itar = 2*(ktar-1)*Nbd + 1;
jtar = 2*ktar*Nbd;
K = [(1:ktar-1) (ktar+1:nvbd)];

D = zeros(2*Nbd,2*Nbd);
for ksou = K % loop over all other walls
  isou = 2*(ksou-1)*Nbd + 1;
  jsou = 2*ksou*Nbd;

  xtar = x(:,ktar); ytar = y(:,ktar);
  xtar = xtar(:,ones(Nbd,1)); 
  ytar = ytar(:,ones(Nbd,1));

  xsou = x(:,ksou); ysou = y(:,ksou);
  xsou = xsou(:,ones(Nbd,1))';
  ysou = ysou(:,ones(Nbd,1))';

  norxtmp = norx(:,ksou); norytmp = nory(:,ksou);
  norxtmp = norxtmp(:,ones(Nbd,1))';
  norytmp = norytmp(:,ones(Nbd,1))';

  satmp = sa(:,ksou);
  satmp = satmp(:,ones(Nbd,1))';

  rho2 = (xtar-xsou).^2 + (ytar-ysou).^2;

  coeff = 1/pi*((xtar-xsou).*norxtmp + ...
      (ytar-ysou).*norytmp).*satmp./rho2.^2;

  D(1:Nbd,:) = 2*pi/Nbd*[coeff.*(xtar-xsou).^2 ...
      coeff.*(xtar-xsou).*(ytar-ysou)];
  D(Nbd+1:end,:) = 2*pi/Nbd*[coeff.*(ytar-ysou).*(xtar-xsou) ...
      coeff.*(ytar-ysou).^2];

  M11(itar:jtar,isou:jsou) = D;
end %end ktar
end %end ksou

% These compute the integral of the density function around each of the
% inner componenents of the geometry
for k = 1:nvbd-1
icol = 3*(k-1)+1;
istart = 2*k*Nbd+1;
iend = istart + Nbd - 1;
M21(icol,istart:iend) = 2*pi/Nbd*sa(:,k+1)';
M21(icol+2,istart:iend) = 2*pi/Nbd*sa(:,k+1)'.*y(:,k+1)';
istart = istart + Nbd;
iend = iend + Nbd;
M21(icol+1,istart:iend) = 2*pi/Nbd*sa(:,k+1)';
M21(icol+2,istart:iend) = -2*pi/Nbd*sa(:,k+1)'.*x(:,k+1)';
end % k


% This is the evaluation of the velocity field due to the stokeslet
% and rotlet terms
for k = 1:nvbd - 1
for ktar = 1:nvbd
  rho2 = (x(:,ktar) - cx(k+1)).^2 + (y(:,ktar) - cy(k+1)).^2;
  istart = (ktar-1)*2*Nbd + 1;
  iend = istart + Nbd - 1;

  icol = 3*(k-1)+1;
  M12(istart:iend,icol) = ...
    M12(istart:iend,icol) + ...
    1/4/pi*(-0.5*log(rho2) + (x(:,ktar)-cx(k+1))./rho2.*...
        (x(:,ktar)-cx(k+1)));
  M12(istart + Nbd:iend + Nbd,icol) = ...
    M12(istart + Nbd:iend + Nbd,icol) + ...
    1/4/pi*((x(:,ktar)-cx(k+1))./rho2.*(y(:,ktar)-cy(k+1)));

  icol = 3*(k-1)+2;
  M12(istart:iend,icol) = ...
    M12(istart:iend,icol) + ...
    1/4/pi*((y(:,ktar)-cy(k+1))./rho2.*(x(:,ktar)-cx(k+1)));
  M12(istart + Nbd:iend + Nbd,icol) = ...
    M12(istart + Nbd:iend + Nbd,icol) + ...
    1/4/pi*(-0.5*log(rho2) + (y(:,ktar)-cy(k+1))./rho2.*...
        (y(:,ktar)-cy(k+1)));

  icol = 3*(k-1)+3;
  M12(istart:iend,icol) = ...
    M12(istart:iend,icol) + ...
    (y(:,ktar)-cy(k+1))./rho2;
  M12(istart + Nbd:iend + Nbd,icol) = ...
    M12(istart + Nbd:iend + Nbd,icol) - ...
    (x(:,ktar)-cx(k+1))./rho2;
end
end

% different combinations of the density functions have to multiply to
% 2*pi multiplied by rotlet or stokeslet terms
M22 = -2*pi*eye(3*(nvbd-1));

% Save the wall2wall interaction matrices if not matrix free
if ~o.matFreeWalls
o.wallDLPandRSmat = [M11 M12; M21 M22];
end

% invert the matrix
Mat = ([M11 M12; M21 M22])\eye(2*nvbd*Nbd + 3*(nvbd-1));

    
end % wallsPrecond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rowVectN,colVectN] = findNewIds(o,newIDs,oldIDs,...
        rowVect,colVect)

% First find the new ids corresponding to rowVect containing row ids
rowVectN = zeros(size(rowVect));
for k = 1 : numel(rowVect)
  rowVectN(k) = newIDs(oldIDs==rowVect(k));
end

% Now, find the new ids corresponding to colVect containing column ids
colVectN = zeros(size(colVect));
for k = 1 : numel(colVect)
  colVectN(k) = newIDs(oldIDs==colVect(k));
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wallsBDprecondConstruct(o,wallsInt,wallsExt)

oc = curve;

[xInt,yInt] = oc.getXY(wallsInt.X);

[noryInt,norxInt] = oc.getXY(wallsInt.xt);
noryInt = -noryInt;
saInt = wallsInt.sa;
[cxInt,cyInt] = oc.getXY(wallsInt.center);


% discretization info
Next = wallsExt.N;
Nint = wallsInt.N;
nv = wallsInt.nv;

% Block corresponding to the external wall
o.invM11Ext = (-1/2*eye(2*Next) + o.wallDLPext+o.wallN0(:,:,1))\eye(2*Next);

% Exactly factorize each of the pillars (but include RS, too)
o.invM11Int = zeros(2*Nint+3,2*Nint+3,wallsInt.nv);
for k = 1 : nv
  M11 = -1/2*eye(2*Nint)+o.wallDLPint(:,:,k);
  M21 = zeros(3,2*Nint); M12 = zeros(2*Nint,3); M22 = -2*pi*eye(3);
  
  M21(1,1:Nint) = 2*pi/Nint*saInt(:,k)';
  M21(3,1:Nint) = 2*pi/Nint*saInt(:,k)'.*yInt(:,k)';
  M21(2,Nint+1:2*Nint) = 2*pi/Nint*saInt(:,k)';
  M21(3,Nint+1:2*Nint) = -2*pi/Nint*saInt(:,k)'.*xInt(:,k)';

  rho2 = (xInt(:,k) - cxInt(k)).^2 + (yInt(:,k) - cyInt(k)).^2;
  M12(1:Nint,1) = 1/4/pi*(-0.5*log(rho2)+(xInt(:,k)-cxInt(k))./rho2.*...
        (xInt(:,k)-cxInt(k)));

  M12(1:Nint,2) = 1/4/pi*((yInt(:,k)-cyInt(k))./rho2.*(xInt(:,k)-cxInt(k)));

  M12(Nint+1:2*Nint,2) = 1/4/pi*(-0.5*log(rho2) + (yInt(:,k)-cyInt(k))./rho2.*...
        (yInt(:,k)-cyInt(k)));

  M12(1:Nint,3) = (yInt(:,k)-cyInt(k))./rho2;

  M12(Nint+1:2*Nint,3) = - (xInt(:,k)-cxInt(k))./rho2;

 o.invM11Int(:,:,k) = [M11 M12;M21 M22]\eye(2*Nint+3);
end


end % wallsBDprecondConstruct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = wallsBDprecondApply(o,z)
% Apply block-diagonal preconditioner for walls

% Get the size information
Next = size(o.invM11Ext,1)/2;
Nint = (size(o.invM11Int,1)-3)/2;
nvbdInt = size(o.invM11Int,3);
    
% separate z into DLP for external wall,posts and rotlets and stokeslets
zExt = z(1:2*Next);
zInt = z(2*Next+1:2*Next+2*Nint*nvbdInt);
zRS  = z(2*Next+2*Nint*nvbdInt+1:end);


% precondition for the exterior wall
valExt = o.invM11Ext*zExt;

% precondition for the posts
valInt = zeros(size(zInt));
valRS = zeros(size(zRS));
for k = 1 : nvbdInt
  istart = (k-1)*2*Nint+1;
  iend = istart+2*Nint-1;
  denInt = zInt(istart:iend);
  RSint = zRS(3*(k-1)+1:3*(k-1)+3);
  zAll = [denInt;RSint];
  valAll = o.invM11Int(:,:,k)*zAll;
  valInt(istart:iend) = valAll(1:2*Nint);
  valRS(3*(k-1)+1:3*(k-1)+3) = valAll(2*Nint+1:end);
end

% precondition for the RS

% put them together
val = [valExt;valInt;valRS];
    
end % wallsBDprecondConstruct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mat = wallsPrecondDiffDisc(o,wallsInt,wallsExt)
% wallsPrecondFD(walls,rhs) computes the inverse of the double-layer potential 
% + rotlets and stokeslets for stokes flow in a bounded domain 
% using Fast Direct Solver. 
% This replaces wallsPrecond when we have a complex geometry and we do not
% compute the exact inverse of the DLP + RS in those cases.

% discretization info
Next = wallsExt.N;
Nint = wallsInt.N;
nv = wallsInt.nv+1;

% if we are going to use HODLR, then we need to order the points using
% Z-curve
% That's why, we change the pointers and build the ordered matrix
if o.fastDirect
  fileNameM11 = [o.wallMatFile 'M11.bin']; % M11
  fileNameM12 = [o.wallMatFile 'M12.bin']; % 
  fileNameM21 = [o.wallMatFile 'M21.bin']; %
  fileNameM22 = [o.wallMatFile 'M22.bin'];  
  % Shuffling order
  ho = o.hodlr;
  idWall2gids = ho.pnts2qtree(wallsInt.X,wallsExt.X,0);


  % Get the old and new ids of M11, M12, M21, M22
  newIDs = [1:Next+Nint*(nv-1)]';
  pntIDs = idWall2gids(newIDs,1);
  wallIDs = idWall2gids(newIDs,2);

  extWallIDs = find(wallIDs==1);
  intWallIDs = find(wallIDs~=1);

  oldIDs = zeros(size(newIDs));
  oldIDs(extWallIDs) = pntIDs(extWallIDs);
  oldIDs(intWallIDs) = pntIDs(intWallIDs)+2*Nint*(wallIDs(intWallIDs)-2)+2*Next;

  oldIDsY = zeros(size(newIDs));
  oldIDsY(extWallIDs) = oldIDs(extWallIDs)+Next;
  oldIDsY(intWallIDs) = oldIDs(intWallIDs)+Nint;

  oldIDs = [oldIDs;oldIDsY];
  newIDs = [2*newIDs-1;2*newIDs];
end % if o.fastDirect

tForm0 = tic;
%if ~o.fastDirect
%message = ['Building and factorizing wall2wall interaction matrix exactly...'];
%o.om.writeMessage(message,'%s\n');
%else
%message = ['Building wall2wall interaction matrix...'];
%o.om.writeMessage(message,'%s\n');    
%end

oc = curve;

[xInt,yInt] = oc.getXY(wallsInt.X);
[xExt,yExt] = oc.getXY(wallsExt.X);

[noryInt,norxInt] = oc.getXY(wallsInt.xt);
noryInt = -noryInt;
saInt = wallsInt.sa;
[cxInt,cyInt] = oc.getXY(wallsInt.center);

[noryExt,norxExt] = oc.getXY(wallsExt.xt);
noryExt = -noryExt;
saExt = wallsExt.sa;

% Allocate space for blocks of matrix that carries the double- layer
% potential, rotlets, and stokeslets to the velocity and the conditions
% in (A4) and (A5) in Rahimian et al.
M11 = zeros(2*Nint*(nv-1)+2*Next);
M12 = zeros(2*Nint*(nv-1)+2*Next,3*(nv-1));
M21 = zeros(3*(nv-1),2*Nint*(nv-1)+2*Next);


% Self interaction terms with the jump coming from the double layer
% potential
if o.fastDirect
  [rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(1:2*Next)',(1:2*Next)');
else
  rowVectN = (1:2*Next)';
  colVectN = (1:2*Next)';
end
M11(rowVectN,colVectN) = M11(rowVectN,colVectN) + o.wallN0(:,:,1);

jump = - 1/2;
M11(rowVectN,colVectN) = M11(rowVectN,colVectN) + ...
jump*eye(2*Next) + o.wallDLPext;

for k = 1:nv-1
  istart = 2*Next + (k-1)*2*Nint+1;
  iend = istart + 2*Nint-1;

  if o.fastDirect
    [rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(istart:iend)',(istart:iend)');
  else
    rowVectN = (istart:iend)';
    colVectN = (istart:iend)';
  end
  M11(rowVectN,colVectN) = M11(rowVectN,colVectN) + ...
    jump*eye(2*Nint) + o.wallDLPint(:,:,k);
end


% -----------------------------------------------------------------------
% EXTERIOR WALL
% Do ktar = 1, ksou = 2 : nv and 
% ktar = 2 : nv, ksou = 1 separately and then the rest
% in the old loop

for ksou = 1:nv-1 % loop over all other walls
  isou = 2*Next + 2*(ksou-1)*Nint + 1;
  jsou = 2*Next + 2*ksou*Nint;

  xtar = xExt(:); ytar = yExt(:);
  xtar = xtar(:,ones(Nint,1)); 
  ytar = ytar(:,ones(Nint,1));

  xsou = xInt(:,ksou); ysou = yInt(:,ksou);
  xsou = xsou(:,ones(Next,1))';
  ysou = ysou(:,ones(Next,1))';

  norxtmp = norxInt(:,ksou); norytmp = noryInt(:,ksou);
  norxtmp = norxtmp(:,ones(Next,1))';
  norytmp = norytmp(:,ones(Next,1))';

  satmp = saInt(:,ksou);
  satmp = satmp(:,ones(Next,1))';

  rho2 = (xtar-xsou).^2 + (ytar-ysou).^2;

  coeff = 1/pi*((xtar-xsou).*norxtmp + ...
    (ytar-ysou).*norytmp).*satmp./rho2.^2;

  if o.fastDirect
    [rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(1:Next)',(isou:jsou)');
  else
    rowVectN = (1:Next)';
    colVectN = (isou:jsou)';
  end
  M11(rowVectN,colVectN) = 2*pi/Nint*[coeff.*(xtar-xsou).^2 ...
    coeff.*(xtar-xsou).*(ytar-ysou)];

  if o.fastDirect
    [rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(Next+1:2*Next)',(isou:jsou)');
  else
    rowVectN = (Next+1:2*Next)';
    colVectN = (isou:jsou)';
  end
  M11(rowVectN,colVectN) = 2*pi/Nint*[coeff.*(ytar-ysou).*(xtar-xsou) ...
    coeff.*(ytar-ysou).^2];

end %end ksou

for ktar = 1:nv-1 % loop over targets
  itar = 2*Next + 2*(ktar-1)*Nint + 1;
  jtar = 2*Next + 2*ktar*Nint;

  D = zeros(2*Nint,2*Next);

  xtar = xInt(:,ktar); ytar = yInt(:,ktar);
  xtar = xtar(:,ones(Next,1)); 
  ytar = ytar(:,ones(Next,1));

  xsou = xExt(:); ysou = yExt(:);
  xsou = xsou(:,ones(Nint,1))';
  ysou = ysou(:,ones(Nint,1))';

  norxtmp = norxExt(:); norytmp = noryExt(:);
  norxtmp = norxtmp(:,ones(Nint,1))';
  norytmp = norytmp(:,ones(Nint,1))';

  satmp = saExt(:);
  satmp = satmp(:,ones(Nint,1))';

  rho2 = (xtar-xsou).^2 + (ytar-ysou).^2;

  coeff = 1/pi*((xtar-xsou).*norxtmp + ...
    (ytar-ysou).*norytmp).*satmp./rho2.^2;

  D(1:Nint,:) = 2*pi/Next*[coeff.*(xtar-xsou).^2 ...
    coeff.*(xtar-xsou).*(ytar-ysou)];
  D(Nint+1:end,:) = 2*pi/Next*[coeff.*(ytar-ysou).*(xtar-xsou) ...
    coeff.*(ytar-ysou).^2];

  if o.fastDirect
    [rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(itar:jtar)',(1:2*Next)');
  else
    rowVectN = (itar:jtar)';
    colVectN = (1:2*Next)';
  end
  M11(rowVectN,colVectN) = D;
end %end ktar

% -----------------------------------------------------------------------
% INTERIOR WALLS

for ktar = 1:nv-1 % loop over targets
  itar = 2*Next + 2*(ktar-1)*Nint + 1;
  jtar = 2*Next + 2*ktar*Nint;
  K = [(1:ktar-1) (ktar+1:nv-1)];

  D = zeros(2*Nint,2*Nint);
  for ksou = K % loop over all other walls
    isou = 2*Next + 2*(ksou-1)*Nint + 1;
    jsou = 2*Next + 2*ksou*Nint;

    xtar = xInt(:,ktar); ytar = yInt(:,ktar);
    xtar = xtar(:,ones(Nint,1)); 
    ytar = ytar(:,ones(Nint,1));

    xsou = xInt(:,ksou); ysou = yInt(:,ksou);
    xsou = xsou(:,ones(Nint,1))';
    ysou = ysou(:,ones(Nint,1))';

    norxtmp = norxInt(:,ksou); norytmp = noryInt(:,ksou);
    norxtmp = norxtmp(:,ones(Nint,1))';
    norytmp = norytmp(:,ones(Nint,1))';

    satmp = saInt(:,ksou);
    satmp = satmp(:,ones(Nint,1))';

    rho2 = (xtar-xsou).^2 + (ytar-ysou).^2;

    coeff = 1/pi*((xtar-xsou).*norxtmp + ...
      (ytar-ysou).*norytmp).*satmp./rho2.^2;

    D(1:Nint,:) = 2*pi/Nint*[coeff.*(xtar-xsou).^2 ...
      coeff.*(xtar-xsou).*(ytar-ysou)];
    D(Nint+1:end,:) = 2*pi/Nint*[coeff.*(ytar-ysou).*(xtar-xsou) ...
      coeff.*(ytar-ysou).^2];
    
    if o.fastDirect
      [rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(itar:jtar)',(isou:jsou)');
    else
      rowVectN = (itar:jtar)';
      colVectN = (isou:jsou)';
    end
    M11(rowVectN,colVectN) = D;
    
  end %end ktar
end %end ksou
% ------------------------------------------------------------------------

% These compute the integral of the density function around each of the
% inner componenents of the geometry
for k = 1:nv-1
  icol = 3*(k-1)+1;
  istart = 2*Next+2*(k-1)*Nint+1;
  iend = istart + Nint - 1;
  
  if o.fastDirect
    [~,colVectN] = o.findNewIds(newIDs,oldIDs,1,(istart:iend)');
  else
    colVectN = (istart:iend)';
  end
  M21(icol,colVectN) = 2*pi/Nint*saInt(:,k)';
  M21(icol+2,colVectN) = 2*pi/Nint*saInt(:,k)'.*yInt(:,k)';
  
  istart = istart + Nint;
  iend = iend + Nint;
  if o.fastDirect
    [~,colVectN] = o.findNewIds(newIDs,oldIDs,1,(istart:iend)');
  else
    colVectN = (istart:iend)';
  end
  M21(icol+1,colVectN) = 2*pi/Nint*saInt(:,k)';
  M21(icol+2,colVectN) = -2*pi/Nint*saInt(:,k)'.*xInt(:,k)';
end % k



% Interior walls on Exterior wall
for k = 1:nv - 1
  rho2 = (xExt(:) - cxInt(k)).^2 + (yExt(:) - cyInt(k)).^2;
  istart = 1;
  iend = Next;

  icol = 3*(k-1)+1;
  if o.fastDirect
    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
  else
    rowVectN = (istart:iend)';
  end
  M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*(-0.5*log(rho2) + (xExt(:)-cxInt(k))./rho2.*...
        (xExt(:)-cxInt(k)));
    
  if o.fastDirect
    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Next:iend+Next)',1);  
  else
    rowVectN = (istart+Next:iend+Next)';
  end  
  M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*((xExt(:)-cxInt(k))./rho2.*(yExt(:)-cyInt(k)));

  icol = 3*(k-1)+2;
  if o.fastDirect
    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
  else
    rowVectN = (istart:iend)';
  end
  M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*((yExt(:)-cyInt(k))./rho2.*(xExt(:)-cxInt(k)));
  
  if o.fastDirect
    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Next:iend+Next)',1);  
  else
    rowVectN = (istart+Next:iend+Next)';
  end
  M12(rowVectN,icol) = M12(rowVectN,icol) + ...
    1/4/pi*(-0.5*log(rho2) + (yExt(:)-cyInt(k))./rho2.*...
      (yExt(:)-cyInt(k)));

  icol = 3*(k-1)+3;
  if o.fastDirect
    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
  else
    rowVectN = (istart:iend)';
  end
  M12(rowVectN,icol) = M12(rowVectN,icol) + (yExt(:)-cyInt(k))./rho2;
  if o.fastDirect
    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Next:iend+Next)',1);  
  else
    rowVectN = (istart+Next:iend+Next)';
  end
  M12(rowVectN,icol) = M12(rowVectN,icol) - (xExt(:)-cxInt(k))./rho2;  
end

% Interior walls on each other
for k = 1:nv - 1
  for ktar = 1:nv-1
    rho2 = (xInt(:,ktar) - cxInt(k)).^2 + (yInt(:,ktar) - cyInt(k)).^2;
    istart = 2*Next+(ktar-1)*2*Nint + 1;
    iend = istart + Nint - 1;

    icol = 3*(k-1)+1;
    if o.fastDirect
      [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
    else
      rowVectN = (istart:iend)';
    end
    M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*(-0.5*log(rho2) + (xInt(:,ktar)-cxInt(k))./rho2.*...
        (xInt(:,ktar)-cxInt(k)));
    
    if o.fastDirect
      [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Nint:iend+Nint)',1);   
    else
      rowVectN = (istart+Nint:iend+Nint)';
    end
    M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*((xInt(:,ktar)-cxInt(k))./rho2.*(yInt(:,ktar)-cyInt(k)));

    icol = 3*(k-1)+2;
    if o.fastDirect
      [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
    else
      rowVectN = (istart:iend)';
    end
    M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*((yInt(:,ktar)-cyInt(k))./rho2.*(xInt(:,ktar)-cxInt(k)));
  
    if o.fastDirect
      [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Nint:iend+Nint)',1);   
    else
      rowVectN = (istart+Nint:iend+Nint)';
    end 
    M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*(-0.5*log(rho2) + (yInt(:,ktar)-cyInt(k))./rho2.*...
        (yInt(:,ktar)-cyInt(k)));

    icol = 3*(k-1)+3;
    if o.fastDirect
      [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
    else
      rowVectN = (istart:iend)';
    end 
    M12(rowVectN,icol) = M12(rowVectN,icol) + (yInt(:,ktar)-cyInt(k))./rho2;
    
    if o.fastDirect
      [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Nint:iend+Nint)',1);
    else
      rowVectN = (istart+Nint:iend+Nint)';
    end 
    M12(rowVectN,icol) = M12(rowVectN,icol) - (xInt(:,ktar)-cxInt(k))./rho2;
  end
end

% different combinations of the density functions have to multiply to
% 2*pi multiplied by rotlet or stokeslet terms
M22 = -2*pi*eye(3*(nv-1));

% save DLP and RS matrix so that we do not compute them in TimeMatVec
% because they do not change with time and computing interior
% walls-to-interior walls potentials take most of the time during the
% simulations.
if ~o.fastDirect % if not fastDirect, then invert the matrix exactly
  Mat = ([M11 M12; M21 M22])\eye(2*Next+2*(nv-1)*Nint + 3*(nv-1));
  %message = ['DONE, it took ' num2str(toc(tForm0),'%2.2e') ' seconds'];
 % o.om.writeMessage(message,'%s\n');
 % o.om.writeMessage(' ','%s\n');
  
  if ~o.matFreeWalls
    o.wallDLPandRSmat = [M11 M12; M21 M22];
  end
  
  if o.saveWallMat
    fileName = [o.wallMatFile 'Inv.bin'];
    fid = fopen(fileName,'w');
    fwrite(fid,Mat(:),'double');
    fclose(fid);
    message = ['Exact factorization was saved on disk'];
    o.om.writeMessage(message,'%s\n');
    o.om.writeMessage(' ','%s\n');
  
    if ~o.matFreeWalls
      fileName = [o.wallMatFile 'DLPandRS.bin'];
      fid = fopen(fileName,'w');
      fwrite(fid,o.wallDLPandRSmat(:),'double');
      fclose(fid);
      message = ['Walls DLP and RS matrix was saved on disk'];
      o.om.writeMessage(message,'%s\n');
      o.om.writeMessage(' ','%s\n');
    end
  end % o.saveWallMat
else % then it is fastDirect
  message = ['DONE, it took ' num2str(toc(tForm0),'%2.2e') ' seconds'];
  o.om.writeMessage(message,'%s\n');
  o.om.writeMessage(' ','%s\n');  
  
  % if we also want to compress wall2wall interactions
  if o.HODLRforW2W
    message = ['Compressing wall2wall interaction matrix via HODLR...'];
    tComp = tic;
    o.om.writeMessage(message,'%s\n');
  
    hoWalls = o.hodlrWalls; 
    hoWalls.L = M11;
    hoWalls = hoWalls.simple_partition(2*Next+2*Nint*(nv-1));
    hoWalls = hoWalls.hodlr_no_FMM(o.wallMatFile); 
    hoWalls.idWall2gids = idWall2gids;
    hoWalls.oldIDs = oldIDs;
    hoWalls.newIDs = newIDs;
    
    o.wallDLPandRSmatM12 = M12;
    o.wallDLPandRSmatM21 = M21;
    o.wallDLPandRSmatM22 = M22;
  
    o.hodlrWalls = hoWalls;       
    message = ['DONE, it took ' num2str(toc(tComp),'%2.2e') ' seconds'];
    o.om.writeMessage(message,'%s\n');
    o.om.writeMessage(' ','%s\n');  
  end

  % save the blocks of the wall2wall interactions on disc
  fid = fopen(fileNameM11,'w');
  fwrite(fid,M11,'double');
  fclose(fid);

  fid = fopen(fileNameM12,'w');
  fwrite(fid,M12,'double');
  fclose(fid);

  fid = fopen(fileNameM21,'w');
  fwrite(fid,M21,'double');
  fclose(fid);

  fid = fopen(fileNameM22,'w');
  fwrite(fid,M22,'double');
  fclose(fid);

  ho.idWall2gids = idWall2gids;

  memmapM11 = memmapfile(fileNameM11,'Format',{'double',size(M11),'M'});
  memmapM12 = memmapfile(fileNameM12,'Format',{'double',size(M12),'M'});
  ho.memmapM21 = memmapfile(fileNameM21,'Format',{'double',size(M21),'M'});
  ho.memmapM22 = memmapfile(fileNameM22,'Format',{'double',size(M22),'M'});

  Mat = [];
  clear M11; 
  clear M12;
  clear M21;
  clear M22;
  
  tSolve0 = tic; 
  message = ['Compressing and factorizing wall2wall interaction matrix with HODLR...'];
  o.om.writeMessage(message,'%s\n');

  ho = ho.solveReqs(memmapM11,memmapM12,...
    2*Next+2*Nint*(nv-1),3*(nv-1),o.wallMatFile,1);
  ho.oldIDs = oldIDs;
  ho.newIDs = newIDs;
  o.hodlr = ho;
  
  tSolve = toc(tSolve0); 
  message = ['DONE, it took ' num2str(tSolve,'%2.2e') ' seconds'];
  o.om.writeMessage(message,'%s\n');
  o.om.writeMessage(' ','%s\n');
end

end % wallsPrecondDiffDisc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compressWall2Wall(o,wallsInt,wallsExt)
% compressWall2Wall(o,wallsInt,wallsExt) compresses wall2wall interaction 
% matrix using HODLR.

% discretization info
Next = wallsExt.N;
Nint = wallsInt.N;
nv = wallsInt.nv+1;

% if we are going to use HODLR, then we need to order the points using
% Z-curve
% That's why, we change the pointers and build the ordered matrix

% Shuffling order
ho = o.hodlrWalls;
idWall2gids = ho.pnts2qtree(wallsInt.X,wallsExt.X,0);


% Get the old and new ids of M11, M12, M21, M22
newIDs = [1:Next+Nint*(nv-1)]';
pntIDs = idWall2gids(newIDs,1);
wallIDs = idWall2gids(newIDs,2);

extWallIDs = find(wallIDs==1);
intWallIDs = find(wallIDs~=1);

oldIDs = zeros(size(newIDs));
oldIDs(extWallIDs) = pntIDs(extWallIDs);
oldIDs(intWallIDs) = pntIDs(intWallIDs)+2*Nint*(wallIDs(intWallIDs)-2)+2*Next;

oldIDsY = zeros(size(newIDs));
oldIDsY(extWallIDs) = oldIDs(extWallIDs)+Next;
oldIDsY(intWallIDs) = oldIDs(intWallIDs)+Nint;

oldIDs = [oldIDs;oldIDsY];
newIDs = [2*newIDs-1;2*newIDs];

tForm0 = tic;
message = ['Building wall2wall interaction matrix...'];
o.om.writeMessage(message,'%s\n');    

oc = curve;

[xInt,yInt] = oc.getXY(wallsInt.X);
[xExt,yExt] = oc.getXY(wallsExt.X);

[noryInt,norxInt] = oc.getXY(wallsInt.xt);
noryInt = -noryInt;
saInt = wallsInt.sa;
[cxInt,cyInt] = oc.getXY(wallsInt.center);

[noryExt,norxExt] = oc.getXY(wallsExt.xt);
noryExt = -noryExt;
saExt = wallsExt.sa;

% Allocate space for blocks of matrix that carries the double- layer
% potential, rotlets, and stokeslets to the velocity and the conditions
% in (A4) and (A5) in Rahimian et al.
M11 = zeros(2*Nint*(nv-1)+2*Next);
M12 = zeros(2*Nint*(nv-1)+2*Next,3*(nv-1));
M21 = zeros(3*(nv-1),2*Nint*(nv-1)+2*Next);


% Self interaction terms with the jump coming from the double layer
% potential
[rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(1:2*Next)',(1:2*Next)');

M11(rowVectN,colVectN) = M11(rowVectN,colVectN) + o.wallN0(:,:,1);
jump = - 1/2;
M11(rowVectN,colVectN) = M11(rowVectN,colVectN) + ...
jump*eye(2*Next) + o.wallDLPext;

for k = 1:nv-1
  istart = 2*Next + (k-1)*2*Nint+1;
  iend = istart + 2*Nint-1;
  [rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(istart:iend)',(istart:iend)');
  
  M11(rowVectN,colVectN) = M11(rowVectN,colVectN) + ...
    jump*eye(2*Nint) + o.wallDLPint(:,:,k);
end


% -----------------------------------------------------------------------
% EXTERIOR WALL
% Do ktar = 1, ksou = 2 : nv and 
% ktar = 2 : nv, ksou = 1 separately and then the rest
% in the old loop

for ksou = 1:nv-1 % loop over all other walls
  isou = 2*Next + 2*(ksou-1)*Nint + 1;
  jsou = 2*Next + 2*ksou*Nint;

  xtar = xExt(:); ytar = yExt(:);
  xtar = xtar(:,ones(Nint,1)); 
  ytar = ytar(:,ones(Nint,1));

  xsou = xInt(:,ksou); ysou = yInt(:,ksou);
  xsou = xsou(:,ones(Next,1))';
  ysou = ysou(:,ones(Next,1))';

  norxtmp = norxInt(:,ksou); norytmp = noryInt(:,ksou);
  norxtmp = norxtmp(:,ones(Next,1))';
  norytmp = norytmp(:,ones(Next,1))';

  satmp = saInt(:,ksou);
  satmp = satmp(:,ones(Next,1))';

  rho2 = (xtar-xsou).^2 + (ytar-ysou).^2;

  coeff = 1/pi*((xtar-xsou).*norxtmp + ...
    (ytar-ysou).*norytmp).*satmp./rho2.^2;

  
  [rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(1:Next)',(isou:jsou)');
  M11(rowVectN,colVectN) = 2*pi/Nint*[coeff.*(xtar-xsou).^2 ...
    coeff.*(xtar-xsou).*(ytar-ysou)];

  [rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(Next+1:2*Next)',(isou:jsou)');
  M11(rowVectN,colVectN) = 2*pi/Nint*[coeff.*(ytar-ysou).*(xtar-xsou) ...
    coeff.*(ytar-ysou).^2];

end %end ksou

for ktar = 1:nv-1 % loop over targets
  itar = 2*Next + 2*(ktar-1)*Nint + 1;
  jtar = 2*Next + 2*ktar*Nint;

  D = zeros(2*Nint,2*Next);

  xtar = xInt(:,ktar); ytar = yInt(:,ktar);
  xtar = xtar(:,ones(Next,1)); 
  ytar = ytar(:,ones(Next,1));

  xsou = xExt(:); ysou = yExt(:);
  xsou = xsou(:,ones(Nint,1))';
  ysou = ysou(:,ones(Nint,1))';

  norxtmp = norxExt(:); norytmp = noryExt(:);
  norxtmp = norxtmp(:,ones(Nint,1))';
  norytmp = norytmp(:,ones(Nint,1))';

  satmp = saExt(:);
  satmp = satmp(:,ones(Nint,1))';

  rho2 = (xtar-xsou).^2 + (ytar-ysou).^2;

  coeff = 1/pi*((xtar-xsou).*norxtmp + ...
    (ytar-ysou).*norytmp).*satmp./rho2.^2;

  D(1:Nint,:) = 2*pi/Next*[coeff.*(xtar-xsou).^2 ...
    coeff.*(xtar-xsou).*(ytar-ysou)];
  D(Nint+1:end,:) = 2*pi/Next*[coeff.*(ytar-ysou).*(xtar-xsou) ...
    coeff.*(ytar-ysou).^2];

  [rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(itar:jtar)',(1:2*Next)');
  M11(rowVectN,colVectN) = D;
end %end ktar

% -----------------------------------------------------------------------
% INTERIOR WALLS

for ktar = 1:nv-1 % loop over targets
  itar = 2*Next + 2*(ktar-1)*Nint + 1;
  jtar = 2*Next + 2*ktar*Nint;
  K = [(1:ktar-1) (ktar+1:nv-1)];

  D = zeros(2*Nint,2*Nint);
  for ksou = K % loop over all other walls
    isou = 2*Next + 2*(ksou-1)*Nint + 1;
    jsou = 2*Next + 2*ksou*Nint;

    xtar = xInt(:,ktar); ytar = yInt(:,ktar);
    xtar = xtar(:,ones(Nint,1)); 
    ytar = ytar(:,ones(Nint,1));

    xsou = xInt(:,ksou); ysou = yInt(:,ksou);
    xsou = xsou(:,ones(Nint,1))';
    ysou = ysou(:,ones(Nint,1))';

    norxtmp = norxInt(:,ksou); norytmp = noryInt(:,ksou);
    norxtmp = norxtmp(:,ones(Nint,1))';
    norytmp = norytmp(:,ones(Nint,1))';

    satmp = saInt(:,ksou);
    satmp = satmp(:,ones(Nint,1))';

    rho2 = (xtar-xsou).^2 + (ytar-ysou).^2;

    coeff = 1/pi*((xtar-xsou).*norxtmp + ...
      (ytar-ysou).*norytmp).*satmp./rho2.^2;

    D(1:Nint,:) = 2*pi/Nint*[coeff.*(xtar-xsou).^2 ...
      coeff.*(xtar-xsou).*(ytar-ysou)];
    D(Nint+1:end,:) = 2*pi/Nint*[coeff.*(ytar-ysou).*(xtar-xsou) ...
      coeff.*(ytar-ysou).^2];
    
    [rowVectN,colVectN] = o.findNewIds(newIDs,oldIDs,(itar:jtar)',(isou:jsou)');
    M11(rowVectN,colVectN) = D;
    
  end %end ktar
end %end ksou
% ------------------------------------------------------------------------

% These compute the integral of the density function around each of the
% inner componenents of the geometry
for k = 1:nv-1
  icol = 3*(k-1)+1;
  istart = 2*Next+2*(k-1)*Nint+1;
  iend = istart + Nint - 1;
  
  [~,colVectN] = o.findNewIds(newIDs,oldIDs,1,(istart:iend)');
  M21(icol,colVectN) = 2*pi/Nint*saInt(:,k)';
  M21(icol+2,colVectN) = 2*pi/Nint*saInt(:,k)'.*yInt(:,k)';
  
  istart = istart + Nint;
  iend = iend + Nint;
  [~,colVectN] = o.findNewIds(newIDs,oldIDs,1,(istart:iend)');
  
  M21(icol+1,colVectN) = 2*pi/Nint*saInt(:,k)';
  M21(icol+2,colVectN) = -2*pi/Nint*saInt(:,k)'.*xInt(:,k)';
end % k



% Interior walls on Exterior wall
for k = 1:nv - 1
  rho2 = (xExt(:) - cxInt(k)).^2 + (yExt(:) - cyInt(k)).^2;
  istart = 1;
  iend = Next;

  icol = 3*(k-1)+1;
  [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
  
  M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*(-0.5*log(rho2) + (xExt(:)-cxInt(k))./rho2.*...
        (xExt(:)-cxInt(k)));

  [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Next:iend+Next)',1);  
  M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*((xExt(:)-cxInt(k))./rho2.*(yExt(:)-cyInt(k)));

  icol = 3*(k-1)+2;
  [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
  M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*((yExt(:)-cyInt(k))./rho2.*(xExt(:)-cxInt(k)));

  [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Next:iend+Next)',1);  
  M12(rowVectN,icol) = M12(rowVectN,icol) + ...
    1/4/pi*(-0.5*log(rho2) + (yExt(:)-cyInt(k))./rho2.*...
      (yExt(:)-cyInt(k)));

  icol = 3*(k-1)+3;
  [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
  
  M12(rowVectN,icol) = M12(rowVectN,icol) + (yExt(:)-cyInt(k))./rho2;
  [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Next:iend+Next)',1);  
  
  M12(rowVectN,icol) = M12(rowVectN,icol) - (xExt(:)-cxInt(k))./rho2;  
end

% Interior walls on each other
for k = 1:nv - 1
  for ktar = 1:nv-1
    rho2 = (xInt(:,ktar) - cxInt(k)).^2 + (yInt(:,ktar) - cyInt(k)).^2;
    istart = 2*Next+(ktar-1)*2*Nint + 1;
    iend = istart + Nint - 1;

    icol = 3*(k-1)+1;
    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
    
    M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*(-0.5*log(rho2) + (xInt(:,ktar)-cxInt(k))./rho2.*...
        (xInt(:,ktar)-cxInt(k)));
    
    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Nint:iend+Nint)',1);   
    M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*((xInt(:,ktar)-cxInt(k))./rho2.*(yInt(:,ktar)-cyInt(k)));

    icol = 3*(k-1)+2;
    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
    
    M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*((yInt(:,ktar)-cyInt(k))./rho2.*(xInt(:,ktar)-cxInt(k)));

    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Nint:iend+Nint)',1);   
    M12(rowVectN,icol) = M12(rowVectN,icol) + ...
      1/4/pi*(-0.5*log(rho2) + (yInt(:,ktar)-cyInt(k))./rho2.*...
        (yInt(:,ktar)-cyInt(k)));

    icol = 3*(k-1)+3;
    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart:iend)',1);
    M12(rowVectN,icol) = M12(rowVectN,icol) + (yInt(:,ktar)-cyInt(k))./rho2;

    [rowVectN,~] = o.findNewIds(newIDs,oldIDs,(istart+Nint:iend+Nint)',1);
    M12(rowVectN,icol) = M12(rowVectN,icol) - (xInt(:,ktar)-cxInt(k))./rho2;
  end
end

% different combinations of the density functions have to multiply to
% 2*pi multiplied by rotlet or stokeslet terms
M22 = -2*pi*eye(3*(nv-1));

% save DLP and RS matrix so that we do not compute them in TimeMatVec
% because they do not change with time and computing interior
% walls-to-interior walls potentials take most of the time during the
% simulations.
message = ['DONE, it took ' num2str(toc(tForm0),'%2.2e') ' seconds'];
o.om.writeMessage(message,'%s\n');
o.om.writeMessage(' ','%s\n');  
  
  
message = ['Compressing wall2wall interaction matrix via HODLR...'];
tComp = tic;
o.om.writeMessage(message,'%s\n');

ho = o.hodlrWalls; 
ho.L = M11;
ho = ho.simple_partition(2*Next+2*Nint*(nv-1));
ho = ho.hodlr_no_FMM(o.wallMatFile); 
ho.idWall2gids = idWall2gids;
ho.oldIDs = oldIDs;
ho.newIDs = newIDs;

o.wallDLPandRSmatM12 = M12;
o.wallDLPandRSmatM21 = M21;
o.wallDLPandRSmatM22 = M22;

o.hodlrWalls = ho;       
message = ['DONE, it took ' num2str(toc(tComp),'%2.2e') ' seconds'];
o.om.writeMessage(message,'%s\n');
o.om.writeMessage(' ','%s\n');  
 

end % compressWall2Wall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = applyWallsPrecondFD(o,prams,rhs)
% applyWallsPrecondFD(prams,rhs) solves for density and rotlets, stokeslets
% of the walls of which the matrix is built in wallsPrecondFD. Here,
% entries of rhs are ordered in the same way the matrix is ordered, this is
% necessary for getting better accuracy from hodlr.m. Then the system is
% solved using the schur complement. 
ho = o.hodlr;

% Order the entries of rhs as well (only ones corresponding to DLP)
Nint = prams.NbdInt;
Next = prams.NbdExt;
nv   = prams.nvbd;
k = 3*(nv-1);
n = 2*Next+2*Nint*(nv-1);
    
rhsNew = zeros(numel(rhs),1);
rhsNew(2*Next+2*Nint*(nv-1)+1:end) = ...
  rhs(2*Next+2*Nint*(nv-1)+1:end);

% Get the order
newIDs = ho.newIDs;
oldIDs = ho.oldIDs;

% order the rhs
rhsNew(newIDs,1) = rhs(oldIDs,1);

% solve the system
lhs = ho.schur_complement(rhsNew,n,k);

% Reorder the result so that the form is [etaX_wall1;etaY_wall1;...
% etaX_wallN;etaY_wallN;RSwall1;RSwall2;...]
u = zeros(2*Nint*(nv-1)+2*Next,1);
u(oldIDs,1) = lhs(newIDs,1);


z = [u;lhs(2*Next+2*Nint*(nv-1)+1:end)];
    
end %applyWallsPrecondFD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = applyWall2WallFD(o,Next,Nint,nv,lhs)
% applyWall2WallFD(o,Next,Nint,nv,lhs) is similar to applyWallsPrecondFD,
% but it applies the compressed wall2wall interaction matrix to the vector
% (lhs) and gets the interactions. This is used when wall2wall interactions
% are handled using the precomputed wall2wall interactions (i.e.
% matFreeWalls = false)

ho = o.hodlrWalls;
n = 2*Next+2*Nint*(nv-1);
    
lhsNew = zeros(size(lhs));
lhsNew(n+1:end) = lhs(n+1:end);

% Get the order (already computed)
newIDs = ho.newIDs;
oldIDs = ho.oldIDs;

% order the lhs
lhsNew(newIDs,1) = lhs(oldIDs,1);

% apply the matrix 
z1 = ho.mult(lhsNew(1:n),3,n,1) + o.wallDLPandRSmatM12*lhsNew(n+1:end);
z2 = o.wallDLPandRSmatM21*lhsNew(1:n)+o.wallDLPandRSmatM22*lhsNew(n+1:end);
z = [z1;z2];


% Reorder the result so that the form is [etaX_wall1;etaY_wall1;...
% etaX_wallN;etaY_wallN;RSwall1;RSwall2;...]
u = zeros(n,1);
u(oldIDs,1) = z(newIDs,1);
z = [u;z(n+1:end)];

end %applyWall2WallFD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = precoIminusD(o,z,vesicle)
% z = o.precoIminusD(z,vesicle) is the block-diagonal preconditioner
% for the system alpha*I + DLP.  This did not significantly reduce the
% number of GMRES steps, so it is currently not being called

N = vesicle.N;
nv = vesicle.nv;
alpha = 0.5*(1+vesicle.viscCont);
val = zeros(2*N*nv,1);

for k = 1:nv
  istart = (k-1)*2*N + 1;
  iend = k*2*N;
  z(istart:iend) = (alpha(k)*eye(2*N) - o.D(:,:,k))\z(istart:iend);
end

end % precoIminusD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerTen(o,z)
% val = preconditionerTen(z) applies the preconditioner to the tension
% term for when we are solving for just the tension and density
% function given a configuration.  Configuration is eliminated using
% the Schur complelent

% number of points per vesicle
nv = size(o.bdiagTen,3); % number of vesicles
N = size(o.bdiagTen,1); 

% number of points per solid wall
nvbd = size(o.wallDLP,3); % number of solid walls
Nbd = size(o.wallDLP,1)/2; 

% part of z correpsonding to the vesicles
zves = z(1:N*nv);
% part of z corresponding to the solid walls
zwall = z(N*nv+1:end-3*(nvbd-1));
% part of z corresonding to rotlets and stokeslets
zrot = z(end-3*(nvbd-1)+1:end);

valVes = zeros(N*nv,1);
for k=1:nv
  valVes((k-1)*N+1:k*N) = o.bdiagTen(:,:,k)*...
    zves((k-1)*N+1:k*N);
end

valWall = o.bdiagWall*[zwall;zrot];

% stack the preconditioned values due to the tension term and
% the solid wall term
val = [valVes;valWall];


end % preconditionerTen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerSpectral(o,vesicle,walls,wallsCoarse,z)
% val = preconditionerSpectral(vesicle,walls,wallsCoarse,z) applies the
% block diagonal preconditioner where each block is the inverse
% corresponding to a circle of the same length as the vesicle

% HAVE NOT BEEN UPDATED (SEPTEMBER 19, 2017) FOR WALLS WITH DIFFERENT DISC.

N = vesicle.N;
nv = vesicle.nv;
if o.confined
  Nbd = walls.N;
  nvbd = walls.nv;
  valWall = zeros(2*Nbd,nvbd);
end
bx = zeros(2*N,nv);
bsig = zeros(N,nv);
for k = 1:nv
  istart = (k-1)*3*N+1;
  iend = istart + 2*N - 1;
  bx(:,k) = z(istart:iend);
  istart = iend + 1;
  iend = istart + N - 1;
  bsig(:,k) = z(istart:iend);
end
% seperate the right hand side terms due to tension and position

if strcmp(o.solver,'method1')
%  fprintf('Preconditioner not yet developed for this solver\n')
  bschur = bsig - vesicle.surfaceDiv(o.invIplusSB(vesicle,bx));

  sig = o.invSchur11(vesicle,bschur);
  % solve the Schur equation using gmres to find the tension

  x = vesicle.tensionTerm(sig);
  for k = 1:nv
    x(:,k) = o.Galpert(:,:,k)*x(:,k);
  end
  x = bx + o.dt * x;
  % build the right hand side term for solving fot the positon

  x = o.invIplusSB(vesicle,x);
  % solve for the position

elseif strcmp(o.solver,'method2')
  bschur = o.invDST(vesicle,bsig);
  bschur = vesicle.tensionTerm(bschur);
  Sbschur = zeros(2*N,nv);
  for k = 1:nv
    Gbschur(:,k) = o.Galpert(:,:,k) * bschur(:,k);
  end
  bschur = bx + o.dt*Gbschur;
  % right-hand side of the Schur equation involving the Schur complement

  x = o.invSchur22(vesicle,bschur);
  % solve the schur equation using gmres to find the position

  Benx = -vesicle.bendingTerm(x);
  for k = 1:nv
    Benx(:,k) = o.Galpert(:,:,k)*Benx(:,k);
  end
  sigRHS = bsig + vesicle.surfaceDiv(Benx);
  % build the right hand side term for solving for the tension
  sig = o.invDST(vesicle,sigRHS);
  % solve for the tension
end
% solved for position and tension


valVes = zeros(3*N*nv,1);
for k = 1:nv
  istart = (k-1)*3*N+1;
  iend = istart + 3*N - 1;
  valVes(istart:iend) = [x(:,k);sig(:,k)];
end
% stack the position and tension as required by the outermost matvec
% routine

if o.confined
  for k = 1:nvbd
    istart = 3*N*nv + (k-1)*2*Nbd + 1;
    iend = istart + 2*Nbd - 1;
    valWall(:,k) = z(istart:iend);
  end
  valLets = z(iend+1:end);
  nvcycle = 1;
  for k = 1:nvcycle
    [valWall,valLets] = o.vcycle(walls,wallsCoarse,valWall,valLets);
  end
  valWall = valWall(:);
else
  valWall = [];
  valLets = [];
end

val = [valVes;valWall;valLets];
% stack the vesicle data followed by the solid wall data followed by
% the rotlets and stokeslets

end % preconditionerSpectral

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [valWalls,valLets] = vcycle(o,walls,wallsCoarse,zWalls,zLets)
% [valWalls,valLets] = vcycle(walls,wallsCoarse,zWalls,zLets) applies
% two-grid V-cycle to the integral equation defined on the solid
% walls.  Post-smoothing seems to not work correctly so use V(1,0)
% cicyles seems best.  Coarse grid solve is done with the precomputed
% double-layer potential on the coarse grid

% HAVE NOT BEEN UPDATED (SEPTEMBER 19, 2017) FOR WALLS WITH DIFFERENT DISC.

N = walls.N;
Ncoarse = wallsCoarse.N; 
nv = walls.nv;
resWalls = zeros(2*N,nv);
valWalls = zeros(2*N,nv);
errWalls = zeros(2*N,nv);
resCoarse = zeros(2*Ncoarse,nv);
errCoarse = zeros(2*Ncoarse,nv);

op = o.op;

valWalls = zeros(2*N,nv);
valLets = zeros(3*(nv-1),1);
% initial guess

for npre = 1:1
  if norm(valWalls) == 0
    valWalls = -2*zWalls;
  else
    valWalls = -2*(zWalls - op.exactStokesDLdiag(walls,o.wallDLP,valWalls));
  end
  % Picard smoother applied to the density function
%  residual = ([zWalls(:);zLets] - ...
%     o.denMatVec([valWalls(:);valLets],walls,[]));
end
% pre-smoothing

residual = ([zWalls(:);zLets] - ...
    o.denMatVec([valWalls(:);valLets],walls,[]));
%clf;
%fprintf('After pre-smoothing\n')

for k = 1:nv
  istart = (k-1)*2*N + 1;
  iend = istart + 2*N - 1;
  resWalls(:,k) = residual(istart:iend);
end
resLets = residual(iend+1:end);
% form the residual

resCoarse = curve.restrict(resWalls,N,Ncoarse);
% restrict the residual

rhsCoarse = [resCoarse(:);resLets];

err = o.bdiagWall*rhsCoarse;
% solve for the coarse grid error

for k = 1:nv
  istart = (k-1)*2*Ncoarse + 1;
  iend = istart + 2*Ncoarse - 1;
  errCoarse(:,k) = err(istart:iend);
end
errLets = err(iend+1:end);

errWalls = curve.prolong(errCoarse,Ncoarse,N);
% prolong the error

valWalls = valWalls + errWalls;
valLets = valLets + errLets;
% add in the error 

for npost = 1:0
  valWalls = -2*(zWalls - ...
      op.exactStokesDLdiag(walls,o.wallDLP,valWalls));
  % Picard smoother applied to the density function
end
% post-smoothing step

end % vcycle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = invSchur22(o,vesicle,z)
% x = invSchur22(vesicle,z) solves the schur complement equation S*x =
% z where S is the operator I + dt*G*Ben - dt*S*Ten*invDST*Div*G*Ben

N = vesicle.N; nv = vesicle.nv;

[z,flag,relres,iter,relresvec] = gmres(@(X) ...
    o.Schur22MatVec(vesicle,X),...
    z(:),[],o.gmresTol,min(o.gmresMaxIter,2*N));

x = zeros(2*N,nv);
for k = 1:nv
  x(:,k) = z((k-1)*2*N+1:k*2*N);
end

end % invSchur22

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = Schur22MatVec(o,vesicle,x)
% val = Schur22MatVec(vesicle,x) applies the Schur complement operator
% to x.  Schur complement is given by I + dt*G*Ben -
% dt*S*Ten*invDST*Div*G*Ben

N = vesicle.N; nv = vesicle.nv;
val = zeros(2*N,nv);
xCols = zeros(2*N,nv);
for k = 1:nv
  xCols(:,k) = x((k-1)*2*N+1:k*2*N); 
end
x = -vesicle.bendingTerm(xCols);
for k = 1:nv
  x(:,k) = o.Galpert(:,:,k) * x(:,k);
end

sig = o.invDST(vesicle,vesicle.surfaceDiv(x));
Tsig = vesicle.tensionTerm(sig);
for k = 1:nv
  val(:,k) = o.Galpert(:,:,k) * Tsig(:,k);
end
val = xCols + o.dt*(x - val);

val = val(:);

end % Schur22MatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = precoSchur22(o,vesicle,alpha,scaling,x)
% val = precoSchur22(vesicle,alpha,scaling,x) solves the equation (I +
% dt*kappa*SLP*Ben)*val = x.  This can be used as a preconditioner when
% doing the Schur complement of the (1,1) block in method1 or when
% doing the Schur complement of the (2,2) block in method2
% TODO: This is wrong.  Should do the analytic preconditioner of the
% actual Schur block.

N = vesicle.N; nv = vesicle.nv;
val = zeros(2*N,nv);

for k = 1:nv
  val(:,k) = x((k-1)*2*N+1:k*2*N);
end
val(1:N,:) = fft(val(1:N,:));
val(N+1:2*N,:) = fft(val(N+1:2*N,:));

for k = 1:nv
  val(3:N-1,k) = val(3:N-1,k).*scaling(3:N-1);
  val(N+3:2*N-1,k) = val(N+3:2*N-1,k).*scaling(3:N-1);
end
% all modes with modulus greater than 1 are diagonal

plusOneMode = [val(2,:) ; val(N+2,:)];
minusOneMode = [val(N,:) ; val(2*N,:)];
% need to save these values since the plus and minus one modes
% communicate ie. isn't diagonal

val(N,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  ((2*alpha^2-4*alpha+1)*val(N,:) - ...
  2*1i*alpha^2*val(2*N,:)) + ...
  alpha/(4*alpha-1) * ...
  (plusOneMode(1,:) + 1i*plusOneMode(2,:));
% -1 mode of the x component

val(2*N,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  (2*1i*alpha^2*val(N,:) + ...
  (2*alpha^2-4*alpha+1)*val(2*N,:)) + ...
  alpha/(4*alpha-1) * ...
  (1i*plusOneMode(1,:) - plusOneMode(2,:));
% -1 mode of the y component

val(2,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  ((2*alpha^2-4*alpha+1)*val(2,:) + ...
  2*1i*alpha^2*val(N+2,:)) + ...
  alpha/(4*alpha-1) * ...
  (minusOneMode(1,:) - 1i*minusOneMode(2,:));
% 1 mode of the x component

val(N+2,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  (-2*1i*alpha^2*val(2,:) + ...
  (2*alpha^2-4*alpha+1)*val(N+2,:)) + ...
  alpha/(4*alpha-1) * ...
  (-1i*minusOneMode(1,:) - minusOneMode(2,:));
% 1 mode of the y component


val(1:N,:) = real(ifft(val(1:N,:)));
val(N+1:2*N,:) = real(ifft(val(N+1:end,:)));
val = val(:);

end % precoSchur22

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = invDST(o,vesicle,z)
% sig = invDST(o,vesicle,z) solves the equation DST*sig = z using pGMRES

N = vesicle.N; nv = vesicle.nv;

[z,flag,relres,iter] = gmres(@(X) o.DSTMatVec(vesicle,X),...
    z(:),[],1e-2*o.gmresTol,min(o.gmresMaxIter,N),...
    @(X) o.precoDST(vesicle,X));
% need a smaller tolerance here to hide that fact that this does not
% make a liner preconditioner.  Krylov stuff breaks down in theory, but
% this shouldn't come up until the error of the outer iteration is
% smaller than the requested tolerance

sig = zeros(N,nv);
for k = 1:nv
  sig(:,k) = z((k-1)*N+1:k*N);
end

end % invDST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = DSTMatVec(o,vesicle,sig)
% sig = DSTMatVec(vesicle,sig) applies the operator Div*SLP*Tension to
% sig

N = vesicle.N; nv = vesicle.nv;
sigCols = zeros(N,nv);
for k = 1:nv
  sigCols(:,k) = sig((k-1)*N+1:k*N); 
end

tension = vesicle.tensionTerm(sigCols);
% compute Tension * sig
for k = 1:nv
  tension(:,k) = o.Galpert(:,:,k)*tension(:,k);
end
% apply the single-layer potential
sig = vesicle.surfaceDiv(tension);
% compute the surface divergence
sig = sig(:);

end % DSTMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = precoDST(o,vesicle,sig)
% sig = precoDST(vesicle,sig) inverts the operator Div*SLP*Tension where
% the geometry is assumed to be a circle that has the same radius as the
% vesicle

N = vesicle.N; nv = vesicle.nv;
rad = vesicle.length/2/pi;
% radius of circle with same length as the vesicles

imodes = 1./[(1:N/2-1)';(N/2:-1:1)'];
for k = 1:nv
  sig((k-1)*N+1:k*N) = fft(sig((k-1)*N+1:k*N));
  sig((k-1)*N+2:k*N) = -4*rad*sig((k-1)*N+2:k*N).*imodes;
  sig((k-1)*N+1:k*N) = ifft(sig((k-1)*N+1:k*N));
end

sig = real(sig);

end % precoDST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = invSchur11(o,vesicle,z)
% sig = invSchur11(vesicle,z) solves the schur complement equation S*x
% = z where S is the operator defined by Schur11MatVec using pGMRES

N = vesicle.N; nv = vesicle.nv;

[z,flag,relres,iter,relresvec] = gmres(@(X) ...
    o.Schur11MatVec(vesicle,X),...
    z(:),[],o.gmresTol,min(o.gmresMaxIter,N),...
    @(X) o.precoSchur11(vesicle,X));

sig = zeros(N,nv);
for k = 1:nv
  sig(:,k) = z((k-1)*N+1:k*N);
end

end % invSchur22

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = Schur11MatVec(o,vesicle,sig)
% val = Schur11MatVec(vesicle,sig) applies the operator
% dt*Div*inv(eye+dt*SLP*Bending)*SLP*TEN which is the Schur complement
% of the (1,1) component
% TODO: THIS IS WRONG.  THE OPERATOR TO SOLVE SHOULD BE DIV*SLP*TEN -
% dt*DIV*SLP*BEN*inv(eye+dt*SLP*BEN)*SLP*TEN

N = vesicle.N; nv = vesicle.nv;
sigCols = zeros(N,nv);

for k = 1:nv
  sigCols(:,k) = sig((k-1)*N+1:k*N); 
end

sigCols = vesicle.tensionTerm(sigCols);
% This changes the size of sigCols to 2*N x nv.  This is no longer the
% tension, but the tension operator applied to the tension
for k = 1:nv
  sigCols(:,k) = o.Galpert(:,:,k) * sigCols(:,k);
end

sigCols = o.invIplusSB(vesicle,sigCols);
sigCols = o.dt*vesicle.surfaceDiv(sigCols);

val = sigCols(:);

end % Schur11MatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = precoSchur11(o,vesicle,sig)
% sig = precoSchur11(vesicle,sig) applies the spectral preconditioner
% of the Schur complement matrix when using the Schur complement of
% method1

N = vesicle.N; nv = vesicle.nv;
modes = [(0:N/2-1)';(-N/2:-1)'];
rad = vesicle.length/2/pi;
% radius of circle with same length as the vesicles

sigColh = zeros(N,nv);
for k = 1:nv
  sigColh(:,k) = sig((k-1)*N+1:k*N);
end
sigColh = fft(sigColh);

alpha = -o.dt*vesicle.kappa/8/rad^3;

scaling = 1./(-1/8*(abs(modes+1)./(1-2*alpha*abs(modes+1).^3) + ...
         abs(modes-1)./(1-2*alpha*abs(modes-1).^3)));
scaling(1) = 1;
scaling(2) = -4*(1-16*alpha);
scaling(3) = 1./(1.25e-1/(2*alpha-1) - 3.75e-1/(1-54*alpha));
scaling(N) = -4*(1-16*alpha);
scaling(N-1) = 1./(1.25e-1/(2*alpha-1) - 3.75e-1/(1-54*alpha));

for k = 1:nv
  sigColh(:,k) = sigColh(:,k) .* scaling;
end

sigColh = real(ifft(sigColh));

sig = rad*sigColh(:)/o.dt;

end % precoSchur11

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = invIplusSB(o,vesicle,z)
% x = invIplusSB(vesicle,z) solves the equation (I + dt*SLP*Bending)x =
% z using pGMRES 

N = vesicle.N; nv = vesicle.nv;

rad = vesicle.length/2/pi;
% radius of circle with same length as the vesicles
alpha = -o.dt*vesicle.kappa/8/rad^3;
% re-occuring constant
scaling = 1./(1-2*alpha*[(0:N/2-1) (N/2:-1:1)]'.^3);
% scaling that happens to the modes from modulus greater than or equal
% to 2 in the preconditioner
[z,flag,relres,iter] = gmres(@(X) o.IplusSBMatVec(vesicle,X),...
    z(:),[],1e-2*o.gmresTol,min(o.gmresMaxIter,N),...
    @(X) o.precoIplusSB(vesicle,alpha,scaling,X));

x = zeros(2*N,nv);
for k = 1:nv
  x(:,k) = z(2*(k-1)*N+1:k*2*N);
end

end % invIplusSB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = IplusSBMatVec(o,vesicle,x)
% x = IplusSBMatVec(vesicle,x) applies the operator identity + dt *
% SLP*Bending where \kappa is absored into the term Bending

N = vesicle.N; nv = vesicle.nv;
xCols = zeros(2*N,nv);
for k = 1:nv
  xCols(:,k) = x(2*(k-1)*N+1:k*2*N); 
end

xCols = -vesicle.bendingTerm(xCols);
% compute Bending * x
for k = 1:nv
  xCols(:,k) = o.Galpert(:,:,k)*xCols(:,k);
end
x = x + o.dt*xCols(:);

end % IplusSBMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = precoIplusSB(o,vesicle,alpha,scaling,x)
% val = precoIplusSB(vesicle,alpha,scaling,x) solves the equation (I +
% dt*kappa*SLP*Ben)*val = x analytically for a circle

N = vesicle.N; nv = vesicle.nv;
val = zeros(2*N,nv);

for k = 1:nv
  val(:,k) = x((k-1)*2*N+1:k*2*N);
end
val(1:N,:) = fft(val(1:N,:));
val(N+1:2*N,:) = fft(val(N+1:2*N,:));

for k = 1:nv
  val(3:N-1,k) = val(3:N-1,k).*scaling(3:N-1);
  val(N+3:2*N-1,k) = val(N+3:2*N-1,k).*scaling(3:N-1);
end
% all modes with modulus greater than 1 are diagonal

plusOneMode = [val(2,:) ; val(N+2,:)];
minusOneMode = [val(N,:) ; val(2*N,:)];
% need to save these values since the plus and minus one modes
% communicate ie. isn't diagonal

val(N,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  ((2*alpha^2-4*alpha+1)*val(N,:) - ...
  2*1i*alpha^2*val(2*N,:)) + ...
  alpha/(4*alpha-1) * ...
  (plusOneMode(1,:) + 1i*plusOneMode(2,:));
% -1 mode of the x component

val(2*N,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  (2*1i*alpha^2*val(N,:) + ...
  (2*alpha^2-4*alpha+1)*val(2*N,:)) + ...
  alpha/(4*alpha-1) * ...
  (1i*plusOneMode(1,:) - plusOneMode(2,:));
% -1 mode of the y component

val(2,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  ((2*alpha^2-4*alpha+1)*val(2,:) + ...
  2*1i*alpha^2*val(N+2,:)) + ...
  alpha/(4*alpha-1) * ...
  (minusOneMode(1,:) - 1i*minusOneMode(2,:));
% 1 mode of the x component

val(N+2,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  (-2*1i*alpha^2*val(2,:) + ...
  (2*alpha^2-4*alpha+1)*val(N+2,:)) + ...
  alpha/(4*alpha-1) * ...
  (-1i*minusOneMode(1,:) - minusOneMode(2,:));
% 1 mode of the y component

val(1:N,:) = real(ifft(val(1:N,:)));
val(N+1:2*N,:) = real(ifft(val(N+1:end,:)));
val = val(:);

end % precoIplusSB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time,pos] = noVesTracers(o,walls,wallsInt,wallsExt,Xtra)
% [time,pos] = noVesTracers(walls,wallsInt,wallsExt,Xtra) 
% computes the position of the
% tracers based only on the density function defined on the solid walls
fprintf('\n Computing the tracers without vesicles \n')
logFile = [o.runName '_tracers.log'];
fid = fopen(logFile,'a');
fprintf(fid,'Computing the tracers without vesicles\n');
fclose(fid);

if ~o.diffDiscWalls
  [eta,~,~,RS] = walls.computeEta([],o);
else
  [~,etaInt,etaExt,RS] = wallsInt.computeEta(wallsExt,o);  
end

% build a structure for tracers which mimics the class capsules.
% Can't use capsules unless there is an even number of tracer
% points (because of FFTs), so we don't want to build the class
% as we do for vesicle and walls
tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;

jump = -1/2;
if ~o.diffDiscWalls
  DLP = o.wallDLP;
  for k = 1:walls.nv
    DLP(:,:,k) = DLP(:,:,k) + jump*eye(2*walls.N);
  end
else
  DLPext = o.wallDLPext + jump*eye(2*wallsExt.N);
  DLPint = o.wallDLPint;
  for k = 1:wallsInt.nv
    DLPint(:,:,k) = DLPint(:,:,k) + jump*eye(2*wallsInt.N);    
  end
end

tfinal = 200;
nTsteps = tfinal/o.deltaTtracers+1;

if ~o.diffDiscWalls
  odefun = @(t,z) o.tracersVel4RK(t,tfinal,z,walls,[],[],...
      tracers,DLP,[],[],eta,[],[],RS,logFile);
else
  odefun = @(t,z) o.tracersVel4RK(t,tfinal,z,[],wallsInt,wallsExt,...
      tracers,[],DLPint,DLPext,[],etaInt,etaExt,RS,logFile);
end
options.RelTol = 1e-8;
options.AbsTol = 1e-10;
[time,pos] = ode15s(odefun,linspace(0,tfinal,nTsteps),Xtra,options);

end % noVesTracers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = tracersVel4RK(o,t,T,Xtra,walls,wallsInt,wallsExt,...
      tracers,DLP,DLPint,DLPext,eta,etaInt,etaExt,RS,logFile)  
%vel = tracersVel4RK(o,t,T,Xtra,walls,wallsInt,wallsExt,...
% tracers,DLP,DLPint,DLPext,eta,etaInt,etaExt,RS,logFile) evaluates the
% right-hand-side for doing time stepping of the tracers when using RK15

fid = fopen(logFile,'a');
message = ['ode15s ' num2str(t/T*100,'%04.1f') ' %% completed'];
nmess = numel(message);
fprintf(repmat('\b',1,nmess));
fprintf(message);
fprintf(fid,[message '\n']);
fclose(fid);

% save TransientTracers t Xtra walls 

if ~o.diffDiscWall
  op = o.opWall;

  tracers.X = Xtra;
  [~,NearW2T] = walls.getZone(tracers,2);
  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
    kernelDirect = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLnewfmm;
    kernelDirect = @op.exactStokesDL;
  end
  
  DLPfun = @(X) op.exactStokesDLdiag(walls,DLP,X);
  Fwall2Tra = op.nearSingInt(walls,eta,DLPfun,[],...
      NearW2T,kernel,kernelDirect,tracers,false,false);
 
  FLets2Tra = zeros(size(Xtra));
  for k = 2:walls.nv
    stokeslet = RS(1:2,k);
    rotlet = RS(3,k);
    FLets2Tra = FLets2Tra + o.RSlets(Xtra,walls.center(:,k),...
      stokeslet,rotlet);
  end

  vel = Fwall2Tra + FLets2Tra;
  indx = find(Xtra(1:end/2) >= max(walls.X(1:end/2,1))*0.9);
  indy = find(abs(Xtra(end/2+1:end)) >= ...
      max(abs(walls.X(end/2+1:end,1)))*0.95);
  vel(indx) = 0;
  vel(indx + tracers.N) = 0;

  vel(indy) = 0;
  vel(indy + tracers.N) = 0;
else
  opInt = o.opWallInt;
  opExt = o.opWallExt;

  tracers.X = Xtra;
  [~,NearWint2T] = wallsInt.getZone(tracers,2);
  [~,NearWext2T] = wallsExt.getZone(tracers,2);

  if ~o.fmmDLP
    kernelInt = @opInt.exactStokesDL;
    kernelExt = @opExt.exactStokesDL;
    kernelDirectInt = @opInt.exactStokesDL;
    kernelDirectExt = @opExt.exactStokesDL;
  else
    kernelInt = @opInt.exactStokesDLnewfmm;
    kernelExt = @opExt.exactStokesDLnewfmm;
    kernelDirectInt = @opInt.exactStokesDL;
    kernelDirectExt = @opExt.exactStokesDL;
  end

  DLPfun = @(X) opInt.exactStokesDLdiag(wallsInt,DLPint,X);
  FwallInt2Tra = opInt.nearSingInt(wallsInt,etaInt,DLPfun,[],...
   NearWint2T,kernelInt,kernelDirectInt,tracers,false,false);

  DLPfun = @(X) opExt.exactStokesDLdiag(wallsExt,DLPext,X);
  FwallExt2Tra = opExt.nearSingInt(wallsExt,etaExt,DLPfun,[],...
   NearWext2T,kernelExt,kernelDirectExt,tracers,false,false);

  Fwall2Tra = FwallInt2Tra + FwallExt2Tra;

  FLets2Tra = zeros(size(Xtra));
  for k = 2:wallsInt.nv+1
    stokeslet = RS(1:2,k);
    rotlet = RS(3,k);
    FLets2Tra = FLets2Tra + o.RSlets(Xtra,wallsInt.center(:,k-1),...
      stokeslet,rotlet);
  end

  vel = Fwall2Tra + FLets2Tra;
  indx = find(Xtra(1:end/2) >= max(wallsExt.X(1:end/2))*0.9);
  indy = find(abs(Xtra(end/2+1:end)) >= ...
      max(abs(wallsExt.X(end/2+1:end)))*0.95);
  vel(indx) = 0;
  vel(indx + tracers.N) = 0;

  vel(indy) = 0;
  vel(indy + tracers.N) = 0;  
    
end % ~o.diffDiscWalls

end % tracersVel4RK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vel,allIdx,rmIdx] = computeVelFieldNoVes(o,prams,...
        Xwalls,XwallsInt,XwallsExt,XtraAll,iRemovePoints)

% build a structure for tracers which mimics the class capsules.
% Can't use capsules unless there is an even number of tracer
% points (because of FFTs), so we don't want to build the class
% as we do for vesicle and walls
tracers.N = numel(XtraAll)/2;
tracers.nv = 1;
tracers.X = XtraAll;

if o.diffDiscWalls
  % Build tracers structure for the remaining points
  xxtra = XtraAll(1:end/2); yytra = XtraAll(end/2+1:end);  
  if iRemovePoints % if there might be points crossing the walls, then remove them
    % Enlarge pillars so that we do not include the points around pillars 
    XwallsIntLarge = zeros(size(XwallsInt));
    for k = 1 : size(XwallsInt,2)
      xpill = XwallsInt(1:end/2,k); ypill = XwallsInt(end/2+1:end,k);
      centxPill = mean(xpill); centyPill = mean(ypill);
      xpilln = 1.05*(xpill-centxPill)+centxPill;
      ypilln = 1.05*(ypill-centyPill)+centyPill;
      XwallsIntLarge(:,k) = [xpilln;ypilln];
    end
    % Also shrink the exterior wall
    centxWall = mean(XwallsExt(1:end/2)); 
    centyWall = mean(XwallsExt(end/2+1:end));
    xExt = 0.99*(XwallsExt(1:end/2)-centxWall)+centxWall;
    yExt = 0.99*(XwallsExt(end/2+1:end)-centyWall)+centyWall;
    XwallsExtSmall = [xExt;yExt];
  
    % We do not need any of the wall matrices that are built in
    % initialConfined to sort points. So if they are already built, let
    % them stay as they are.
    wallsInt = capsules(XwallsIntLarge,[],[],zeros(prams.nvbdInt,1),...
        zeros(prams.nvbdInt,1),o.antiAlias);
    wallsExt = capsules(XwallsExtSmall,[],[],0,0,o.antiAlias);

    % Get the near structure
    [~,NearWext2T] = wallsExt.getZone(tracers,2);
    [~,NearWint2T] = wallsInt.getZone(tracers,2);
  
    % indices of all points
    allIdx = 1:numel(XtraAll)/2;
  
    % Remove the points outside the exterior wall
    inOut = wallsExt.sortPts(XtraAll,o.fmm,NearWext2T,o.opWallExt);
    rmIdxOut1 = find(inOut==0);
  
    % Remove the points inside the pillars
    inOut = wallsInt.sortPts(XtraAll,o.fmm,NearWint2T,o.opWallInt);
    rmIdxOut2 = find(inOut~=0);

    % keep the indices of removed points
    rmIdx = [rmIdxOut1;rmIdxOut2];

    % Find the indices of remaining points
    for idx = 1 : numel(rmIdx)
      jd = find(allIdx~=rmIdx(idx));
      if ~isempty(jd)
        allIdx = allIdx(jd);
      end
    end
    Xtra = [xxtra(allIdx);yytra(allIdx)];
    % Since we built wall matrices with enlarged walls, we need to remove
    % those matrices as they are not correct or replace them with the
    % correct ones if they were already built.
    
  else
    % we know all the points are in the domain and no crossing.
    
    Xtra = XtraAll;
    allIdx = 1:numel(XtraAll)/2;
    rmIdx = [];
  end
  tracers.N = numel(Xtra)/2;
  tracers.nv = 1;
  tracers.X = Xtra;
  Ntra = tracers.N;

  % Build walls structures using the correct pillar shapes (not enlarged)
  [~,wallsInt,wallsExt] = o.initialConfined(prams,[],...
      XwallsInt,XwallsExt); 
  nvbd = wallsInt.nv+1;
  % Get the near structure for them, too
  [~,NearWint2T] = wallsInt.getZone(tracers,2);
  [~,NearWext2T] = wallsExt.getZone(tracers,2);

  % compute density
  [~,etaInt,etaExt,RS] = wallsInt.computeEta(wallsExt,o);
  
  % compute the velocity on the tracers due to the
  % solid walls
  jump = -1/2;

  potWallInt = o.opWallInt; 
  if ~o.fmmDLP
    kernel = @potWallInt.exactStokesDL;
    kernelDirect = @potWallInt.exactStokesDL;
  else
    kernel = @potWallInt.exactStokesDLnewfmm;
    kernelDirect = @potWallInt.exactStokesDL;
  end

  DLP = @(X) jump*X + potWallInt.exactStokesDLdiag(wallsInt,...
    o.wallDLPint,X);
  FwallInt2Tra = potWallInt.nearSingInt(wallsInt,etaInt,DLP,[],...
    NearWint2T,kernel,kernelDirect,tracers,false,false);

  potWallExt = o.opWallExt;
  if ~o.fmmDLP
    kernel = @potWallExt.exactStokesDL;
    kernelDirect = @potWallExt.exactStokesDL;
  else
    kernel = @potWallExt.exactStokesDLnewfmm;
    kernelDirect = @potWallExt.exactStokesDL;
  end

  DLP = @(X) jump*X + potWallExt.exactStokesDLdiag(wallsExt,...
    o.wallDLPext,X);
  FwallExt2Tra = potWallExt.nearSingInt(wallsExt,etaExt,DLP,[],...
    NearWext2T,kernel,kernelDirect,tracers,false,false);
  
  
  Fwall2Tra = FwallInt2Tra + FwallExt2Tra;

  % compute velocity due to rotlets and stokeslets on the vesicles
  FLets2Tra = zeros(2*Ntra,1);
  for k = 2:nvbd
    stokeslet = RS(1:2,k);
    rotlet = RS(3,k);
    FLets2Tra = FLets2Tra + o.RSlets(tracers.X,wallsInt.center(:,k-1),...
        stokeslet,rotlet);
  end

  vel = Fwall2Tra + FLets2Tra;
  % END OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY
else
  [walls,~,~] = o.initialConfined(prams,Xwalls,[],[]); 
  nvbd = walls.nv+1; % Number of components to walls

  % HERE, WE SUPPOSE ALL THE TRACER POINTS ARE IN DOMAIN

  Xtra = XtraAll;
  tracers.N = numel(Xtra)/2;
  tracers.nv = 1;
  tracers.X = Xtra;
  Ntra = tracers.N;

  % Get the near field structure
  [~,NearW2T] = walls.getZone(tracers,2);
  
  % Get the density
  [eta,~,~,RS] = walls.computeEta([],o);

  % compute the velocity on the tracers due to the
  % solid walls
  jump = -1/2;

  potWall = o.opWall; 
  if ~o.fmmDLP
    kernel = @potWall.exactStokesDL;
    kernelDirect = @potWall.exactStokesDL;
  else
    kernel = @potWall.exactStokesDLnewfmm;
    kernelDirect = @potWall.exactStokesDL;
  end

  DLP = @(X) jump*X + potWall.exactStokesDLdiag(walls,o.wallDLP,X);
  Fwall2Tra = potWall.nearSingInt(walls,eta,DLP,[],...
    NearW2T,kernel,kernelDirect,tracers,false,false);

  % compute velocity due to rotlets and stokeslets on the vesicles
  FLets2Tra = zeros(2*Ntra,1);
  for k = 2:nvbd
    stokeslet = RS(1:2,k);
    rotlet = RS(3,k);
    FLets2Tra = FLets2Tra + o.RSlets(tracers.X,walls.center(:,k),...
        stokeslet,rotlet);
  end

  vel = Fwall2Tra + FLets2Tra;
  % END OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY  
    
end % o.diffDiscWalls

end % computeVelFieldNoVes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vel,allIdx,rmIdx] = computeVelFieldMatFree(o,prams,XwallsInt,...
        XwallsExt,XtraAll,scaling,iRemovePoints)

% build a structure for tracers which mimics the class capsules.
% Can't use capsules unless there is an even number of tracer
% points (because of FFTs), so we don't want to build the class
% as we do for vesicle and walls
tracers.N = numel(XtraAll)/2;
tracers.nv = 1;
tracers.X = XtraAll;

NbdInt = prams.NbdInt;
NbdExt = prams.NbdExt;
nvbdInt = prams.nvbdInt;
nvbd = nvbdInt+1;

% Build tracers structure for the remaining points
xxtra = XtraAll(1:end/2); yytra = XtraAll(end/2+1:end);  
if iRemovePoints % if there might be points crossing the walls, then remove them
  % Enlarge pillars so that we do not include the points around pillars 
  XwallsIntLarge = zeros(size(XwallsInt));
  for k = 1 : size(XwallsInt,2)
    xpill = XwallsInt(1:end/2,k); ypill = XwallsInt(end/2+1:end,k);
    centxPill = mean(xpill); centyPill = mean(ypill);
    xpilln = 1.05*(xpill-centxPill)+centxPill;
    ypilln = 1.05*(ypill-centyPill)+centyPill;
    XwallsIntLarge(:,k) = [xpilln;ypilln];
  end
  % Also shrink the exterior wall
  centxWall = mean(XwallsExt(1:end/2)); 
  centyWall = mean(XwallsExt(end/2+1:end));
  xExt = 0.99*(XwallsExt(1:end/2)-centxWall)+centxWall;
  yExt = 0.99*(XwallsExt(end/2+1:end)-centyWall)+centyWall;
  XwallsExtSmall = [xExt;yExt];
  
  % We do not need any of the wall matrices that are built in
  % initialConfined to sort points. So if they are already built, let
  % them stay as they are.
  wallsInt = capsules(XwallsIntLarge,[],[],zeros(prams.nvbdInt,1),...
      zeros(prams.nvbdInt,1),o.antiAlias);
  wallsExt = capsules(XwallsExtSmall,[],[],0,0,o.antiAlias);
  % Get the near structure
  [~,NearWext2T] = wallsExt.getZone(tracers,2);
  [~,NearWint2T] = wallsInt.getZone(tracers,2);
  
  % indices of all points
  allIdx = 1:numel(XtraAll)/2;
  
  % Remove the points outside the exterior wall
  inOut = wallsExt.sortPts(XtraAll,o.fmm,NearWext2T,o.opWallExt);
  rmIdxOut1 = find(inOut==0);
  
  % Remove the points inside the pillars
  inOut = wallsInt.sortPts(XtraAll,o.fmm,NearWint2T,o.opWallInt);
  rmIdxOut2 = find(inOut~=0);

  % keep the indices of removed points
  rmIdx = [rmIdxOut1;rmIdxOut2];

  % Find the indices of remaining points
  for idx = 1 : numel(rmIdx)
    jd = find(allIdx~=rmIdx(idx));
    if ~isempty(jd)
      allIdx = allIdx(jd);
    end
  end
  Xtra = [xxtra(allIdx);yytra(allIdx)];
  % Since we built wall matrices with enlarged walls, we need to remove
  % those matrices as they are not correct or replace them with the
  % correct ones if they were already built.
else
  % we know all the points are in the domain and no crossing.
    
  Xtra = XtraAll;
  allIdx = 1:numel(XtraAll)/2;
  rmIdx = [];
end
  tracers.N = numel(Xtra)/2;
  tracers.nv = 1;
  tracers.X = Xtra;
  Ntra = tracers.N;

  % BUILD WALLS AND MATRICES using the correct pillar shapes (not enlarged)
  % but do not build the W2W interaction matrix
  potWallInt = o.opWallInt;
  potWallExt = o.opWallExt;
  
  % velocity on the solid walls from the BCs
  [uwallsExt,uwallsInt] = o.farField(XwallsExt,XwallsInt);
  % build the walls
  wallsInt = capsules(XwallsInt,[],uwallsInt,...
      zeros(prams.nvbdInt,1),zeros(prams.nvbdInt,1),o.antiAlias);
  wallsExt = capsules(XwallsExt,[],uwallsExt,0,0,o.antiAlias);
  % set the upsampling rate for walls
  if o.antiAlias
    wallsInt.setUpRate(potWallInt);
    wallsExt.setUpRate(potWallExt);
  end
    
  % build the double layer potential matrix for walls and save on memory
  if isempty(o.wallDLPint)
    o.wallDLPint = potWallInt.stokesDLmatrix(wallsInt);
  end
  if isempty(o.wallDLPext)
    o.wallDLPext = potWallExt.stokesDLmatrix(wallsExt);
  end
  
  % build the DLP for walls without any correction if FMM is on
  if o.fmmDLP 
    if isempty(o.wallDLPintNoCorr)
      o.wallDLPintNoCorr = potWallInt.stokesDLmatrixNoCorr(wallsInt);
    end
    if isempty(o.wallDLPextNoCorr)
      o.wallDLPextNoCorr = potWallExt.stokesDLmatrixNoCorr(wallsExt);    
    end
  else
    o.wallDLPintNoCorr = [];
    o.wallDLPextNoCorr = [];
  end
  
  % N0 is computed only for the exterior wall
  if isempty(o.wallN0)
    o.wallN0 = potWallExt.stokesN0matrix(wallsExt);
  end
  
  % NOW, COMPUTE DENSITY USING GMRES 
  
  % Get the near structure for them, too
  [~,NearWint2T] = wallsInt.getZone(tracers,2);
  [~,NearWext2T] = wallsExt.getZone(tracers,2);

  % BUILD THE PRECONDITIONER
  useWallPreco = true;
  if useWallPreco
    if isempty(o.invM11Ext)
      disp('Constructing preconditioner')
      o.wallsBDprecondConstruct(wallsInt,wallsExt);    
    end
  end
  
  if isempty(scaling) && isempty(o.etaExt) % then we want to compute density
    % compute density using GMRES
    rhs = zeros(2*NbdInt*nvbdInt+2*NbdExt+3*(nvbd-1),1);
    rhs(1:2*NbdExt) = uwallsExt;
    for k = 1 : nvbdInt
      istart = 2*NbdExt + 1 + (k-1)*2*NbdInt;
      iend = istart + 2*NbdInt-1;
      rhs(istart:iend) = uwallsInt(:,k);
    end
    % GMRES
    maxGMRESiter = (2*NbdInt*nvbdInt+2*NbdExt+3*(nvbd-1));
    if useWallPreco
    [etaRS,iflag,relRes,I,~] = gmres(@(X) o.computeEtaMatVec(X,wallsInt,wallsExt),...
        rhs,[],1e-5,maxGMRESiter,@o.wallsBDprecondApply);
    else
    [etaRS,iflag,relRes,I,~] = gmres(@(X) o.computeEtaMatVec(X,wallsInt,wallsExt),...
        rhs,[],1e-5,maxGMRESiter);
    end
  

    disp(['GMRES took ' num2str(I(2)) ' iterations'])
    % Unstack
    etaInt = zeros(2*NbdInt,nvbdInt);
    RS = zeros(3,nvbd); % rotlets and stokeslets
    etaExt = etaRS(1:2*NbdExt); % assuming one outer wall
    for k = 1:nvbdInt
      istart = 2*NbdExt+(k-1)*2*NbdInt+1;
      iend = istart+2*NbdInt-1;
      etaInt(:,k) = etaRS(istart:iend);
    end
    sizeSoFar = 2*nvbdInt*NbdInt+2*NbdExt;
    
    for k = 2:nvbd
      istart = sizeSoFar + 3*(k-2) + 1;
      iend = istart + 2;
      RS(:,k) = etaRS(istart:iend);
    end
    % store them, we will need it if we scale the right hand side
    o.etaExt = etaExt;
    o.etaInt = etaInt;
    o.RS = RS;
  elseif ~isempty(scaling) && ~isempty(o.etaExt)
    % scaling Ufar just scales eta so do not compute again  
    etaExt = scaling*o.etaExt; o.etaExt = etaExt;
    etaInt = scaling*o.etaInt; o.etaInt = etaInt;
    RS = scaling*o.RS; o.RS = RS;
  elseif isempty(scaling) && ~isempty(o.etaExt)
    etaExt = o.etaExt;
    etaInt = o.etaInt;
    RS = o.RS;
  end
  
  idebug = 0;
  if idebug
    % compute density and RS exactly
    if isempty(o.bdiagWall)  
      o.bdiagWall = o.wallsPrecondDiffDisc(wallsInt,wallsExt);
    end
    [~,exactEtaInt,exactEtaExt,exactRS] = wallsInt.computeEta(wallsExt,o);
    errorInetaExt = norm(exactEtaExt-etaExt)/norm(exactEtaExt)
    errorInetaInt = norm(exactEtaInt-etaInt)/norm(exactEtaInt)
    errorInRS = norm(exactRS-RS)/norm(exactRS)
    pause
  end

  % compute the velocity on the tracers due to the
  % solid walls
  jump = -1/2;

  potWallInt = o.opWallInt; 
  if ~o.fmmDLP
    kernel = @potWallInt.exactStokesDL;
    kernelDirect = @potWallInt.exactStokesDL;
  else
    kernel = @potWallInt.exactStokesDLnewfmm;
    kernelDirect = @potWallInt.exactStokesDL;
  end

  DLP = @(X) jump*X + potWallInt.exactStokesDLdiag(wallsInt,...
    o.wallDLPint,X);
  FwallInt2Tra = potWallInt.nearSingInt(wallsInt,etaInt,DLP,[],...
    NearWint2T,kernel,kernelDirect,tracers,false,false);

  potWallExt = o.opWallExt;
  if ~o.fmmDLP
    kernel = @potWallExt.exactStokesDL;
    kernelDirect = @potWallExt.exactStokesDL;
  else
    kernel = @potWallExt.exactStokesDLnewfmm;
    kernelDirect = @potWallExt.exactStokesDL;
  end

  DLP = @(X) jump*X + potWallExt.exactStokesDLdiag(wallsExt,...
    o.wallDLPext,X);
  FwallExt2Tra = potWallExt.nearSingInt(wallsExt,etaExt,DLP,[],...
    NearWext2T,kernel,kernelDirect,tracers,false,false);
  
  
  Fwall2Tra = FwallInt2Tra + FwallExt2Tra;

  % compute velocity due to rotlets and stokeslets on the vesicles
  FLets2Tra = zeros(2*Ntra,1);
  for k = 2:nvbd
    stokeslet = RS(1:2,k);
    rotlet = RS(3,k);
    FLets2Tra = FLets2Tra + o.RSlets(tracers.X,wallsInt.center(:,k-1),...
        stokeslet,rotlet);
  end

  vel = Fwall2Tra + FLets2Tra;
  % END OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY


end % computeVelFieldNoVes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xn = computeEtaMatVec(o,X,wallsInt,wallsExt)

NbdInt = wallsInt.N;
NbdExt = wallsExt.N;
nvbdInt = wallsInt.nv;
nvbdExt = 1;
nvbd = nvbdInt + 1;
    
% Unstack X
etaInt = zeros(2*NbdInt,nvbdInt);
RS = zeros(3,nvbd);
etaExt = X(1:2*NbdExt); % assuming one outer wall
for k = 1:nvbdInt
  istart = 2*NbdExt+(k-1)*2*NbdInt+1;
  iend = istart+2*NbdInt-1;
  etaInt(:,k) = X(istart:iend);
end
sizeSoFar = 2*nvbdInt*NbdInt+2*nvbdExt*NbdExt;
RS = X(2*NbdExt+2*nvbdInt*NbdInt+1:end);

% EVALUATE WALL-TO-WALL INTERACTIONS
potWallInt = o.opWallInt;
potWallExt = o.opWallExt;
if ~o.fmmDLP
  kernel = @potWallExt.exactStokesDL;
  [~,FDLPwallExt2wallInt] = kernel(wallsExt,etaExt,[],wallsInt.X,1);

  kernel = @potWallInt.exactStokesDL;
  FDLPwallInt2wallInt = kernel(wallsInt,etaInt,[]);
  [~,FDLPwallInt2wallExt] = kernel(wallsInt,etaInt,[],...
      wallsExt.X,1:wallsInt.nv);
else
  kernel = @potWallInt.exactStokesDLnewfmm;
  FDLPwallInt2wallInt = kernel(wallsInt,etaInt,o.wallDLPintNoCorr);

  [~,FDLPwallExt2wallInt] = kernel(wallsExt,etaExt,[],wallsInt.X,1);
  [~,FDLPwallInt2wallExt] = kernel(wallsInt,...
      etaInt,[],wallsExt.X,1:wallsInt.nv);
end % o.fmmDLP

% EVALUATE POTENTIAL DUE TO STOKESLETS AND ROTLETS
LetsWallsInt = zeros(2*NbdInt,nvbdInt);
LetsWallsExt = zeros(2*NbdExt,nvbdExt);
for k = 2:nvbd
  stokeslet = RS(3*(k-2)+1:3*(k-2)+2);
  rotlet = RS(3*(k-1));

  % compute velocity due to rotlets and stokeslets on the solid walls
  LetsWallsInt = LetsWallsInt + o.RSlets(wallsInt.X,...
      wallsInt.center(:,k-1),stokeslet,rotlet);
  LetsWallsExt = LetsWallsExt + o.RSlets(wallsExt.X,...
      wallsInt.center(:,k-1),stokeslet,rotlet);
end
% Integral constraints on the density function eta related
% to the weights of the stokeslets and rotlets
valLets = o.letsIntegrals(RS,[],etaInt,[],wallsInt);

% EVALUATE VELOCITY ON WALLS
valWallsInt = zeros(2*NbdInt,nvbdInt);
valWallsExt = zeros(2*NbdExt,nvbdExt);

valWallsInt = valWallsInt - 1/2*etaInt + ...
  potWallInt.exactStokesDLdiag(wallsInt,o.wallDLPint,etaInt);

valWallsExt = valWallsExt - 1/2*etaExt + ...
  potWallExt.exactStokesDLdiag(wallsExt,o.wallDLPext,etaExt);
valWallsExt(:,1) = valWallsExt(:,1) + ...
  potWallExt.exactStokesN0diag(wallsExt,o.wallN0,etaExt); 

% velocity on walls due to all other walls
valWallsInt = valWallsInt + FDLPwallExt2wallInt;

% velocity on walls due to all other walls
valWallsInt = valWallsInt + FDLPwallInt2wallInt;

% velocity on walls due to the rotlets and stokeslets  
valWallsInt = valWallsInt + LetsWallsInt;

% velocity on walls due to all other walls
valWallsExt = valWallsExt + FDLPwallInt2wallExt;

% velocity on walls due to the rotlets and stokeslets
valWallsExt = valWallsExt + LetsWallsExt;

% pad everythin
Xn = [valWallsExt(:);valWallsInt(:);valLets];

end % computeEtaMatVec
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = tracersVel(o,X,sigma,u,kappa,viscCont,...
    walls,wallsInt,wallsExt,eta,etaInt,etaExt,RS,Xtra)
% vel = tracersVel(o,X,sigma,u,kappa,viscCont,...
%    walls,wallsInt,wallsExt,eta,etaInt,etaExt,RS,Xtra) computes the
% velocity at the set of points Xtra due to the vesicle, viscosity
% contrast, solid walls, and stokeslets and rotlets.  The particles Xtra
% are treated passively

% build a structure for tracers which mimics the class capsules.
% Can't use capsules unless there is an even number of tracer
% points (because of FFTs), so we don't want to build the class
% as we do for vesicle and walls
vesicle = capsules(X,sigma,[],kappa,viscCont,o.antiAlias);
tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;

% build near-singular integration structures for vesicle to tracer and
% wall to tracer interactions
[~,NearV2T] = vesicle.getZone(tracers,2);

if o.confined
if~o.diffDiscWalls
  % Need wall-to-tracers interactions
  [~,NearW2T] = walls.getZone(tracers,2);


  NearWint2T = [];
  NearWext2T = [];
else
  % When walls are discretized separately, we cannot call getZone 
  % for all the walls. That function is not adjusted for that. Instead
  % this if statement is added and near structures are formed
  % separately.
  [~,NearWint2T] = wallsInt.getZone(tracers,2);
  [~,NearWext2T] = wallsExt.getZone(tracers,2);
  NearW2T = [];
end
else
NearW2T = []; NearWint2T = []; NearWext2T = [];
end

Ntra = size(Xtra,1)/2; % Number of tracers
N = size(X,1)/2; % Number of points
nv = size(X,2); % Number of vesicles
if o.confined
  if ~o.diffDiscWalls
    Nbd = walls.N; % Number of points on walls
    NbdInt = 0; NbdExt = 0;
    nvbd = walls.nv; % Number of components to walls  
  else
    NbdInt = wallsInt.N; % Number of points on interior walls
    NbdExt = wallsExt.N; % Number of points on exterior walls
    nvbd = wallsInt.nv + wallsExt.nv; % Number of components to walls  
    Nbd = 0;
  end
else
  Nbd = 0; NbdInt = 0; NbdExt = 0;
  nvbd = 0;
end

% compute traction jump
f = vesicle.tracJump(X,sigma);

% poten classes
op = poten(N,o.fmmPrecision,1);
if o.confined
  if ~o.diffDiscWalls
    opWall = poten(Nbd,o.fmmPrecision,1);
  else
    opWallInt = poten(NbdInt,o.fmmPrecision,1);
    opWallExt = poten(NbdExt,o.fmmPrecision,1);
  end
end

%% START OF VELOCITY DUE TO VESICLES
if ~o.fmm
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
  kernelDirect = @op.exactStokesSL;
end

% need single-layer potential matrix for doing near-singular
% integration
o.Galpert = op.stokesSLmatrix(vesicle);
SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);

% Evaulate velocity due to vesicle on tracers with near-singular
% integration
Fves2Tra = op.nearSingInt(vesicle,f,SLP,[],...
   NearV2T,kernel,kernelDirect,tracers,false,false);
% END OF VELOCITY DUE TO VESICLES

% START OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY
if o.confined
  if ~o.diffDiscWalls  
    if ~o.fmmDLP
      kernel = @opWall.exactStokesDL;
      kernelDirect = @opWall.exactStokesDL;
    else
      kernel = @opWall.exactStokesDLnewfmm;
      kernelDirect = @opWall.exactStokesDL;
    end
    % compute the velocity on the tracers due to the
    % solid walls
    DLP = o.wallDLP;
    jump = -1/2;
    for k = 1:nvbd
      DLP(:,:,k) = DLP(:,:,k) + jump*eye(2*Nbd);
    end
    DLPfun = @(X) opWall.exactStokesDLdiag(walls,DLP,X);
    Fwall2Tra = opWall.nearSingInt(walls,eta,DLPfun,[],...
        NearW2T,kernel,kernelDirect,tracers,false,false);

    % velocity due to the stokeslets and rotlets
    FLets2Tra = zeros(2*Ntra,1);
    for k = 2:nvbd
      stokeslet = RS(1:2,k);
      rotlet = RS(3,k);
      FLets2Tra = FLets2Tra + o.RSlets(Xtra,walls.center(:,k),...
        stokeslet,rotlet);
    end
  else
    if ~o.fmmDLP
     kernelInt = @opWallInt.exactStokesDL;
     kernelExt = @opWallExt.exactStokesDL;
     kernelDirectInt = @opWallInt.exactStokesDL;
     kernelDirectExt = @opWallExt.exactStokesDL;
    else
     kernelInt = @opWallInt.exactStokesDLnewfmm;
     kernelExt = @opWallExt.exactStokesDLnewfmm;
     kernelDirectInt = @opWallInt.exactStokesDL;
     kernelDirectExt = @opWallExt.exactStokesDL;
    end

    jump = -1/2;
    DLP = o.wallDLPext + jump*eye(2*NbdExt);
    DLPfun = @(X) opWallExt.exactStokesDLdiag(wallsExt,DLP,X);
    FwallExt2Tra = opWallExt.nearSingInt(wallsExt,etaExt,DLPfun,[],...
        NearWext2T,kernelExt,kernelDirectExt,tracers,false,false);

    DLP = o.wallDLPint;
    for k = 1:nvbd-1
      DLP(:,:,k) = DLP(:,:,k) + jump*eye(2*NbdInt);
    end
    DLPfun = @(X) opWallInt.exactStokesDLdiag(wallsInt,DLP,X);
    FwallInt2Tra = opWallInt.nearSingInt(wallsInt,etaInt,DLPfun,[],...
        NearWint2T,kernelInt,kernelDirectInt,tracers,false,false);

    % compute the velocity on the tracers due to the
    % solid walls
    Fwall2Tra = FwallExt2Tra + FwallInt2Tra;
    
    % velocity due to the stokeslets and rotlets    
    FLets2Tra = zeros(2*Ntra,1);
    for k = 2:nvbd
      stokeslet = RS(1:2,k);
      rotlet = RS(3,k);
      FLets2Tra = FLets2Tra + o.RSlets(Xtra,wallsInt.center(:,k-1),...
        stokeslet,rotlet);
    end  
  end % if ~o.diffDiscWalls
  
else
  % if no solid walls, velocity field comes from background
  % velocity.  No stokeslet or rotlets so they induce no
  % velocity  
  Fwall2Tra = o.farField(Xtra,[]);
  FLets2Tra = zeros(2*Ntra,1);
end
% END OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY

% START OF VELOCITY DUE TO VISCOSITY CONTRAST
if any(viscCont ~= 1)
  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
    kernelDirect = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLnewfmm;
    kernelDirect = @op.exactStokesDL;
  end

  % need double-layer potential matrix for doing near-singular
  % integration
  o.D = op.stokesDLmatrix(vesicle);
  
  DLP = o.D;
  % Get the DLP for inside&outside
  DLP_out = DLP;
  DLP_in  = DLP;
  for k=1:nv
    % if tracer is outside a vesicle, jump:  
    jump_out = 1/2*(1-viscCont(k));
    % if tracer is inside a vesicle, jump:
    jump_in  = -jump_out;
      
    % Based on jump, DLP is computed for inside/outside
    DLP_out(:,:,k) = DLP_out(:,:,k) + jump_out*eye(2*N);
    DLP_in (:,:,k) = DLP_in (:,:,k) + jump_in *eye(2*N);
      
    DLPfun_out = @(X) op.exactStokesDLdiag(vesicle,DLP_out,X);
    DLPfun_in  = @(X) op.exactStokesDLdiag(vesicle,DLP_in ,X);
      
    % TODO: THIS JUMP NEEDS TO BE NEGATED FOR TRACERS INSIDE
    % THE VESICLE BECAUSE OF THE DISCONTINUITY OF THE DOUBLE-
    % LAYER POTENTIAL.
  end
  

  
  
  % Initialize Fvisc2Tra
  Fvisc2Tra = zeros(2*Ntra,1);
  % Determine if the points inside or outside
  InOut = vesicle.sortPts(Xtra,true);
  InOut_ext = [InOut;InOut]; 
  % since size(InOut) = size(Xtra)/2 (only for one component) 
    
  % OUTSIDE: build tracers and compute Fvisc2Tra
  %----------------------------------------------------------------------
  Xtra_out = Xtra(InOut_ext == 0);
  tracers.N = numel(Xtra_out)/2;
  tracers.X = Xtra_out;
  [~,NearV2T_out] = vesicle.getZone(tracers,2);
  
  % Evaulate velocity due to vesicle's viscosity contrast on 
  % tracers with near-singular integration
  Fvisc2Tra(InOut_ext==0) = op.nearSingInt(vesicle,u,DLPfun_out,[],...
    NearV2T_out,kernel,kernelDirect,tracers,false,false);
  
    
  % INSIDE: build tracers and compute Fvisc2Tra
  %----------------------------------------------------------------------
  Xtra_in = Xtra(InOut_ext == 1);
  tracers.N = numel(Xtra_in)/2;
  tracers.X = Xtra_in;
  [~,NearV2T_in] = vesicle.getZone(tracers,2);
    
  % Evaulate velocity due to vesicle's viscosity contrast on 
  % tracers with near-singular integration
  Fvisc2Tra(InOut_ext==1) = op.nearSingInt(vesicle,u,DLPfun_in,[],...
    NearV2T_in,kernel,kernelDirect,tracers,false,false);
    
  % Back to original
  tracers.N = numel(Xtra)/2;
  tracers.X = Xtra;
else
  % no velocity due to viscosity contrast  
  Fvisc2Tra = zeros(2*Ntra,1);
end
% END OF VELOCITY DUE TO VISCOSITY CONTRAST


% velocity is the sum of velocities due to vesicles, solid walls,
% rotlet and stokeslets, and viscosity contrast
vel = Fves2Tra + Fwall2Tra + FLets2Tra + Fvisc2Tra;


%vel = Fwall2Tra + FLets2Tra;
% DEBUG: FOR TESTING TRACERS WITH NO VESICLES
%r2 = Xtra(1).^2 + Xtra(2).^2;
%speed = -1/3 + 400/3/r2;
%velTrue = speed*[-Xtra(2);Xtra(1)];
% This is the true velocity for the couette apparatus when both
% boundaries are centered at the origin, the inner boundary rotates
% once every 2*pi time units and the outer boundary is fixed

end % tracersVel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = tracersVelDLD(o,X,sigma,u,kappa,viscCont,...
    wallsInt,wallsExt,etaInt,etaExt,RS,Xtra,NearWint2T,NearWext2T)
% vel = tracersVel(o,X,sigma,u,kappa,viscCont,...
%    walls,wallsInt,wallsExt,eta,etaInt,etaExt,RS,Xtra) computes the
% velocity at the set of points Xtra due to the vesicle, viscosity
% contrast, solid walls, and stokeslets and rotlets.  The particles Xtra
% are treated passively

% build a structure for tracers which mimics the class capsules.
% Can't use capsules unless there is an even number of tracer
% points (because of FFTs), so we don't want to build the class
% as we do for vesicle and walls
vesicle = capsules(X,sigma,[],kappa,viscCont,o.antiAlias);
vesicle.setUpRate();
tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;

% build near-singular integration structures for vesicle to tracer and
% wall to tracer interactions
[~,NearV2T] = vesicle.getZone(tracers,2);

Ntra = size(Xtra,1)/2; % Number of tracers
N = size(X,1)/2; % Number of points
nv = size(X,2); % Number of vesicles

NbdInt = wallsInt.N; % Number of points on interior walls
NbdExt = wallsExt.N; % Number of points on exterior walls
nvbd = wallsInt.nv + wallsExt.nv; % Number of components to walls  
Nbd = 0;


% compute traction jump
f = vesicle.tracJump(X,sigma);

% poten classes
op = o.op;
opWallInt = o.opWallInt;
opWallExt = o.opWallExt;

%% START OF VELOCITY DUE TO VESICLES
if ~o.fmm
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
  kernelDirect = @op.exactStokesSL;
end

% need single-layer potential matrix for doing near-singular
% integration
o.Galpert = op.stokesSLmatrix(vesicle);
SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);

% Evaulate velocity due to vesicle on tracers with near-singular
% integration
Fves2Tra = op.nearSingInt(vesicle,f,SLP,[],...
   NearV2T,kernel,kernelDirect,tracers,false,false);
% END OF VELOCITY DUE TO VESICLES

% START OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY
if ~o.fmmDLP
 kernelInt = @opWallInt.exactStokesDL;
 kernelExt = @opWallExt.exactStokesDL;
 kernelDirectInt = @opWallInt.exactStokesDL;
 kernelDirectExt = @opWallExt.exactStokesDL;
else
 kernelInt = @opWallInt.exactStokesDLnewfmm;
 kernelExt = @opWallExt.exactStokesDLnewfmm;
 kernelDirectInt = @opWallInt.exactStokesDL;
 kernelDirectExt = @opWallExt.exactStokesDL;
end

jump = -1/2;
DLP = o.wallDLPext + jump*eye(2*NbdExt);
DLPfun = @(X) opWallExt.exactStokesDLdiag(wallsExt,DLP,X);
FwallExt2Tra = opWallExt.nearSingInt(wallsExt,etaExt,DLPfun,[],...
    NearWext2T,kernelExt,kernelDirectExt,tracers,false,false);

DLP = o.wallDLPint;
for k = 1:nvbd-1
  DLP(:,:,k) = DLP(:,:,k) + jump*eye(2*NbdInt);
end
DLPfun = @(X) opWallInt.exactStokesDLdiag(wallsInt,DLP,X);
FwallInt2Tra = opWallInt.nearSingInt(wallsInt,etaInt,DLPfun,[],...
    NearWint2T,kernelInt,kernelDirectInt,tracers,false,false);

% compute the velocity on the tracers due to the
% solid walls
Fwall2Tra = FwallExt2Tra + FwallInt2Tra;

% velocity due to the stokeslets and rotlets    
FLets2Tra = zeros(2*Ntra,1);
for k = 2:nvbd
  stokeslet = RS(1:2,k);
  rotlet = RS(3,k);
  FLets2Tra = FLets2Tra + o.RSlets(Xtra,wallsInt.center(:,k-1),...
    stokeslet,rotlet);
end  

  
% END OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY


% START OF VELOCITY DUE TO VISCOSITY CONTRAST
if any(viscCont ~= 1)
  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
    kernelDirect = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLnewfmm;
    kernelDirect = @op.exactStokesDL;
  end

  % need double-layer potential matrix for doing near-singular
  % integration
  o.D = op.stokesDLmatrix(vesicle);
  
  DLP = o.D;
  % Get the DLP for inside&outside
  DLP_out = DLP;
  DLP_in  = DLP;
  for k=1:nv
    % if tracer is outside a vesicle, jump:  
    jump_out = 1/2*(1-viscCont(k));
    % if tracer is inside a vesicle, jump:
    jump_in  = -jump_out;
      
    % Based on jump, DLP is computed for inside/outside
    DLP_out(:,:,k) = DLP_out(:,:,k) + jump_out*eye(2*N);
    DLP_in (:,:,k) = DLP_in (:,:,k) + jump_in *eye(2*N);
      
    DLPfun_out = @(X) op.exactStokesDLdiag(vesicle,DLP_out,X);
    DLPfun_in  = @(X) op.exactStokesDLdiag(vesicle,DLP_in ,X);
  end
  
  % Initialize Fvisc2Tra
  Fvisc2Tra = zeros(2*Ntra,1);
  % Determine if the points inside or outside
  InOut = vesicle.sortPts(Xtra,true);
  InOut_ext = [InOut;InOut]; 
  % since size(InOut) = size(Xtra)/2 (only for one component) 
    
  % OUTSIDE: build tracers and compute Fvisc2Tra
  %----------------------------------------------------------------------
  Xtra_out = Xtra(InOut_ext == 0);
  tracers.N = numel(Xtra_out)/2;
  tracers.X = Xtra_out;
  [~,NearV2T_out] = vesicle.getZone(tracers,2);
  
  % Evaulate velocity due to vesicle's viscosity contrast on 
  % tracers with near-singular integration
  Fvisc2Tra(InOut_ext==0) = op.nearSingInt(vesicle,u,DLPfun_out,[],...
    NearV2T_out,kernel,kernelDirect,tracers,false,false);
  
    
  % INSIDE: build tracers and compute Fvisc2Tra
  %----------------------------------------------------------------------
  Xtra_in = Xtra(InOut_ext == 1);
  tracers.N = numel(Xtra_in)/2;
  tracers.X = Xtra_in;
  [~,NearV2T_in] = vesicle.getZone(tracers,2);
    
  % Evaulate velocity due to vesicle's viscosity contrast on 
  % tracers with near-singular integration
  Fvisc2Tra(InOut_ext==1) = op.nearSingInt(vesicle,u,DLPfun_in,[],...
    NearV2T_in,kernel,kernelDirect,tracers,false,false);
    
  % Back to original
  tracers.N = numel(Xtra)/2;
  tracers.X = Xtra;
else
  % no velocity due to viscosity contrast  
  Fvisc2Tra = zeros(2*Ntra,1);
end
% END OF VELOCITY DUE TO VISCOSITY CONTRAST


% velocity is the sum of velocities due to vesicles, solid walls,
% rotlet and stokeslets, and viscosity contrast
vel = Fves2Tra + Fwall2Tra + FLets2Tra + Fvisc2Tra;

%vel = Fwall2Tra + FLets2Tra;
% DEBUG: FOR TESTING TRACERS WITH NO VESICLES
%r2 = Xtra(1).^2 + Xtra(2).^2;
%speed = -1/3 + 400/3/r2;
%velTrue = speed*[-Xtra(2);Xtra(1)];
% This is the true velocity for the couette apparatus when both
% boundaries are centered at the origin, the inner boundary rotates
% once every 2*pi time units and the outer boundary is fixed

end % tracersVelDLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = tracersVelCouette(o,X,sigma,u,kappa,viscCont,walls,eta,...
        RS,Xtra,NearW2T)
% vel = tracersVel(o,X,sigma,u,kappa,viscCont,...
%    walls,wallsInt,wallsExt,eta,etaInt,etaExt,RS,Xtra) computes the
% velocity at the set of points Xtra due to the vesicle, viscosity
% contrast, solid walls, and stokeslets and rotlets.  The particles Xtra
% are treated passively

% build a structure for tracers which mimics the class capsules.
% Can't use capsules unless there is an even number of tracer
% points (because of FFTs), so we don't want to build the class
% as we do for vesicle and walls
vesicle = capsules(X,sigma,[],kappa,viscCont,o.antiAlias);
vesicle.setUpRate();
tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;

% build near-singular integration structures for vesicle to tracer and
% wall to tracer interactions
[~,NearV2T] = vesicle.getZone(tracers,2);

Ntra = size(Xtra,1)/2; % Number of tracers
N = size(X,1)/2; % Number of points
nv = size(X,2); % Number of vesicles

nvbd = walls.nv;
Nbd = walls.N;

% compute traction jump
f = vesicle.tracJump(X,sigma);

% poten classes
op = o.op;
opWall = o.opWall;


%% START OF VELOCITY DUE TO VESICLES
if ~o.fmm
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
  kernelDirect = @op.exactStokesSL;
end

% need single-layer potential matrix for doing near-singular
% integration
o.Galpert = op.stokesSLmatrix(vesicle);
SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);

% Evaulate velocity due to vesicle on tracers with near-singular
% integration
Fves2Tra = op.nearSingInt(vesicle,f,SLP,[],...
   NearV2T,kernel,kernelDirect,tracers,false,false);
% END OF VELOCITY DUE TO VESICLES

% START OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY
if ~o.fmmDLP
 kernel = @opWall.exactStokesDL;
else
 kernel = @opWall.exactStokesDLnewfmm;
end
kernelDirect = @opWall.exactStokesDL;   

% compute the velocity on the tracers due to the
% solid walls
jump = -1/2;
DLP = o.wallDLP + jump*eye(2*Nbd);
DLPfun = @(X) opWall.exactStokesDLdiag(walls,DLP,X);
Fwall2Tra = opWall.nearSingInt(walls,eta,DLPfun,[],...
    NearW2T,kernel,kernelDirect,tracers,false,false);

% velocity due to the stokeslets and rotlets    
FLets2Tra = zeros(2*Ntra,1);
for k = 2:nvbd
  stokeslet = RS(1:2,k);
  rotlet = RS(3,k);
  FLets2Tra = FLets2Tra + o.RSlets(Xtra,walls.center(:,k),...
    stokeslet,rotlet);
end  

  
% END OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY


% START OF VELOCITY DUE TO VISCOSITY CONTRAST
if any(viscCont ~= 1)
  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
    kernelDirect = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLnewfmm;
    kernelDirect = @op.exactStokesDL;
  end

  % need double-layer potential matrix for doing near-singular
  % integration
  o.D = op.stokesDLmatrix(vesicle);
  
  DLP = o.D;
  % Get the DLP for inside&outside
  DLP_out = DLP;
  DLP_in  = DLP;
  for k=1:nv
    % if tracer is outside a vesicle, jump:  
    jump_out = 1/2*(1-viscCont(k));
    % if tracer is inside a vesicle, jump:
    jump_in  = -jump_out;
      
    % Based on jump, DLP is computed for inside/outside
    DLP_out(:,:,k) = DLP_out(:,:,k) + jump_out*eye(2*N);
    DLP_in (:,:,k) = DLP_in (:,:,k) + jump_in *eye(2*N);
      
    DLPfun_out = @(X) op.exactStokesDLdiag(vesicle,DLP_out,X);
    DLPfun_in  = @(X) op.exactStokesDLdiag(vesicle,DLP_in ,X);
  end
  
  % Initialize Fvisc2Tra
  Fvisc2Tra = zeros(2*Ntra,1);
  % Determine if the points inside or outside
  InOut = vesicle.sortPts(Xtra,true);
  InOut_ext = [InOut;InOut]; 
  % since size(InOut) = size(Xtra)/2 (only for one component) 
    
  % OUTSIDE: build tracers and compute Fvisc2Tra
  %----------------------------------------------------------------------
  Xtra_out = Xtra(InOut_ext == 0);
  tracers.N = numel(Xtra_out)/2;
  tracers.X = Xtra_out;
  [~,NearV2T_out] = vesicle.getZone(tracers,2);
  
  % Evaulate velocity due to vesicle's viscosity contrast on 
  % tracers with near-singular integration
  Fvisc2Tra(InOut_ext==0) = op.nearSingInt(vesicle,u,DLPfun_out,[],...
    NearV2T_out,kernel,kernelDirect,tracers,false,false);
  
    
  % INSIDE: build tracers and compute Fvisc2Tra
  %----------------------------------------------------------------------
  Xtra_in = Xtra(InOut_ext == 1);
  tracers.N = numel(Xtra_in)/2;
  tracers.X = Xtra_in;
  [~,NearV2T_in] = vesicle.getZone(tracers,2);
    
  % Evaulate velocity due to vesicle's viscosity contrast on 
  % tracers with near-singular integration
  Fvisc2Tra(InOut_ext==1) = op.nearSingInt(vesicle,u,DLPfun_in,[],...
    NearV2T_in,kernel,kernelDirect,tracers,false,false);
    
  % Back to original
  tracers.N = numel(Xtra)/2;
  tracers.X = Xtra;
else
  % no velocity due to viscosity contrast  
  Fvisc2Tra = zeros(2*Ntra,1);
end
% END OF VELOCITY DUE TO VISCOSITY CONTRAST


% velocity is the sum of velocities due to vesicles, solid walls,
% rotlet and stokeslets, and viscosity contrast
vel = Fves2Tra + Fwall2Tra + FLets2Tra + Fvisc2Tra;

%vel = Fwall2Tra + FLets2Tra;
% DEBUG: FOR TESTING TRACERS WITH NO VESICLES
%r2 = Xtra(1).^2 + Xtra(2).^2;
%speed = -1/3 + 400/3/r2;
%velTrue = speed*[-Xtra(2);Xtra(1)];
% This is the true velocity for the couette apparatus when both
% boundaries are centered at the origin, the inner boundary rotates
% once every 2*pi time units and the outer boundary is fixed

end % tracersVelCouette
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vInf,vSec] = bgFlow(o,X,varargin)
% [vInf,vSec] = bgFlow(X,varargin) computes the velocity field at the 
% points X. vInf is either background or wall velocity. vSec is outputted
% if there are two sets of walls with different discretizations. Then, 
% vSec is the velocity for the inner walls, vInf for the outer walls.
% Default flow is shear with magnitude 1.  Other flows are relaxation,
% extensional, parabolic, Taylor-Green, cubic, invParabolic.  Some flows
% have an option to input a strength of the flow (ie. shear rate,
% extensional rate, etc).  Flows are given by:
%     cubic:          (y^3,x)
%     relaxation:     (0,0)
%     extensional:    (-x,y)
%     parabolic:      (k(1-y/r)^2,0)
%     invParabolic:   (ky^2,0)
%     rotate:         (2y,-x)
%     taylorGreen:    (sin(x)cos(y),-cos(x)sin(y))
%     shear:          (ky,0)
%     choke:          poeusille-like flow at intake and outtake
%     doublechoke:    same as choke
%     couette:        rotating boundary
%     doubleCouette   two rotating boundaries
%     quadCouette     four rotation boundaries
%     doubleFlower    same as doubleCouette
%     figureEight     tangential boundary condition
%     tube            poiseuille flow in a tube 
%     DLD             a period of a DLD device
% k is a 'flow strength', and r controls the width of the parabola
% in the parabolic flow

N = size(X,1)/2; % number of points per vesicle
nv = size(X,2); % number of vesicles

oc = curve;

% Separate out x and y coordinates
[x,y] = oc.getXY(X);

% speed of the background velocity
speed = varargin{find(strcmp(varargin,'Speed'))+1};    

if any(strcmp(varargin,'vortexSize'))
vortexSize = varargin{find(strcmp(varargin, 'vortexSize'))+1};
else
vortexSize = 1;
end
% usually we do not have diffDiscWalls, so set vSec empty
vSec = [];

if any(strcmp(varargin,'cubic'))
  vInf = [y.^3;zeros(N,nv)];

elseif any(strcmp(varargin,'relaxation'))
  vInf = zeros(2*N,nv); 

elseif any(strcmp(varargin,'extensional'))
  vInf = [-x;y];

elseif any(strcmp(varargin,'parabolic'))
  R = find(strcmp(varargin,'R'));
  if isempty(R)
%    R = 100;
      R = 10;
  else
    R = varargin{R+1};
  end
  UM = find(strcmp(varargin,'Umax'));
  if isempty(UM)
%    UM = 1e5; % default value for strength of flow
    UM = 10;
  else
    UM = varargin{UM+1};
  end
  vInf = [UM*(1-(y/R).^2);zeros(N,nv)];

elseif any(strcmp(varargin,'invParabolic'))
  k = find(strcmp(varargin,'k'));
  if isempty(k)
    k = 1; % default value for strength of flow
  else
    k = varargin{k+1};
  end
  vInf = [k*y.^2;zeros(N,nv)];

elseif any(strcmp(varargin,'rotate'))
  vInf = [y;-x];

elseif any(strcmp(varargin,'freeCouette'))
  % defined in 10<= r <= 20
  velx = y/3.*(1-400./(x.^2+y.^2));
  vely = -x/3.*(1-400./(x.^2+y.^2));
  r = sqrt(x.^2+y.^2);
  velx(r<=10) = 0; vely(r<=10) = 0;
  velx(r>=20) = 0; vely(r>=20) = 0;
   
  vInf = [velx;vely]; 

elseif any(strcmp(varargin,'freeCouetteLarge'))
  % defined in 5<= r <= 25
  % v = 10 on the inner cylinder  
  velx = 2*y/24.*(1-625./(x.^2+y.^2));
  vely = -2*x/24.*(1-625./(x.^2+y.^2));

  vInf = [velx;vely];

elseif any(strcmp(varargin,'taylorGreen'))
  vInf = vortexSize*[sin(x/vortexSize * pi).*cos(y/vortexSize * pi);-cos(x/vortexSize * pi).*sin(y/vortexSize * pi)];

elseif any(strcmp(varargin,'shear'))
  k = find(strcmp(varargin,'k'));
  if isempty(k)
    k = 1; % default value for strength of flow
  else
    k = varargin{k+1};
  end
  vInf = [k*y;zeros(N,nv)];
  
elseif any(strcmp(varargin,'random_shear'))
  rng('shuffle');
  vInf = varargin{find(strcmp(varargin,'vInf'))+1}; % this might be given
  if isempty(vInf)
    speed = varargin{find(strcmp(varargin,'Speed'))+1};
    numRandFuncs = varargin{find(strcmp(varargin,'numRandFuncs'))+1};
    freqRange = varargin{find(strcmp(varargin,'freqRange'))+1};
    if isempty(o.randFreqs)
      o.randFreqs = (freqRange(1)+diff(freqRange)*rand(numRandFuncs,1))/o.finalTime;
      o.randAmps = (numRandFuncs:-1:1)'-rand(numRandFuncs,1);
      o.randAmps = o.randAmps/sum(o.randAmps);
    end
    vInfRand = @(t,y,freq,A) A*sin(2*pi*freq*t);%.*cos(2*pi*freq*y);

    ushear = y;
    vshear = zeros(size(y));
    urand = zeros(size(y));
    vrand = zeros(size(y));
    for nf = 1 : numRandFuncs
      urand = urand + ...
          vInfRand(o.currentTime,y,o.randFreqs(nf),o.randAmps(nf));      
      vrand = vrand + ...
          vInfRand(o.currentTime,zeros(size(y)),o.randFreqs(nf),0.5*o.randAmps(nf));
    end

    vInf = [ushear+urand;vshear+vrand];
%     figure(2);clf;hold on;
%     plot(x,y,'r')
%     quiver(x,y,vInf(1:end/2),vInf(end/2+1:end))
%     axis equal
    
  end % if isempty(vInf)
  
  
elseif any(strcmp(varargin,'oscil_shear1'))
  speed = 1;
  u_oscil = @(t,n) sin(2*pi*n*t*3/o.finalTime);
  ushear = y;
  vshear = zeros(size(y));
  Nmodes = 4;
  uosc = zeros(size(ushear));

  for nm = 1 : Nmodes
    uosc = uosc + ushear.*2/Nmodes.*u_oscil(o.currentTime,nm);
%     vshear = vshear + rand/5*1/Nmodes*u_oscil(o.currentTime,nm);
  end
  vInf = [uosc;vshear];    
  
elseif any(strcmp(varargin,'oscil_shear3'))
  speed = 5;
  u_oscil = @(t,n) sin(2*pi*n*t*2/o.finalTime);
  ushear = y;
  vshear = zeros(size(y));
  Nmodes = 2;
  uosc = zeros(size(ushear));
  
  for nm = 1 : Nmodes
    uosc = uosc + ushear.*2/Nmodes.*u_oscil(o.currentTime,nm);
%     vshear = vshear + rand/5*1/Nmodes*u_oscil(o.currentTime,nm);
  end
  vInf = [uosc;vshear];  
  
elseif (any(strcmp(varargin,'choke')) || ...
      any(strcmp(varargin,'doublechoke')) || ...
      any(strcmp(varargin,'choke2')))
  vInf = zeros(2*N,nv);
  ind = abs(x)>7;
  vx = exp(1./((y(ind)/max(y)).^2-1))/exp(-1);
  % typical mollifer so that velocity decays smoothly to 0
  vx(vx==Inf) = 0;
  vInf(ind,:) = vx;
elseif (any(strcmp(varargin,'tube')))
  vInf = zeros(2*N,nv);
  ind = abs(x)>0.8*max(x);
  vx = exp(1./((y(ind)/max(y)).^2-1))/exp(-1);
  % typical mollifer so that velocity decays smoothly to 0
  vx(vx==Inf) = 0;
  vInf(ind,:) = vx;
    
  
elseif any(strcmp(varargin,'couette'));
  vInf = [zeros(2*N,1) 1*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))]];
  
elseif any(strcmp(varargin,'couetteOuter'));
  vInf = [1*[-y(:,1)+mean(y(:,1));x(:,1)-mean(x(:,1))] zeros(2*N,1)];  

elseif any(strcmp(varargin,'couette10'));
  vInf = [zeros(2*N,1) 10*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))]];
  
elseif any(strcmp(varargin,'couette100'));
  vInf = [zeros(2*N,1) 100*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))]];
  
elseif (any(strcmp(varargin,'doubleCouette')) || ...
      any(strcmp(varargin,'doubleFlower')));
  vInf = [zeros(2*N,1) 1*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))] ...
      -[y(:,3)-mean(y(:,3));-x(:,3)+mean(x(:,3))]];

elseif (any(strcmp(varargin,'quadCouette')));
  vInf = [zeros(2*N,1) ...
      -[y(:,2)-mean(y(:,2));-x(:,2)+mean(x(:,2))] ...
      -[y(:,3)-mean(y(:,3));-x(:,3)+mean(x(:,3))] ...
      +[y(:,4)-mean(y(:,4));-x(:,4)+mean(x(:,4))] ...
      +[y(:,5)-mean(y(:,5));-x(:,5)+mean(x(:,5))]];

elseif any(strcmp(varargin,'cylinder'))
%  theta = (0:N-1)'*2*pi/N;
%  vInf = [cos(10*theta);sin(2*theta)];
  vInf = 1*[-y+mean(y);x-mean(x)];

elseif any(strcmp(varargin,'figureEight'))
  oc = curve;
  [~,vInf,~] = oc.diffProp([x;y]);
%  vInf(1:end/2,:) = 3*vInf(1:end/2,:);

  sup = find(abs(x)<=1 & y>0);
  sdown = find(abs(x)<=1 & y<0);
  omega = linspace(-1,1,numel(sup)+2)';
  omega = omega(2:end-1);
  mollifier = 4*exp(1./(omega.^2 - 1))+1;
  vInf(sup,:) = vInf(sup,:) .* mollifier;
  vInf(sdown,:) = vInf(sdown,:) .* mollifier;
  % increase the velocity in a smooth fashion near the middle
  % of the solid walls

elseif any(strcmp(varargin,'shear'))
  vInf = [y;zeros(N,nv)];

elseif any(strcmp(varargin,'diffuser'));
  vInf = zeros(2*N,nv);
  ind = abs(x(:,1))>9;
  vx = exp(1./((y(ind,1)/max(y(ind,1))).^2-1))/exp(-1);
  % typical mollifer so that velocity decays smoothly to 0
  vx(vx==Inf) = 0;
  vInf(ind,1) = vx;

elseif any(strcmp(varargin,'microfluidic'));
  oc = curve;
  [~,tangent,~] = oc.diffProp(X); 
  vInf = tangent;
  vInf(:,1) = 0*vInf(:,1);
  vInf(:,2) = +1*vInf(:,2);
  vInf(:,3) = -1*vInf(:,3);
  vInf(:,4) = -1*vInf(:,4);
  vInf(:,5) = +1*vInf(:,5);
  
elseif any(strcmp(varargin,'DLD'))
  % interior walls
  XwallsInt = varargin{find(strcmp(varargin,'intWalls'))+1};
  
  idebug = 0;
  if idebug
    figure(3);clf;hold on;
    plot(x,y,'k','linewidth',2);
    axis equal;
    plot(XwallsInt(1:end/2,:),XwallsInt(end/2+1:end,:),'k','linewidth',2)
  end
  
  vSec = zeros(size(XwallsInt)); % 0 velocity on the interior walls
  % number of rows in DLD
  nrow = varargin{find(strcmp(varargin,'nrow'))+1};    
  % number of columns in DLD
  ncol = varargin{find(strcmp(varargin,'ncol'))+1};
  % Gap in x direction
  Dx = varargin{find(strcmp(varargin,'GapX'))+1};
  % Gap in y direction
  Dy = varargin{find(strcmp(varargin,'GapY'))+1};
  % row-shift fraction
  epsilon = varargin{find(strcmp(varargin,'epsilon'))+1};
  % period of a device 
  periodN = ceil(1/epsilon);
  
  Dpostx = varargin{find(strcmp(varargin,'Dpostx'))+1};
  Dposty = varargin{find(strcmp(varargin,'Dposty'))+1};
  
  
  % row-shift
  delLat = (Dy+Dposty)*epsilon;
  
  % left hand side of DLD
  pnrow = 2*nrow;
  intervals = zeros(pnrow+1,2);
  intervals(:,1) = Dposty*linspace(pnrow/2,-(pnrow/2),pnrow+1)';
  intervals(:,2) = intervals(:,1);
  intervals(:,1) = intervals(:,1)+linspace((pnrow+1)/2,-...
      (pnrow+1)/2+1,pnrow+1)'*Dy;
  intervals(:,2) = intervals(:,1)-Dy;
  
  intervalsIn = intervals - delLat;
  
  vInf = zeros(2*N,nv);
  indx = (x<=0.99*min(x));
  for jint = 1 : pnrow+1
    hmax = intervalsIn(jint,1); hmin = intervalsIn(jint,2);
    indy = y<=hmax & y>=hmin;
    ind = indx & indy;
    
    ymean = (hmax+hmin)/2;
    vy = 1-((y(ind)-ymean)/(hmax-ymean)).^2;
    vInf(ind,1) = vy;
    
    if idebug
      figure(3);
      plot(x(ind),y(ind),'bo','markersize',5)    
      figure(4);clf;
      plot(vy,y(ind),'g','linewidth',2)
      xlabel('velocity')
      ylabel('y')
      pause
    end
    
  end
  
  % right end
%   moveUp = ncol-periodN;
  moveUp = ncol;

  intervalsOut = intervals + moveUp*delLat;
  
  % out for rotated
  %   intervalsOut = intervals + epsilon*L;

  indx = (x>=0.99*max(x));
  %   intervals = intervals+L*epsilon;
  for jint = 1 : pnrow+1
    hmax = intervalsOut(jint,1); hmin = intervalsOut(jint,2);
    indy = y<=hmax & y>=hmin;
    ind = indx & indy;
    
    ymean = (hmax+hmin)/2;
    vy = 1-((y(ind)-ymean)/(hmax-ymean)).^2;
    vInf(ind,1) = vy;
    
    if idebug
      figure(3);
      plot(x(ind),y(ind),'bo','markersize',5)    
      figure(4);clf;
      plot(vy,y(ind),'g','linewidth',2)
      xlabel('velocity')
      ylabel('y')
      pause
    end
    
  end
  
% ===================== ROTATED DLD FOR OPTIMIZATION ======================  
elseif any(strcmp(varargin,'rotDLDiter'))
  % This is for small rotated DLD device with 5 rows.
  XwallsInt = varargin{find(strcmp(varargin,'intWalls'))+1};
  vSec = zeros(size(XwallsInt)); % 0 velocity on the interior walls
  idebug = 1;
  
  if idebug
    figure(3);clf;hold on;
    plot([x;x(1)],[y;y(1)],'k','linewidth',2);
    plot([XwallsInt(1:end/2,:);XwallsInt(1,:)],...
        [XwallsInt(end/2+1:end,:);XwallsInt(end/2+1,:)],'k','linewidth',2)
    axis equal;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'ztick',[]);

    set(gca,'xcolor','w');
    set(gca,'ycolor','w');
    set(gca,'zcolor','w');
    set(gca,'visible','off');
  end
  
  % Velocity is already given for top and bottom walls
  if any(strcmp(varargin,'velTop'))
    velTop = varargin{find(strcmp(varargin,'velTop'))+1};
    Xtop = varargin{find(strcmp(varargin,'velTop'))+2}; %[x;y] for top
    Xbot = varargin{find(strcmp(varargin,'velTop'))+3}; %[x;y] for bottom
    % Velocity is already given for right and left walls
    velSide = varargin{find(strcmp(varargin,'velSide'))+1};
    Xleft = varargin{find(strcmp(varargin,'velSide'))+2};%[x;y] for left
    Xright = varargin{find(strcmp(varargin,'velSide'))+3}; %[x;y] for right
  else
    velTop = [];
    velSide = [];
    % number of rows in DLD
    nrow = varargin{find(strcmp(varargin,'nrow'))+1};    
    % number of columns in DLD
    ncol = varargin{find(strcmp(varargin,'ncol'))+1};
    % Gap in y direction
    Dx = varargin{find(strcmp(varargin,'GapX'))+1};
    % Gap in y direction
    Dy = varargin{find(strcmp(varargin,'GapY'))+1};
    % row-shift fraction
    epsilon = varargin{find(strcmp(varargin,'epsilon'))+1};

  end
  
  
  if isempty(velTop)    
    % Find Dpostx and Dposty, Dpost is useless
    Dpostx = varargin{find(strcmp(varargin,'Dpostx'))+1};
    Dposty = varargin{find(strcmp(varargin,'Dposty'))+1};
    
    % row-shift
    delLat = (Dy+Dposty)*epsilon;
  
    % left hand side of DLD
    pnrow = nrow+2;
    nGaps = pnrow-1;
    intervals = zeros(nGaps,2);
    intervals(:,1) = Dposty*linspace((nGaps-1)/2,-(nGaps-1)/2,nGaps)';
    intervals(:,1) = intervals(:,1)+linspace(nGaps/2,1-nGaps/2,nGaps)'*Dy;
    intervals(:,2) = intervals(:,1)-Dy;
  
    intervalsIn = intervals - delLat;
  
    vInf = zeros(2*N,nv);
    indx = (x<=0.99*min(x));
    for jint = 1 : nGaps
      hmax = intervalsIn(jint,1); hmin = intervalsIn(jint,2);
      indy = y<=hmax & y>=hmin;
      ind = indx & indy;
    
      ymean = (hmax+hmin)/2;
      vy = 1-((y(ind)-ymean)/(hmax-ymean)).^2;
      vInf(ind,1) = vy;
    
      if idebug
        figure(3);
        plot(x(ind),y(ind),'bo','markersize',6,'markerfacecolor','b')    
        figure(4);clf;
        plot(vy,y(ind),'r','linewidth',2)
        xlabel('velocity')
        ylabel('y')
        title('inlet')
        pause
      end
    
    end
  
    % right end
    moveUp = ncol;

    intervalsOut = intervals + moveUp*delLat;

    indx = (x>=0.99*max(x));
    %   intervals = intervals+L*epsilon;
    for jint = 1 : nGaps
      hmax = intervalsOut(jint,1); hmin = intervalsOut(jint,2);
      indy = y<=hmax & y>=hmin;
      ind = indx & indy;
    
      ymean = (hmax+hmin)/2;
      vy = 1-((y(ind)-ymean)/(hmax-ymean)).^2;
      vInf(ind,1) = vy;
    
      if idebug
        figure(3);
        plot(x(ind),y(ind),'bo','markersize',6,'markerfacecolor','b')    
        figure(4);clf;
        plot(vy,y(ind),'r','linewidth',2)
        xlabel('velocity')
        ylabel('y')
        title('outlet')
        pause
      end
    
    end %for jint
  else %top, bottom and side walls have velocity given
    % no need to scale speed b/c we interpolate from already scaled vel.
    %speed = 1;
    
    
    idebug = 0;
    
    xtop = Xtop(1:end/2); ytop = Xtop(end/2+1:end);
    xbot = Xbot(1:end/2); ybot = Xbot(end/2+1:end);
    
    ylft = Xleft(end/2+1:end); yrgt = Xright(end/2+1:end);
    
    % Points on the left, right, top and bottom of the exterior wall
    idcsTop = N/8+1:N/4+N/8;
    idcsLeft = N/4+N/8+1:N/2+N/8;
    idcsBot = N/2+N/8+1:3*N/4+N/8;
    idcsRight = [3*N/4+N/8+1:N 1:N/8];

%     idcsLeft = abs(x-min(Xleft(1:end/2)))<=1e-1;
%     idcsRight = abs(x-max(Xright(1:end/2)))<=1e-1;
%     
%     % Points that are not included in left and right boundaries are
%     % included in top and bottom boundaries
%     idcsTop = (y>=max(y(idcsLeft)) & ~idcsLeft & ~idcsRight);
%     idcsBot = (y<=min(y(idcsRight)) & ~idcsLeft & ~idcsRight);

    % interpolate velocity for those points
    vInf = zeros(2*N,nv);
    
    % top and bottom walls
    interpMeth = 'pchip';
    velxTop = interp1(xtop,velTop(1:end/2),x(idcsTop),interpMeth);
    velyTop = interp1(xtop,velTop(end/2+1:end),x(idcsTop),interpMeth);
%     ids = find(idcsTop==1);
    ids = idcsTop;
    vInf(ids) = velxTop; vInf(ids+N) = velyTop;
    
    velxBot = interp1(xbot,velTop(1:end/2),x(idcsBot),interpMeth);
    velyBot = interp1(xbot,velTop(end/2+1:end),x(idcsBot),interpMeth);
%     ids = find(idcsBot==1);
    ids = idcsBot;
    vInf(ids) = velxBot; vInf(ids+N) = velyBot;
    
    
    % right and left walls
    velxRgt = interp1(yrgt,velSide(1:end/2),y(idcsRight),interpMeth);
    velyRgt = interp1(yrgt,velSide(end/2+1:end),y(idcsRight),interpMeth);
%     ids = find(idcsRight==1);
    ids = idcsRight;
    vInf(ids) = velxRgt; vInf(ids+N) = velyRgt;

    velxLft = interp1(ylft,velSide(1:end/2),y(idcsLeft),interpMeth);
    velyLft = interp1(ylft,velSide(end/2+1:end),y(idcsLeft),interpMeth);
%     ids = find(idcsLeft==1);
    ids = idcsLeft;
    vInf(ids) = velxLft; vInf(ids+N) = velyLft;

    
    
    % ADD DEBUGGING, PLOT RIGHT,LEFT,TOP,BOTTOM WALL INDICES
    % PLOT VELOCITIES THERE
    
    if idebug
      disp('Maximum velocity magnitude in the interpolated boundary conditions on side walls')
      max(sqrt(velxRgt.^2 + velyRgt.^2))
      disp('Maximum velocity magnitude in the given boundary conditions on side walls')
      max(sqrt(velSide(1:end/2).^2+velSide(end/2+1:end).^2))
      disp('Maximum velocity magnitude in the interpolated boundary conditions on top and bottom walls')
      max(sqrt(velxTop.^2+velyTop.^2))
      disp('Maximum velocity magnitude in the given boundary conditions on top and bottom walls')
      max(sqrt(velTop(1:end/2).^2+velTop(end/2+1:end).^2))  
      
      figure(7); clf;hold on
      plot([x;x(1)],[y;y(1)],'k-o','linewidth',2,'markersize',6,'markerfacecolor','k')
      plot([XwallsInt(1:end/2,:);XwallsInt(1,:)],...
          [XwallsInt(end/2+1:end,:);XwallsInt(end/2+1,:)],'k','linewidth',2)
      
      axis equal
      set(gca,'xtick',[]);
      set(gca,'ytick',[]);
      set(gca,'ztick',[]);

      set(gca,'xcolor','w');
      set(gca,'ycolor','w');
      set(gca,'zcolor','w');
      set(gca,'visible','off');
      
      plot(x(idcsTop),y(idcsTop),'g*','markersize',6,'markerfacecolor','g')
      plot(x(idcsBot),y(idcsBot),'g*','markersize',6,'markerfacecolor','g')
      plot(x(idcsRight),y(idcsRight),'md','markersize',6,'markerfacecolor','m')
      plot(x(idcsLeft),y(idcsLeft),'md','markersize',6,'markerfacecolor','m')
      
      plot(Xleft(1:end/2),Xleft(end/2+1:end),'c','linewidth',3)
      plot(Xright(1:end/2),Xright(end/2+1:end),'c','linewidth',3)
      plot(xtop,ytop,'c','linewidth',3)
      plot(xbot,ybot,'c','linewidth',3)
      
      pause
      
      
      % Plot the given velocity first
      
      quiver(xtop,ytop,velTop(1:end/2),velTop(end/2+1:end),'r','AutoScale','off')
      quiver(xbot,ybot,velTop(1:end/2),velTop(end/2+1:end),'r','AutoScale','off');
      quiver(Xleft(1:end/2),ylft,velSide(1:end/2),velSide(end/2+1:end),'r','AutoScale','off')
      quiver(Xright(1:end/2),yrgt,velSide(1:end/2),velSide(end/2+1:end),'r','AutoScale','off')
      
      % Plot the interpolated velocity now
      quiver(x(idcsTop),y(idcsTop),velxTop,velyTop,'b','AutoScale','off')
      quiver(x(idcsBot),y(idcsBot),velxBot,velyBot,'b','AutoScale','off')
      quiver(x(idcsRight),y(idcsRight),velxRgt,velyRgt,'b','AutoScale','off')
      quiver(x(idcsLeft),y(idcsLeft),velxLft,velyLft,'b','AutoScale','off')
      pause
    end
    
  end
 
% ===================== ROTATED DLD FOR OPTIMIZATION ======================  
elseif any(strcmp(varargin,'rotDLD'))
  % This is for small rotated DLD device with 5 rows and 4 columns
  % BC is given
  XwallsInt = varargin{find(strcmp(varargin,'intWalls'))+1};
  vSec = zeros(size(XwallsInt)); % 0 velocity on the interior walls
  
  vInf = varargin{find(strcmp(varargin,'velG'))+1};
  
  speed = 1;
else
  vInf = [y;zeros(N,nv)];
  % default flow is shear
end


vInf = vInf * speed;
end % bgFlow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GLpts = gaussLobatto(o,orderGL)
% GLpts = gaussLobatto(orderGL) loads the Guass- Lobatto points for the
% desired number of points 

if (orderGL == 2 || orderGL == -2)
  GLpts = [-1 1];
elseif (orderGL == 3 || orderGL == -3)
  GLpts = [-1 0 1];
elseif orderGL == 4
  GLpts = [-1 -0.447213595499958 0.447213595499958 1];
elseif orderGL == -4
  GLpts = [-1 -0.333333333333333 0.333333333333333 1];
elseif orderGL == 5
  GLpts = [-1 -0.654653670707977 0 0.654653670707977 1];
elseif orderGL == -5
  GLpts = [-1 -0.5 0 0.5 1];
elseif orderGL == 6
  GLpts = [-1 -0.765055323929465 -0.285231516480645 ...
      0.285231516480645 0.765055323929465 1];
elseif orderGL == -6
  GLpts = [-1 -0.6 -0.2 0.2 0.6 1];
elseif orderGL == 7
  GLpts = [-1 -0.830223896278567 -0.468848793470714 ...
    0 0.468848793470714 0.830223896278567 1];
elseif orderGL == -7
  GLpts = [-1 -0.6666666666666666  -0.3333333333333333 ...
    0 0.3333333333333333 0.6666666666666667 1];
else
  fprintf('**************************************\n')
  fprintf('NO GAUSS-LOBATTO POINTS FOR THIS ORDER\n')
  fprintf('**************************************\n')
  pause
end

end % gaussLobatto

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fint = lobattoInt(o,f)
% fint = lobattoInt(f) returns the integral of f from [-1,t] where t
% ranges from -1 to 1.  If orderGL > 0, the interior points are
% Gauss-Lobato points.  Otherwise, they are equispaced points.  The
% quadrature rules were found by first interpolating f with a
% polynomial and then integrating the polynomial exactly.  All
% coefficients are precomputed f is of size (N,nv,order) where N is the
% number of points that need to be integrated, nv is a second index
% which in this case corresponds to number of vesicles, and order is
% the number of time steps

nv = size(f,2);
% number of vesicles
fint = zeros(size(f));
% initialize the integral to zero

if (o.orderGL == 2 || o.orderGL == -2)
  t = [-1 1];

  for n = 1:2
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*0.25*(-f(:,k,1)+f(:,k,2));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*0.5*(f(:,k,1)+f(:,k,2));

      fint(:,k,n) = fint(:,k,n) + ...
          (0.75*f(:,k,1) + 0.25*f(:,k,2));
    end
  end
  % order 1 Gauss-Lobatto or equispaced


elseif (o.orderGL == 3 || o.orderGL == -3)
  t = [-1 0 1];

  for n = 1:3
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.1666666666666667*(f(:,k,1)+f(:,k,3)) - ...
          0.3333333333333333*f(:,k,2));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.25*(f(:,k,1)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*f(:,k,2);

      fint(:,k,n) = fint(:,k,n) + ...
          0.4166666666666667*f(:,k,1) + ...
          0.6666666666666667*f(:,k,2) - ...
          0.0833333333333333*f(:,k,3);
    end
  end
  % order 2 Gauss-Lobatto or equispaced


elseif (o.orderGL == 4)
  t = [-1 -0.447213595499958 0.447213595499958 1];

  for n = 1:4
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.15625*(f(:,k,1)-f(:,k,4)) + ...
          0.3493856214843422*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.2083333333333333*(f(:,k,1)+f(:,k,4)) - ...
          0.2083333333333333*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.0625*(f(:,k,1)-f(:,k,4)) - ...
          0.6987712429686845*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(-0.125*(f(:,k,1)+f(:,k,4)) + ...
          0.625*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          (0.1770833333333333*f(:,k,1) + ...
          0.7660522881510090*f(:,k,2) + ...
          0.0672810451823247*f(:,k,3) - ...
          0.0104166666666666*f(:,k,4));
    end
  end
  % order 3 Gauss-Lobatto

elseif (o.orderGL == -4)
  t = [-1 -0.333333333333333 0.333333333333333 1];

  for n = 1:4
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.140625*(f(:,k,1)-f(:,k,4)) + ...
          0.421875*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.1875*(f(:,k,1)+f(:,k,4)) - ...
          0.1875*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.03125*(f(:,k,1)-f(:,k,4)) - ...
          0.84375*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(-0.0625*(f(:,k,1)+f(:,k,4)) + ...
          0.5625*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          (0.234375*f(:,k,1) + ...
          0.796875*f(:,k,2) - ...
          0.046875*f(:,k,3) + ...
          0.015625*f(:,k,4));
    end
  end
  % order 3 equispaced

elseif (o.orderGL == 5)
  t = [-1 -0.654653670707977 0 0.654653670707977 1];

  for n = 1:5
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.175*(f(:,k,1)+f(:,k,5)) - ...
          0.4083333333333333*(f(:,k,2)+f(:,k,4)) + ...
          0.4666666666666666*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.21875*(f(:,k,1)-f(:,k,5)) + ...
          0.3341461444238633*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.125*(f(:,k,1)+f(:,k,5)) + ...
          0.6805555555555555*(f(:,k,2)+f(:,k,4)) - ...
          1.1111111111111111*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.1875*(f(:,k,1)-f(:,k,5)) - ...
          0.6682922888477265*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*f(:,k,3);

      fint(:,k,n) = fint(:,k,n) + ...
          0.08125*f(:,k,1) + ...
          0.6063683666460855*f(:,k,2) + ...
          0.3555555555555555*f(:,k,3) - ...
          0.0619239222016411*f(:,k,4) + ...
          0.01875*f(:,k,5);
    end
  end
  % order 4 Gauss-Lobatto

elseif (o.orderGL == -5)
  t = [-1 -0.5 0 0.5 1];

  for n = 1:5
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.1333333333333333*(f(:,k,1)+f(:,k,5)) - ...
          0.53333333333333333*(f(:,k,2)+f(:,k,4)) + ...
          0.8*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.16666666666666667*(f(:,k,1)-f(:,k,5)) + ...
          0.3333333333333333*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.0555555555555556*(f(:,k,1)+f(:,k,5)) + ...
          0.8888888888888889*(f(:,k,2)+f(:,k,4)) - ...
          1.6666666666666667*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.0833333333333333*(f(:,k,1)-f(:,k,5)) - ...
          0.6666666666666667*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*f(:,k,3);

      fint(:,k,n) = fint(:,k,n) + ...
          0.1611111111111111*f(:,k,1) + ...
          0.6888888888888889*f(:,k,2) + ...
          0.1333333333333333*f(:,k,3) + ...
          0.0222222222222222*f(:,k,4) - ...
          0.0055555555555556*f(:,k,5);
    end
  end
  % order 4 equi-spaced points

elseif (o.orderGL == 6)
  t = [-1 -0.765055323929465 -0.285231516480645 ...
      0.285231516480645 0.765055323929465 1];

  for n = 1:6
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.21875*(f(:,k,1)-f(:,k,6)) + ...
          0.5212094304495727*(f(:,k,2)-f(:,k,5)) - ...
          0.6310805056491861*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.2625*(f(:,k,1)+f(:,k,6)) - ...
          0.4785048595772274*(f(:,k,2)+f(:,k,5)) + ...
          0.2160048595772274*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.21875*(f(:,k,1)-f(:,k,6)) - ...
          0.8454202131918329*(f(:,k,2)-f(:,k,5)) + ...
          1.500687022042463*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.2916666666666666*(f(:,k,1)+f(:,k,6)) + ...
          0.8623909800799931*(f(:,k,2)+f(:,k,5)) - ...
          0.5707243134133265*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.03125*(f(:,k,1)-f(:,k,6)) + ...
          0.1272121350349483*(f(:,k,2)-f(:,k,5)) - ...
          1.108132527137368*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(0.0625*(f(:,k,1)+f(:,k,6)) - ...
          0.1946486423538423*(f(:,k,2)+f(:,k,5)) + ...
          0.6321486423538425*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          0.06458333333333333*f(:,k,1) + ...
          0.3862361258562357*f(:,k,2) + ...
          0.5159551992618339*f(:,k,3) + ...
          0.03890317777365201*f(:,k,4) - ...
          0.007761169558388661*f(:,k,5) + ...
          0.002083333333333333*f(:,k,6);
    end
  end
  % order 5 Gauss-Lobatto

elseif (o.orderGL == -6)
  t = [-1 -0.6 -0.2 0.2 0.6 1];

  for n = 1:6
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.1356336805555556*(f(:,k,1)-f(:,k,6)) + ...
          0.6781684027777778*(f(:,k,2)-f(:,k,5)) - ...
          1.3563368055555556*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.1627604166666667*(f(:,k,1)+f(:,k,6)) - ...
          0.48828125*(f(:,k,2)+f(:,k,5)) + ...
          0.3255208333333333*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.0813802083333333*(f(:,k,1)-f(:,k,6)) - ...
          1.0579427083333333*(f(:,k,2)-f(:,k,5)) + ...
          2.7669270833333333*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.108506944444444*(f(:,k,1)+f(:,k,6)) + ...
          0.8463541666666667*(f(:,k,2)+f(:,k,5)) - ...
          0.7378472222222222*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.005859375*(f(:,k,1)-f(:,k,6)) + ...
          0.0813802083333333*(f(:,k,2)-f(:,k,5)) - ...
          1.46484375*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(0.01171875*(f(:,k,1)+f(:,k,6)) - ...
          0.09765625*(f(:,k,2)+f(:,k,5)) + ...
          0.5859375*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          0.1260850694444444*f(:,k,1) + ...
          0.5588107638888889*f(:,k,2) + ...
          0.2278645833333333*f(:,k,3) + ...
          0.1193576388888889*f(:,k,4) - ...
          0.0379774305555556*f(:,k,5) + ...
          0.005859375*f(:,k,6);
    end
  end
  % order 5 equi-spaced points

elseif (o.orderGL == 7)
  t = [-1.000000000000000 -0.830223896278567 ...
      -0.468848793470714 0 0.468848793470714 ...
      0.830223896278567 1.000000000000000];

  for n = 1:7
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^7*(0.29464285714285718*(f(:,k,1)+f(:,k,7)) - ...
          0.71040995801242226*(f(:,k,2)+f(:,k,6)) + ...
          0.88719567229813686*(f(:,k,3)+f(:,k,5)) - ...
          0.94285714285714356*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.34375*(f(:,k,1)-f(:,k,7)) + ...
          0.68809921051219412*(f(:,k,2)-f(:,k,6)) - ...
          0.48528739061765716*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(-0.375*(f(:,k,1)+f(:,k,7)) + ...
          1.21320038050366995*(f(:,k,2)+f(:,k,6)) - ...
          2.09820038050367082*(f(:,k,3)+f(:,k,5)) + ...
          2.52000000000000176*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.46875*(f(:,k,1)-f(:,k,7)) - ...
          1.25903493358549612*(f(:,k,2)-f(:,k,6)) + ...
          1.22967339607367386*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.10416666666666660*(f(:,k,1)+f(:,k,7)) - ...
          0.36437739881046465*(f(:,k,2)+f(:,k,6)) + ...
          1.42687739881046537*(f(:,k,3)+f(:,k,5)) - ...
          2.3333333333333346*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.15625*(f(:,k,1)-f(:,k,7)) + ...
          0.45377223563440987*(f(:,k,2)-f(:,k,6)) - ...
          1.00348462029437623*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(1.0*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          0.55059523809523700e-1*f(:,k,1) + ...
          0.25557651111967517*f(:,k,2) + ...
          0.47497130544329094*f(:,k,3) + ...
          0.24380952380952357*f(:,k,4) - ...
          0.4322592423342813e-1*f(:,k,5) + ...
          0.2124953624189091e-1*f(:,k,6) - ...
          0.7440476190476161e-2*f(:,k,7);
    end
  end
  % order 6 Gauss-Lobatto

elseif (o.orderGL == -7)
  t = [-1 -0.6666666666666666  ...
      -0.3333333333333333 0 0.3333333333333333 ...
      0.6666666666666667 1];

  for n = 1:7
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^7*(0.1446428571428571*(f(:,k,1)+f(:,k,7)) - ...
          0.8678571428571429*(f(:,k,2)+f(:,k,6)) + ...
          2.1696428571428571*(f(:,k,3)+f(:,k,5)) - ...
          2.8928571428571429*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.16875*(f(:,k,1)-f(:,k,7)) + ...
          0.675*(f(:,k,2)-f(:,k,6)) - ...
          0.84375*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(-0.1125*(f(:,k,1)+f(:,k,7)) + ...
          1.35*(f(:,k,2)+f(:,k,6)) - ...
          4.3875*(f(:,k,3)+f(:,k,5)) + ...
          6.3*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.140625*(f(:,k,1)-f(:,k,7)) - ...
          1.125*(f(:,k,2)-f(:,k,6)) + ...
          1.828125*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.0166666666666667*(f(:,k,1)+f(:,k,7)) - ...
          0.225*(f(:,k,2)+f(:,k,6)) + ...
          2.25*(f(:,k,3)+f(:,k,5)) - ...
          4.0833333333333333*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.025*(f(:,k,1)-f(:,k,7)) + ...
          0.225*(f(:,k,2)-f(:,k,6)) - ...
          1.125*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(1.0*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          0.1019345238095238*f(:,k,1) + ...
          0.4821428571428571*f(:,k,2) + ...
          0.1727678571428571*f(:,k,3) + ...
          0.3238095238095238*f(:,k,4) - ...
          0.1084821428571429*f(:,k,5) + ...
          0.03214285714285714*f(:,k,6) - ...
          0.004315476190476190*f(:,k,7);
    end
  end
  % order 6 equi-spaced points
end % o.orderGL

end % lobattoInt

end % methods

end % classdef
