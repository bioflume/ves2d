%clear all; %clc

disp('One elliptical vesicles in a couette apparatus.');
disp('First-order semi-implicit time stepping.');
disp('Implicit vesicle-boundary interactions.');

%format long
prams.N = 192;                  % points per vesicle
prams.nv = 2;                  % number of vesicles
prams.T = 100;                % time horizon
prams.m = 1e6;                  % number of time steps
prams.Nbd = 384;                % number of points on solid wall
prams.Nbdcoarse = prams.Nbd;                % number of points on solid wall
prams.nvbd = 1;                % number of components to solid walls
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = 1e0;          % viscosity contrast
prams.gmresTol = 1e-10;         % GMRES tolerance
prams.errorTol = 1e-1;
% Maximum relative error in area and length before the simulation
% is stopped
prams.rtolArea = 5e-1;
prams.rtolLength = 5e-1;


options.farField = 'diffuser';    % Constricted domain
% background velocity or solid walls (which comes with a default
% boundary condition)
options.order = 1;          % time stepping order
options.inextens = 'method1'; 
% Inextensibility condition can be treated as either
% method1 - div x^{N+1} = div x^{N}
% method2 - div(SLP*tractionjump) = 0;
options.vesves = 'implicit';
% Discretization of vesicle-vesicle and vesicle-boundary 
% intearctions.  Either 'explicit' or 'implicit'
% Boundary-boundary interactions are always handled 
% implicitly.  The preconditioner for the boundary-boundary
% interactions is exact
options.near = true;        % near-singular integration
options.fmm = true;
options.fmmDLP = true;
% use FMM to compute single-layer potential
options.confined = true;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.usePlot = true;     % View vesicles during run
options.track = false;      % Tracker points on vesicle
options.quiver = false;      % draw velocity vectors on vesicle
options.axis = [-8 0 -3 3];
options.axis = [-12 12 -3 3];
%options.axis = [-8 2 -3 3];
% Axis for plots if usePlot = true
options.saveData = true;
% Save vesicle information and create a log file
options.logFile = 'output/diffuser.log';
% Name of log file
options.dataFile = 'output/diffuserData.bin';
% Name of binary data file for storing vesicle information
options.verbose = true;
% Decides how much info is written to console
options.profile = false;    % Profile code
options.collision = false;   % Collision detection
options.tracers = false;      % trace passive particles in flow
options.pressure = false;     % compute pressure at certain fixed points
% For now, there will be a bug when using near-singular integration to
% compute the pressure at a point interior to a vesicle
options.timeAdap = true;
options.orderGL = 3;
% number of Gauss-Lobatto points to use in the substep.  Letting
% orderGL=1 corresponds to not doing any intermediate time step sizes
options.nsdc = 1;
options.periodic = true;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
Xwalls = oc.initConfig(prams.Nbd,options.farField);
% Build solid walls

scale = 0.1;
%[X,prams.nv] = oc.initConfig(prams.N,'volFrac',0.1,Xwalls,...
%    options.fmm,'scale',scale,'reducedArea',0.65);
[X,prams.nv] = oc.initConfig(prams.N,'volFrac',0.02,Xwalls,...
    options.fmm,'scale',scale,'reducedArea',0.65);
ind = find(max(X(1:end/2,:)) <= 8);
X = X(:,ind); 
%X = X(:,1);
prams.nv = size(X,2);



%addpath output
%[posx,posy,ten,wallx,wally,ea,el,time,n,nv] = loadFile('output/diffuser2StuckData.bin');
%X = [posx(:,:,end);posy(:,:,end)];
%prams.nv = size(X,2);
%prams.m = round(prams.T/(time(end) - time(end-1)));


Xfinal = Ves2D(X,Xwalls,prams,options);
% Run vesicle code


