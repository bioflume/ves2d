clear all; %clc

disp('One elliptical vesicles in a couette apparatus.');
disp('First-order semi-implicit time stepping.');
disp('Implicit vesicle-boundary interactions.');

format short e
prams.N = 64;                  % points per vesicle
prams.nv = 42;                  % number of vesicles
prams.T = 10;                % time horizon
prams.m = 286*4;                  % number of time steps
prams.Nbd = 128;                % number of points on solid wall
prams.Nbdcoarse = 128;          % number of points on solid wall
prams.nvbd = 2;                % number of components to solid walls
pams.kappa = 1e-1;            % bending coefficient
prams.gmresTol = 1e-6;         % GMRES tolerance
prams.gmresMaxIter = 100;
% maximum number of iterations at any one GMRES call
prams.errorTol = 1e-1;
% Maximum relative error in area and length before the simulation
% is stopped
prams.rtolArea = 1e-2;
prams.rtolLength = 1e-2;


options.farField = 'couette';    % Constricted domain
% background velocity or solid walls (which comes with a default
% boundary condition)
options.order = 2;          % time stepping order
option.inextens = 'method1'; 
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
options.usePlot = false;     % View vesicles during run
options.axis = [-20 20 -20 20];
% Axis for plots if usePlot = true
options.saveData = true;
% Save vesicle information and create a log file
options.logFile = 'output/couette42.log';
% Name of log file
options.dataFile = 'output/couette42Data.bin';
% Name of binary data file for storing vesicle information
options.verbose = true;
% Decides how much info is written to console
options.timeAdap = false;
options.orderGL = 2;
% number of Gauss-Lobatto points to use in the substep.
options.nsdc = 0;
options.correctShape = false;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
if options.confined
  Xwalls = oc.initConfig(prams.Nbd,...
      options.farField,'nv',prams.nvbd,'center',[[0;0] [-5;0]]);
end
% Build solid walls if they exist

load couette42Init.dat;
x = couette42Init(:,1);
y = couette42Init(:,2);
x = reshape(x,64,42);
y = reshape(y,64,42);
%x = interpft(x,128);
%y = interpft(y,128);
X = [x(:,1:prams.nv);y(:,1:prams.nv)];


Xfinal = Ves2D(X,Xwalls,prams,options,[],[]);
% Run vesicle code


