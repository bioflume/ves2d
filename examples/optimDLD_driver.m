% Driver for unconstrained nonlinear optimization of post shape for DLD
clear; clc;
addpath ../src/
addpath ./output/

% DLD device with a 
prams.period = 5;

% Vesicle properties (first vesicle is wanted to displace, the second one 
% wanted to zig-zag), we know that higher VCs and lower kappas zig-zag, so
% choose parameters based on that
prams.VCs = [5;10]; % viscosity contrasts of the cells
% If vesicles with different properties will be in the same simulation,
% then we cannot run it for different kappas. That requires separating
% cells with different kappas or modifying the code such that it allows
% different kappas in one simulation
prams.Ufar = 1.2; % maximum velocity 
prams.kappas = 1e-2;  % bending stiffness of the cells 
 
% optimize for pillar shape only, 
% Dx and Dy come out of the bounds for pillars.
nParams = 16;

% to get the function evaluation faster
options.useFMM = true; % FMM
options.useHODLR = false; % HODLR to compress&factorize walls precond.
options.repulsion = true; % repulsion
options.HODLRforW2W = false; % HODLR to compress wall2wall interactions

% FILE NAMES
prams.runName = ['sDLD'];
prams.saveFolder = ['./output/OptimRuns/VCs' ...
    num2str(prams.VCs(1)) 'a' num2str(prams.VCs(2)) '_Per' ...
    num2str(prams.period) '/'];

prams.fileData = [prams.saveFolder prams.runName '_OptimData.bin'];
prams.logFile = [prams.saveFolder prams.runName '_Info.log'];
prams.fileParams = [prams.saveFolder prams.runName '_params.bin'];

% DISCRETIZATION
% prams.Next = 1024; % # of points for the exterior wall
prams.Next = 1280; % FOR PUT-BACK (3 rows x 2 columns model)
prams.Nint = 64; % # of points per an interior wall
prams.Nves = 32;

% we assume that a pillar lives in a box which has equal dimensions in x-y
% wbox = Dx+Dpostx (Dx: horizontal spacing, Dpostx: post's max. horiz. dim)
% wbox = Dy+Dposty (Dy: vertical spacing, Dposty: post's max. vertc. dim)

% for a given post shape, we compute Dx and Dy using the above equations.
% if Dx or Dy is less than 0.25 (dim. of a cell), then we penalize that
% post shape. So Dx, Dy depend on post shape and are not parameters for
% optimization
prams.wbox = 2.5; % based on our initial guess (Dx=Dy=1 and Dpost = 1.5)
prams.spaLowLim = 0.7; % if spacing is less than this, then penalize the shape
% 0.25 is a dimension of a cell.

% Initial guess 
prams.Dpost = 1.5; % initial post diameter
prams.Dx = 1; % initial gap size in horizontal direction
prams.Dy = 1; % initial gap size in vertical direction

if 1
% initial guess is a circle (D = 1.5)
% 8 control points for B-splines are used
points = [[1 0];0.5*[sqrt(2) sqrt(2)];[0 1];...
    0.5*[-sqrt(2) sqrt(2)];[-1 0];0.5*[-sqrt(2) -sqrt(2)];...
    [0 -1];0.5*[sqrt(2) -sqrt(2)]]*1.5/1.7581;
stdDev = 0.15; % standard deviation for the multivar. normal dist.    
end

% load warmStartData
% points = optimPoints;
% stdDev = optimStd;

%load newOptPosts
%points = bsplinePoints(:,bestIds(1));
%stdDev = 0.10;

% OPTIMIZATION
% build the optimization class
% ooptim = optimDLD_setUp(prams,options);
ooptim = optimDLD2Ves_setUp(prams,options);

% keep number of iterations
global niter;
niter = 0;

% Optimization options
opts = cmaes('defaults');
opts.TolX = 1e-3; opts.TolFun = 1e-2; 
opts.SaveVariables = false;

% if parallel CMAES then the objective function admits NxM input (N: number
% of parameters, M: number of parallel evaluations)
opts.EvalParallel = 'yes';
opts.EvalInitialX = 'yes';
opts.StopOnStagnation = 'on';
%nParEvals = 4+floor(3*log(nParams)); % number of parallel evaluations
nParEvals = 32;
opts.PopSize = nParEvals;

% Note: Initial guess is just one set of parameters, it evaluates loss
% function and in the next iterations it samples nParEvals sets of
% parameters and evaluates them in parallel.

% the initial parameters
xstart = points(:);

% run the optimizer
%parprofile = parcluster('local');
parprofile = parallel.cluster.MJS;
delete(findJob(parprofile));

ooptim.parCluster = parprofile;
xend = cmaes('DLD_objFunc',xstart,stdDev,opts,ooptim);


