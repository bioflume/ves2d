% Driver for unconstrained nonlinear optimization of post shape for DLD
clear; clc;
addpath ../src/
addpath ./output/
addpath ../examples/

% DLD device with a 
prams.period = 6;

% Vesicle properties (first vesicle is wanted to displace, the second one 
% wanted to zig-zag), we know that higher VCs and lower kappas zig-zag, so
% choose parameters based on that
prams.VCs = 1; % viscosity contrasts of the cells
prams.volFrac = 0.25; % density of vesicle suspension
prams.RAparticle = 0.9; % reduced area of the particle
prams.lenParticle = 1.0; % arc-length of the particle, 1.5 or 2.5

prams.Ufar = 80; % maximum velocity 
prams.kappas = 1;  % bending stiffness of the cells 
 
% optimize for pillar shape only, 
% Dx and Dy come out of the bounds for pillars.
nParams = 16;

% to get the function evaluation faster
options.useFMM = true; % FMM
options.useHODLR = false; % HODLR to compress&factorize walls precond.
options.repulsion = false; % repulsion
options.HODLRforW2W = false; % HODLR to compress wall2wall interactions

% FILE NAMES
prams.runName = 'lowFid';
prams.saveFolder = ['/workspace/gokberk/denseOptimRuns/len1_RAp' ...
    num2str(prams.RAparticle*100) '_VFp' num2str(100*prams.volFrac) '_Per' ...
    num2str(prams.period) '/'];

prams.fileData = [prams.saveFolder prams.runName '_OptimData'];
prams.logFile = [prams.saveFolder prams.runName '_Info.log'];

% DISCRETIZATION
% prams.Next = 1024; % # of points for the exterior wall
prams.Next = 1280; % FOR PUT-BACK (3 rows x 2 columns model)
prams.Nint = 48; % # of points per an interior wall
prams.Nves = 32;

% we assume that a pillar lives in a box which has equal dimensions in x-y
% wbox = Dx+Dpostx (Dx: horizontal spacing, Dpostx: post's max. horiz. dim)
% wbox = Dy+Dposty (Dy: vertical spacing, Dposty: post's max. vertc. dim)

% for a given post shape, we compute Dx and Dy using the above equations.
% if Dx or Dy is less than 0.25 (dim. of a cell), then we penalize that
% post shape. So Dx, Dy depend on post shape and are not parameters for
% optimization
prams.wbox = 36; % based on our initial guess (Dx=Dy=10 and Dpost = 15)
prams.spaLowLim = 14; % if spacing is less than this, then penalize the shape
% 0.25 is a dimension of a cell.

% Initial guess 
prams.Dpost = 20; % initial post diameter
prams.Dx = 16; % initial gap size in horizontal direction
prams.Dy = 16; % initial gap size in vertical direction
prams.DLDscale = 0.05; % this multiplies all the sizes

if 1
% initial guess is a circle (D = 1.5)
% 8 control points for B-splines are used
points = [[1 0];0.5*[sqrt(2) sqrt(2)];[0 1];...
    0.5*[-sqrt(2) sqrt(2)];[-1 0];0.5*[-sqrt(2) -sqrt(2)];...
    [0 -1];0.5*[sqrt(2) -sqrt(2)]]*(prams.Dpost*prams.DLDscale)/1.7581;
stdDev = 2*(prams.DLDscale); % standard deviation for the multivar. normal dist.    
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
ooptim = optimDLDdense_setUp(prams,options);

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
nParEvals = 40;
opts.PopSize = nParEvals;

% Note: Initial guess is just one set of parameters, it evaluates loss
% function and in the next iterations it samples nParEvals sets of
% parameters and evaluates them in parallel.

% the initial parameters
xstart = points(:);

% run the optimizer
%parprofile = parcluster('local');
parprofile = parallel.cluster.MJS('Name','luke-mjs');
delete(findJob(parprofile));

ooptim.parCluster = parprofile;
xend = cmaes('DLD_objFunc',xstart,stdDev,opts,ooptim);


