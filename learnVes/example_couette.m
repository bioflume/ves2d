% COUETTE FLOW EXAMPLE

clear; clc;
addpath ../src/
addpath ../examples/

% PARAMETERS
%-------------------------------------------------------------------------
prams.bgFlow = 'couette'; % 'rotation' or 'couette' (w/ solid boundaries)
prams.speed = 100; % 70 for rotation, 100 for couette 
prams.Th = 1.5; % time horizon

prams.N = 32; % num. points per vesicle
prams.nv = 81; % num. of vesicles
prams.kappa = 2e-1; % bending stiffness
prams.dt = 5E-5; % time step size

prams.Nbd = 512; % num. points per solid wall
prams.nvbd = 2; % num. of solid walls

% OPTIONS
%-------------------------------------------------------------------------
options.exactSolve = false; % either exact solve or MLARM
options.plotOnFly = false;
options.jiggle = true;
prams.fmm = true; % use FMM for ves2ves
prams.fmmDLP = true; % use FMM for ves2walls
options.nCompRelax = 32; % number of PCA components for relaxation problem
options.fileName = ['./output/alternTstep_MLARMcouette_nves' num2str(prams.nv) ...
    '_Nves' num2str(prams.N) '_Kappa' num2str(prams.kappa) ...
    '_Dt' num2str(prams.dt) '_Speed' num2str(prams.speed) '.mat'];

% VESICLES and WALLS:
% -------------------------------------------------------------------------
% There are several ICs already initialized:
% VF35_81VesIC.mat
% VF30_70VesIC.mat
% ellipseVF20pN16
% ellipseVF10pN16
% ellipse15vesN24closeIC (VF = 4.4%, closely initialized)
% ellipse10vesN24closeIC (VF = 2.2%, closely initialized)
initialVesFile = 'VF35_81VesIC'; % from another run % probVF20IC: fails in 2it
volFrac = 0.3; % leave initialVesFile empty if you want to initialize new IC

% initialize vesicles and walls
[X,Xwalls,prams] = initializeVesiclesAndWalls(initialVesFile,...
    prams,volFrac); % see below for the code
% -------------------------------------------------------------------------

% RUNNING A SIMULATION:
Xhist = MLARM(X,Xwalls,prams,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Xwalls,prams] = initializeVesiclesAndWalls(initialVesFile,...
    prams,volFrac)
oc = curve;
thet = (0:prams.Nbd-1)'*2*pi/prams.Nbd;
Xwalls = [ [2.2*cos(thet); 2.2*sin(thet)] [cos(-thet); sin(-thet)] ];


if ~isempty(initialVesFile)
  disp('loading an IC file')
  load(initialVesFile)   
  if numel(X(:,1))/2 ~= prams.N
    X = [interpft(X(1:end/2,:),prams.N);interpft(X(end/2+1:end,:),prams.N)];
  end 
else
  X0 = oc.initConfig(N,'ellipse');
  [~,~,len] = oc.geomProp(X0);
  X0 = X0./len;    
  op = poten(prams.N);
  X = oc.fillCouetteArea(X0,Xwalls,volFrac,false,op);
end % if initType

% number of vesicles
[~,prams.area0,prams.len0] = oc.geomProp(X);
prams.nv = numel(prams.area0);

end % initializeVesiclesAndWalls

