clear; clc;
fprintf('Simple elliptical vesicle in a free Couette flow.\n');
restoredefaultpath;
rng('shuffle');
% volume fraction
volFrac = 0.2;

% Save vesicle information and create a log file
prams.runName = 'freeCouetteVF20per_ideal';

options.farField = 'freeCouette'; % background velocity
options.farFieldSpeed = 1; % scaling of background velocity
prams.nv = 1;    % number of vesicles

options.saveVinf = ~true;

% Spatial-resolution
prams.N = 64;    % points per vesicle

% physical properties of vesicles
prams.kappa = 1e-2;   % bending coefficient
prams.viscCont = 1;   % viscosity contrast

% parameters for numerics
options.fmm = false;  % fmm for single-layer potentials
options.fmmDLP = false; % fmm for double-layer potentials

prams.gmresTol = 1e-10;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

% Low-Resolution Correction Algorithms
options.antiAlias = true; % upsampling for anti-aliasing
options.repulsion = false; % repulsion
options.reparameterization = true; % reparametrization of membrane
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA

% Temporal Resolution (parameters for new implementation)
prams.T = 1;   % time horizon (or 100*lenScale/options.farFieldSpeed) etc.
options.saveVtotal = true;

options.timeAdap = ~true; % time adaptivity
options.nsdc = 0;  % number of sdc corrections
prams.areaLenTol = 5e-3; % tolerance for errors in area-length for time-stepping


options.usePlot = ~true; % Plot on-the-fly
options.track = false; % trackers on membrane


% Name of log file for saving messages
options.logFile = ['output/' prams.runName '.log'];
% Name of binary data file for storing vesicle information
options.dataFile = ['output/' prams.runName 'Data.bin'];


% Set options and parameters that the user doesn't
% Also add src to path
options.axis = [-21 21 -21 21];
[options,prams] = initVes2D(options,prams);


% time scale for simulation
% # of time steps (if constant time stepping) or Dt for the first step
% (adaptive)
prams.m = prams.T/0.001;
prams.dtMax = 5*1/options.farFieldSpeed; % maximum allowable time-step size
prams.dtMin = 1e-4*1/options.farFieldSpeed; % minimum allowable time-step size

om = monitor([],options,prams);
tt = tstep(options,prams,om);
% Initial configuration of vesicles
oc = curve;
if 1 % based on the volume fraction
  X = oc.fillCouetteArea(volFrac,prams.N,tt,[10 20]);
  prams.nv = numel(X(1,:));
  prams.viscCont = prams.viscCont*ones(prams.nv,1);

else % a single vesicle
  load relaxed64.mat
  X = [interpft(X(1:end/2),prams.N);interpft(X(end/2+1:end),prams.N)];
  
  X = 3.57*[X(1:end/2)-mean(interpft(X(1:end/2),256));...
      X(end/2+1:end)-mean(interpft(X(end/2+1:end),256))];
  IA = oc.getIncAngle(X);  
  Xnew = zeros(numel(X),prams.nv);

  for i = 1 : prams.nv
    IArand = rand*2*pi;
    % Decenter and rotate
    cr = 12+6*rand; ctheta = 2*pi*rand;
    cx = cr*cos(ctheta); cy = cr*sin(ctheta);
    
    Xnew(1:end/2,i) = cos(IArand-IA)*X(1:end/2)-sin(IArand-IA)*X(end/2+1:end)+cx;
    Xnew(end/2+1:end,i) = sin(IArand-IA)*X(1:end/2)+cos(IArand-IA)*X(end/2+1:end)+cy;
  end
  X = Xnew;
  prams.viscCont = prams.viscCont*ones(nv,1);
end
om = monitor(X,options,prams);
tt.om = om;
Xfinal = Ves2D(X,[],[],[],prams,options,tt,[]);
% Run vesicle code

