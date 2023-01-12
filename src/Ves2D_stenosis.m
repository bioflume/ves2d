function Xfinal = Ves2D_stenosis(X,Xwalls,XwallsInt,XwallsExt,prams,options,tt,...
    Xtra,pressTar)
% Ves2D does time stepping on the intial configuration X
% with parameters and options defined in prams and options.
% Also can pass a set of initial tracer locations (Xtra) and
% target points where one wants the pressure and stress (pressTar)

% build an object of class tstep with required options and parameters
% and an object for outputting
vesWidth = max(X(1:end/2))-min(X(1:end/2));

if isempty(tt)
  om = monitor(X,options,prams);
  tt = tstep(options,prams,om);
else
  om = tt.om;
end

oc = curve;
if nargin == 7
  Xtra = [];       % Xtra is positions of tracers for visualization
  pressTar = [];   % points to compute pressure for postprocessing
elseif nargin == 8
  pressTar = [];
end

% GLOBAL and LOCAL COUNTERS
global matvecs ;  % number of matvecs
global derivs  ;  % number of times we compute differential
                  % operators for preconditioning
global fmms       % number of fmm calls

global repuls  ;  % number of times we activate repulsion

global ncollVes; % number of tsteps vesicle-vesicle collision occurs

global ncollWal; % number of tsteps vesicle-wall collision occurs

matvecs  = 0; % counter for the total number of time steps
derivs   = 0; % counter for the total number of times derivative
              % operators are computed
fmms     = 0; % counter for the total number of FMM calls
repuls   = 0; % counter for the number of times we activate repulsion

ncollVes = 0; % counter for the number of tsteps (rejected) ves-ves
              % collision occurs

ncollWal = 0; % counter for the number of tsteps (rejected) ves-wall
              % collision occurs

nreject = 0; % count for the number of rejected time steps
naccept = 0; % count for the number of accepted time steps
countGMRES = 0; % count for the number of total GMRES iterations

% For small DLD we put vesicle back to between first two pillars, and let
% vesicle move over 2 pillars. So we want to do it nPutBackMax = periodN+1
% times at most.
if options.putBackDLD
  if isempty(prams.nPutBackMax)
    nPutBackMax = ceil(prams.periodN)+1;
  else
    nPutBackMax = prams.nPutBackMax;
  end
  nPutBacks = zeros(prams.nv,1);
end

% GET # of VESICLES, WALLS and THEIR DISCRETIZATIONS
N = prams.N; % Number of points per vesicle
nv = prams.nv; % Number of vesicles

Nbd = prams.Nbd; % number of points per solid wall
nvbd = prams.nvbd; % number of solid wall components
NbdInt = prams.NbdInt;
NbdExt = prams.NbdExt;
nvbdInt = prams.nvbdInt;
nvbdExt = prams.nvbdExt;
if options.diffDiscWalls
  % then we have two sets of walls with different disc.
  % so add the number of exterior and interior walls to find total walls
  nvbd = nvbdInt+nvbdExt;
end

% Turn on profiler
if options.profile
  profile off; profile on -timer real;
%  profile off; profile on
end
% Start a timer
tstart = tic;
% get current time in case tic and toc give funny anwers
c = clock;

% PRINT INFORMATION ABOUT THE RUN AND KEEP IT IN THE LOG-FILE, TOO
message = ['The simulation has been started on '];
om.writeMessage(message,'%s\n')
message = ['Initial Month is ' num2str(c(2),'%d')];
om.writeMessage(message,'%s\n')
message = ['Initial Day is ' num2str(c(3),'%d')];
om.writeMessage(message,'%s\n')
message = ['Initial Hour is ' num2str(c(4),'%d')];
om.writeMessage(message,'%s\n')
message = ['Initial Minute is ' num2str(c(5),'%d')];
om.writeMessage(message,'%s\n')
message = ['Initial Second is ' num2str(round(c(6)),'%d')];
om.writeMessage(message,'%s\n')
message = ' ';
om.writeMessage(message,'%s\n')
% Also print out some information about the run



% END OF PRINTING INFORMATION ABOUT THE RUN


% Find the repulsion strength scale
if options.repulsion
  vesicle1 = capsules(X(:,1),[],[],prams.kappa,prams.viscCont(1),...
      options.antiAlias);
  prams.repStrength = vesicle1.repulStrengthScale(prams.minDist,tt,...
      options.farFieldSpeed);

  message = ['Repulsion strength is ' num2str(prams.repStrength,'%d')];
  om.writeMessage(message,'%s\n');
  message = ' ';
  om.writeMessage(message,'%s\n');
end

% INITIALIZE
% initial velocity on the vesicles is set to the background velocity.
% If flow is unbounded, there is no density function eta.  If bounded,
% compute a structure for the solid walls

sigma = zeros(N,nv);
if ~options.confined
  u = tt.farField(X,[]);
  walls = [];
  wallsInt = [];
  wallsExt = [];
else
  u = zeros(2*N,nv);
  if ~options.diffDiscWalls
    walls = tt.initialConfined(prams,Xwalls,[],[]);
    wallsInt = [];
    wallsExt = [];
  else
    [~,wallsInt,wallsExt] = tt.initialConfined(prams,[],...
        XwallsInt,XwallsExt);
    walls = [];
  end
end

% []store keeps information for all vesicles, this is used for streaming
% and freezing in DLD usually, for other cases []store = [] basically
% This streaming is implemented for DLD only, but can be used for other
% geometries after a couple of easy adjustments

% keep the GLOBAL ids of remaning vesicles (needed for freezing and streaming)
rem_ids = [1:nv]';

if options.streaming
  % assign viscosity contrasts randomly for the cells to be streamed
  nves2 = ceil(prams.totnv*prams.ves2Concent); % number of 2nd type of vesicles
  ves2Ids = randperm(prams.totnv,nves2); % global ids of 2nd type vesicles
  % keep the viscosity contrasts of all cells to be streamed
  allVesViscConts = ones(prams.totnv,1)*prams.vesViscCont;
  allVesViscConts(ves2Ids) = prams.ves2ViscCont;
  % assign viscosity contrasts of initialized vesicles
  prams.viscCont = allVesViscConts(1:prams.nv);

  % keeps all vesicles (frozen and remaining in the simulation)
  Xstore = zeros(2*N,prams.totnv);
  Xstore(:,1:nv) = X;
  sigStore = zeros(N,prams.totnv);
  uStore = zeros(2*N,prams.totnv);

  % also keep the ids of vesicles generated so far
  genSoFar_ids = [1:nv]';

  % the following is needed to keep track of saved files, saving is
  % different when we stream.
  % keep the number of instants saved in the simulation
  nInstant = 0;

else
  Xstore = X;
  sigStore = sigma;
  uStore = zeros(2*N,prams.nv);
  genSoFar_ids = [];
  allVesViscConts = prams.viscCont;
  nInstant = [];
end

% For higher-order methods (only 2nd order for now), need to initially
% take smaller time steps so that we have two initial conditions
if ~options.diffDiscWalls

  [X,sigma,u,eta,RS,Xtra] = tt.firstSteps(options,prams,...
      X,walls,Xtra,pressTar);
  etaInt = [];
  etaExt = [];
else
  [Xstore,sigStore,uStore,etaInt,etaExt,RS,Xtra] = ...
      tt.firstStepsDLD(prams,Xstore,sigStore,uStore,wallsInt,...
      wallsExt,rem_ids,Xtra,options.streaming,allVesViscConts);
  X = Xstore(:,rem_ids);
  sigma = sigStore(:,rem_ids);
  u = uStore(:,rem_ids);
  eta = [];
end

% initial time.  firstSteps took the first time steps so that there is
% enough data to use the time stepping order that is desired
time = (options.order-1)*tt.dt;

% if we will save velocity due to fixed tension, compute fixed tension
if options.saveVtotal
  sigFixed = computeSigBasedOnRelaxed(N,prams.kappa,prams.viscCont(1));
  sigFixed = repmat(sigFixed,1,nv);
end

accept = true;
% MAIN TIME STEPPING LOOP
while time < prams.T - 1e-10
  % make sure we land exactly on the time horizon
  if time+tt.dt > prams.T
    tt.dt = prams.T - time;
  end

  dt = tt.dt;
  tt.currentTime = time;
  time = time + dt; % current time

  %------------------------------------------------------------------------
  % 1) REPARAMETERIZE
  %------------------------------------------------------------------------
  if options.reparameterization && accept

    % reparameterize so that Fourier coeff of X decays at a desired rate
    % If there is a collision after reparametrization, do not do that
    [Xreparam,niter] = oc.reparametrize(X,u*dt,6,prams.maxReparamIter);
    message = ['Max. number of reparameterization iterations ' ...
        num2str(max(niter))];
    om.writeMessage(message,'%s\n')


    % align center and inclination angle of the reparametrized vesicles
    % with the those of the ones given by the time stepping
    if options.alignCenterAngle
      [Xreparam] = oc.alignCenterAngle(X,Xreparam);
    end
    % check if reparametrization leads to collision, if so undo
    if ~options.diffDiscWalls
      X = oc.checkCollisionForLRCA(Xreparam,X,walls,options.fmm,tt.op,om,1);
    else
      X = oc.checkCollisionForLRCA(Xreparam,X,wallsInt,options.fmm,tt.op,om,1);
    end
  end

  if options.equiDistArcLength && accept
    % reparameterize such that points are distributed equally in arc-length
    Xreparam = X;
    for it = 1 : 10
      Xreparam = oc.redistributeArcLength(Xreparam);
    end
    niter = 10;
    message = ['Max. number of reparameterization iterations ' ...
        num2str(max(niter))];
    om.writeMessage(message,'%s\n')

    % align center and inclination angle of the reparametrized vesicles
    % with the those of the ones given by the time stepping
    if options.alignCenterAngle
      [Xreparam] = oc.alignCenterAngle(X,Xreparam);
    end
    % check if reparametrization leads to collision, if so undo
    if ~options.diffDiscWalls
      X = oc.checkCollisionForLRCA(Xreparam,X,walls,options.fmm,tt.op,om,1);
    else
      X = oc.checkCollisionForLRCA(Xreparam,X,wallsInt,options.fmm,tt.op,om,1);
    end
  end
  %------------------------------------------------------------------------

  %------------------------------------------------------------------------
  % 2) TIME STEPPING
  %------------------------------------------------------------------------
  if options.saveVtotal
%     [~,~,~,~] = computeTotalVelocityExplicit(X,prams,options,tt);
    [~,~,~,~] = computeTotalVelocityPreComp(X,sigFixed,prams,options,tt);
  end

  [X,sigma,u,eta,etaInt,etaExt,RS,iter,accept,dtScale,res,iflag,...
      collUprate] = tt.timeStepGL(X,sigma,u,eta,etaInt,etaExt,RS,...
      prams.kappa,prams.viscCont,walls,wallsInt,wallsExt);
  countGMRES = countGMRES + iter;

  %------------------------------------------------------------------------

  %------------------------------------------------------------------------
  % 3) AREA-LENGTH CORRECTION
  %------------------------------------------------------------------------
  if options.correctShape && accept
    [Xcorr] = oc.correctAreaAndLength(X,prams.areaLenTol,om);

    %Align center and angle
    if options.alignCenterAngle
      [Xcorr] = oc.alignCenterAngle(X,Xcorr);
    end

    % Collision detection after correction (if there is any crossing, do
    % not correct the shape)

    if ~options.diffDiscWalls
      X = oc.checkCollisionForLRCA(Xcorr,X,walls,options.fmm,tt.op,om,collUprate);
    else
      X = oc.checkCollisionForLRCA(Xcorr,X,wallsInt,options.fmm,tt.op,om,collUprate);
    end
  end
  % correct length and area
  %------------------------------------------------------------------------


  %------------------------------------------------------------------------
  % 4) CHECK IF ERROR, DT ARE ACCEPTABLE TO GO ON, and
  % IF ALL THE VESICLES HAVE BEEN STREAMED AND COMPLETED SIMULATION.
  %------------------------------------------------------------------------
  if ~accept
    time = time - dt;
  end
  % go back to old time

  if accept

    % if confined Poiseuille flow and putBackOrigin
    if options.putBackOrigin
      X = [X(1:end/2)-mean(X(1:end/2));X(end/2+1:end)];
    end


    % if DLD case, then check whether we freeze and stream in this iteration
    % then
    % update the positions, tension, and velocity field of
    % the vesicles, and the density function and rotlet and
    % stokeslets
    if options.freezing && options.diffDiscWalls % this is for DLD example
      [X,sigma,u,Xstore,sigStore,uStore,rem_ids,genSoFar_ids,prams,om,nInstant] = ...
          oc.freezeAndStream(om,tt,Xstore,sigStore,uStore,X,sigma,u,prams,...
          XwallsExt,rem_ids,genSoFar_ids,options.streaming,...
          allVesViscConts,time,nInstant);
    elseif options.putBackDLD
      % if DLD and putBackDLD (usually for small DLD, move back to the 1st
      % column), if vesicle has already zig-zagged or nPutBacks exceeds
      % nPutBackMax, freeze it.
      [X,sigma,u,Xstore,sigStore,uStore,rem_ids,nPutBacks,prams,om] = ...
          oc.freezeAndPutBack(om,Xstore,sigStore,uStore,X,sigma,u,prams,...
          rem_ids,allVesViscConts,nPutBacks,nPutBackMax);
    else
      Xstore = X;
      sigStore = sigma;
      uStore = u;
    end


    % TERMINATION CRITERIA
    % ---------------------------------------------------------------------
    % check if we have violated the error in area or length
    % also plot and save the current solution to dat file.
    % Print information to log file and console
    terminate = om.outputInfo(Xstore,sigStore,uStore,Xwalls,XwallsInt,...
        XwallsExt,Xtra,time,iter,dtScale,res,rem_ids,dt,allVesViscConts,iflag);


    % ---------------------------------------------------------------------
    naccept = naccept + 1;

  else
    nreject = nreject + 1;
    terminate = false;
    if dt <= prams.dtMin
      message = ['Time step shrinks too much, it is not feasible to complete it'];
      om.writeMessage(message,'%s\n')
      message = ['I AM STOPPING'];
      om.writeMessage(message,'%s\n')
      terminate = true;
    end

  end % if accept
  % write a final summary with number of matvecs, accepts/rejects,
  % and the CPU time
  % Check the location
  if max(X(1:end/2)) >= max(Xwalls(1:end/2))-vesWidth
      terminate = true;
      message = ['Reached the wall on the right'];
      om.writeMessage(message,'%s\n')
  end
  if terminate
    Xfinal = X;
    totTime = toc(tstart);
    message = ['Final Month is ' num2str(c(2),'%d')];
    om.writeMessage(message,'%s\n')
    message = ['Final Day is ' num2str(c(3),'%d')];
    om.writeMessage(message,'%s\n')
    message = ['Final Hour is ' num2str(c(4),'%d')];
    om.writeMessage(message,'%s\n')
    message = ['Final Minute is ' num2str(c(5),'%d')];
    om.writeMessage(message,'%s\n')
    message = ['Final Second is ' num2str(round(c(6)),'%d')];
    om.writeMessage(message,'%s\n')
    message = ' ';
    om.writeMessage(message,'%s\n')
    om.summary(matvecs,derivs,fmms,repuls,ncollVes,ncollWal,...
        naccept,nreject,countGMRES,totTime);
    return;
  end

end
% end of main

% Total time doing time stepping
totTime = toc(tstart);

% get current time in case tic and toc give funny anwers
c = clock;
message = ['Final Month is ' num2str(c(2),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Day is ' num2str(c(3),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Hour is ' num2str(c(4),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Minute is ' num2str(c(5),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Second is ' num2str(round(c(6)),'%d')];
om.writeMessage(message,'%s\n')
message = ' ';
om.writeMessage(message,'%s\n')

% Save profile
if options.profile
  profile off;
  p = profile('info');
  filename = [options.logFile(1:end-4) 'Profile'];
  save([filename '.mat'],'p');
%  profview
  profsave(profile('info'),filename);
end

% write a final summary with number of matvecs, accepts/rejects,
% and the CPU time
om.summary(matvecs,derivs,fmms,repuls,ncollVes,ncollWal,...
    naccept,nreject,countGMRES,totTime);

% final configuation
Xfinal = X;


message = 'SIMULATION SUCCESSFULLY COMPLETED';
om.writeMessage(message,'%s\n')
om.writeStars
