classdef optimDLDdense_setUp
% This class defiens the functions required to optimize the shape of the posts
% used in the (D)eterministic (L)ateral (D)isplacement device to separate
% particles. It optimizes the shape such that one of the vesicles having 
% different capillary numbers (Cas) or viscosity contrasts (VCs) zig-zag and
% the other one displaces. 

properties
runName    % a common phrase in all related files
Next       % number of points to discretize the outer wall
Nint       % number of points to discretize the posts
Nves       % number of points per vesicle
Ufar       % maximum velocities for two vesicles (different Cas)
kappas     % bending stiffness
VCs        % viscosity contrasts for two vesicles
useFMM     % whether we use FMM or not
useHODLR   % whether we use HODLR to compress&factorize the matrices for wall2wall ints.
HODLRforW2W% whether we use HODLR to compress wall2wall interactions
repulsion  % whether we use repulsion or not
fileData   % file to save post shapes, spacings and objective function
logFile    % log file of the optimization run

saveFolder % folder in which we save
RAparticle % RA of circular particle
lenParticle% length of circular particle
DLDscale   % scaling of DLD device
period     % period
theta      % tilt angle
initFlux   % initial Dy*Ufar, scale Ufar as Dy changes
volFrac 

% size of the box containing pillar
% i.e. wbox = Dpostx+Dx = Dposty+Dy
wbox

% lower limit for spacing, if Dx or Dy is less than this, then penalize the
% shape
spaLowLim

% parallel profile
parCluster

end % end properties

methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = optimDLDdense_setUp(prams,options)
    
o.runName = prams.runName;
o.period = prams.period;
o.theta = atan(1/prams.period);
o.Next = prams.Next;
o.Nint = prams.Nint;
o.Nves = prams.Nves;
o.Ufar = prams.Ufar;
o.VCs = prams.VCs;
o.kappas = prams.kappas;
o.RAparticle = prams.RAparticle;
o.lenParticle = prams.lenParticle;
o.initFlux = prams.Ufar*prams.Dy*prams.DLDscale;
o.volFrac = prams.volFrac;

o.useFMM = options.useFMM;
o.useHODLR = options.useHODLR;
o.HODLRforW2W = options.HODLRforW2W;
o.repulsion = options.repulsion;

o.fileData = prams.fileData;
o.logFile = prams.logFile;
o.saveFolder = prams.saveFolder;
o.DLDscale = prams.DLDscale;
o.wbox = prams.wbox;
o.spaLowLim = prams.spaLowLim;

% initialize the files
Nint = o.Nint; wbox = o.wbox; initFlux = o.initFlux;
RApart = o.RAparticle; lenPart = o.lenParticle;
save([o.fileData '_SETUP'],'Nint','wbox','initFlux','RApart','lenPart')

message = ['Optimization for post shape and wbox is ' num2str(o.wbox) '.'];
paramSize = 16;
o.saveLog(o.logFile,message,0);

message = ['Initial flux is ' num2str(o.initFlux)];
o.saveLog(o.logFile,message,0);
message = ['Bending stiffness is ' num2str(o.kappas)];
o.saveLog(o.logFile,message,1);

end % end shapeDLDoptimPar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function Fvals = DLD_objFunc(xorig,o)
% Start timing one set of iterations
tIter = tic;

% Get the number of parameters and parallel evaluations
% Input might be a multiple sets of parameters but also a single set of 
% parameters. The latter occurs when there are NaN's in the multiple sets
% of parameters. In those cases, it resamples and evaluates the objective
% function one by one for each set until it has not NaN PopSize function
% evaluations.
[~,nParEvals] = size(xorig);

% Increase the number of iterations whenever it is called
global niter;

% Run Names for Ves2D
localIters = zeros(nParEvals,1);
runNames = cell(nParEvals,1);
for j = 1 : nParEvals
  niter = niter + 1;
  runNames{j} = [o.runName '_Iter' num2str(niter)];
  localIters(j) = niter;
end

% initialize Fvals
Fvals = zeros(1,nParEvals);
nZigZagsAll = zeros(nParEvals,1);
term1time = zeros(nParEvals,1);
term2meanRBC = zeros(nParEvals,1);
term3rigid = zeros(nParEvals,1);
term4nZZ = zeros(nParEvals,1);
      
% Now we will submit nParEvals batches
parprofile = o.parCluster;

% EXTRACT THE OPTIMIZATION PARAMETERS (SAVE THEM) and SUBMIT BATCHES
%------------------------------------------------------------------------
oc = curve;

% flags storing whether post shape is smooth and self-intersecting.
flgSmooth = ones(nParEvals,1); % 1: means, it is good to go
flgIntersect = zeros(nParEvals,1); % 0: means, it is good to go
flgSpacing = ones(nParEvals,1); % 1: means, it is good to go
spacings = zeros(2,nParEvals); % keeps Dxs and Dys
Dpostxs = zeros(nParEvals,1); % keeps Dpostxs
Dpostys = zeros(nParEvals,1); % keeps Dpostys
Xposts = zeros(2*o.Nint,nParEvals); % keeps shapes

job = [];

for j = 1 : nParEvals
  % unzip the optimization parameters for each set of parameters
  % we need to get the spline control points to construct a post shape
  % and then calculate Dx and Dy

  % build the post shape using spline control points 
  points = reshape(xorig(:,j),8,2);
  [Xpost,Dpostx,Dposty] = oc.samplePost(o.Nint,points);
  Dpostxs(j,1) = Dpostx; 
  Dpostys(j,1) = Dposty;
  Xposts(:,j) = Xpost;
  % calculate Dx and Dy: spacings between posts in x and y directions
  Dx = o.DLDscale*o.wbox-Dpostx; spacings(1,j) = Dx;
  Dy = o.DLDscale*o.wbox-Dposty; spacings(2,j) = Dy;
  % Also adjust Ufar based on given Dy so that we still have the same flux,
  % and capillary number
  Ufar = o.initFlux/Dy;
  
  % check if Dx and Dy are greater than the lower limit for spacing
  if Dx < o.spaLowLim*o.DLDscale || Dy < o.spaLowLim*o.DLDscale
    flgSpacing(j) = false;
  end
  
  % Check if there is any self-intersection in the post shape
  [xIntersect,~,~]=oc.selfintersect(Xpost);
  if ~isempty(xIntersect)
    flgIntersect(j) = true;
  end
  
  % build a shape object
  shape = capsules(Xpost,[],[],1,1,1);
  % compute the upsampling rate to avoid aliasing
  shape.setUpRate([]);
  % If upsampling rate is higher than 4, the shape is not physical
  if shape.uprate > 6
    flgSmooth(j) = false;
  end
  
  % We want to run DLD simulations only with the posts which have no 
  % intersection (i.e. flgIntersect = false), and are smooth 
  % (i.e. flgSmooth = true), and Dx and Dy are greater than the lower 
  % limiw (i.e. flgSpacing = true). Otherwise, we
  % penalize the set of parameters.
 
  % Submit batches now
  disp('Submitting a batch')
  %submitDLDrunManyVes(flgSmooth(j),flgIntersect(j),...
  %    flgSpacing(j),runNames{j},o.theta,Xpost,Dx,Dy,Dpostx,Dposty,o.VCs,...
  %    o.kappas,o.RAparticle,o.lenParticle,o.volFrac,Ufar,o.Nint,o.Nves,...
  %    o.Next,o.useFMM,o.saveFolder)
  job{j} = batch(parprofile,@submitDLDrunManyVes_Explicit,0,{flgSmooth(j),flgIntersect(j),...
     flgSpacing(j),runNames{j},o.theta,Xpost,Dx,Dy,Dpostx,Dposty,o.VCs,...
      o.kappas,o.RAparticle,o.lenParticle,o.volFrac,Ufar,o.Nint,o.Nves,...
      o.Next,o.useFMM,o.saveFolder});

  %job
end % nParEvals
%------------------------------------------------------------------------

% NOW HERE WAIT FOR THE LAST RUN TO BE COMPLETED and DELETE THE JOBS
%------------------------------------------------------------------------
disp('Waiting for batches to end...')
for j = 1 : nParEvals
  %findJob(parprofile)
  wait(job{j},'finished');
end
delete(findJob(parprofile))
%------------------------------------------------------------------------

% EXTRACT DATA AND EVALUATE OBJECTIVE FUNCTION
%------------------------------------------------------------------------
disp('Evaluating objective function...')
for j = 1 : nParEvals
    
  if ~flgIntersect(j) && flgSmooth(j) && flgSpacing(j)
    % check if file exists. If not, the simulation did not start somehow.
    % So assign NaNs and penalize

    if exist([o.saveFolder runNames{j} '.mat'],'file') == 2
      load([o.saveFolder runNames{j} '.mat']);
      
      % loaded is it: time step, XhistStore{it}(2*N,nv): history of vesicle
      % shapes, initCys and finalCys: y-coordinates of initial and final
      % vesicle configurations, rem_ids: the ids of remaining vesicles,
      % frozen_ids: ids of the frozen vesicles, nPutBacks: number of
      % putBacks of the circular particle, nputLater: number of vesicles to
      % be initialized later but somehow could not initialized, nZigZags:
      % number of vesicles that zig-zagged and were frozen, putBackLimit:
      % the maximum number of putBacks for the circular particle. 
      % zigzag_ids: ids of zig-zagged vesicles
    
      % Get some information about vesicles' dynamics. We want them to
      % zig-zag.
      % Zigzagged vesicles are the ones that have finalCys-initCys <=
      % Dposty 
      
      deltaCys_RBCs = ...
        finalCys(rem_ids(2:numel(rem_ids)))-initCys(rem_ids(2:numel(rem_ids)));
      deltaCy_rigid = mean(XhistStore{it}(end/2+1:end,1)) - initCys(1) + ...
          +(Dpostx+Dx)*tan(o.theta);
      DelYscale = (Dposty+Dy);
      Lscale = ((putBackLimit-1)*4+3)*(Dpostx+Dx);
      timeScale = Lscale/(o.initFlux/Dy);
      processTime = time(it);
      nvesInit = numel(XhistStore{1}(1,:))-1;
      Fvals(j) = 0.5*processTime/timeScale + 0.5*mean(deltaCys_RBCs)/DelYscale - ...
          10*deltaCy_rigid/DelYscale - nZigZags/nvesInit;
      
      term1time(j) = processTime/timeScale;
      term2meanRBC(j) = mean(deltaCys_RBCs)/DelYscale;
      term3rigid(j) = -deltaCy_rigid/DelYscale;
      term4nZZ(j) = -nZigZags/nvesInit;
      
      
%       deltaCys = finalCys-initCys;
%       nvesTotal = numel(rem_ids)+nZigZags;
%       DelYscaling = 4*(Dpostx+Dx)*tan(o.theta)/2+2*(Dposty+Dy);
%       Yend = mean(XhistStore{it}(end/2+1:end,1));
%       fval = -nZigZags/nvesTotal/time(it)+...
%           mean(deltaCys(rem_ids(2:end)))/DelYscaling/time(it)-...
%           (Yend-initCys(1)+(Dpostx+Dx)*tan(o.theta))/Dy;
    
      % goal: make circular one to displace, increase # of zig-zagging
      % vesicles, increase the vertical displacement of vesicles in the -y
      % direction, increase the vertical displacement of circular one in
      % the +y direction.
      
      nZigZagsAll(j) = nZigZags;
      
      % EVALUATE THE OBJECTIVE FUNCTION
      if nZigZags == 0 % simulation crashes
        Fvals(j) = 1E+4;
      elseif ~iPartZZ
        if abs(mean(XhistStore{it}(1:end/2,1))-mean(XwallsInt(1:end/2,1)))>Dx ...
                && nPutBacks<putBackLimit % simulation crashes
           Fvals(j) = 1E+4;
        end
      end
      
      % If the simulation blown up? or any vesicle shape gone crazy?
      % then nPutBack is smaller than putBackLimit, so not a desired shape
    else
      % file does not exist, so assign NaN, now and below penalize that  
      Fvals(j) = 1E+4;
      message = [num2str(localIters(j)) ...
          'th iteration file does not exist, so penalize the shape'];
      o.saveLog(o.logFile,message,1)      
    end % if exist( )
    
    % save the log info
    message = [num2str(localIters(j)) 'th Iteration gives the loss ' ...
      num2str(Fvals(j),'%2.2e') '.'];
    o.saveLog(o.logFile,message,1);
    %----------------------------------------------------------------------
  
  else 
    % penalize it if the shape is self-intersecting or not smooth enough,
    % or spacing is smaller than the lower limit.
    % Here, we do not resample because the shape is way off. Since we do not 
    % want to resample, we need to penalize it with a very large number.
    % The velocity is known to be in
    % the order of 0.1, so we penalize with 1E+4.
    Fvals(j) = 1E+4;
    
    if flgIntersect(j) || ~flgSmooth(j)
    message = ['The shape is not smooth and/or self-intersecting, so '...
          num2str(localIters(j)) 'th iteration has the loss ' ... 
           num2str(Fvals(j),'%2.2e') '.'];
    end
    
    if ~flgSpacing(j) 
    message = ['The spacing is smaller than the lower limit, so '...
          num2str(localIters(j)) 'th iteration has the loss ' ... 
           num2str(Fvals(j),'%2.2e') '.'];
    end

    
    o.saveLog(o.logFile,message,1); 
  end %if ~flgIntersect(j) && flgSmooth(j) && flgSpacing(j)
      
  % save the data: post shape, Dx, Dy and value of objective function
  save([o.fileData '_nIter' num2str(localIters(end))],'Xposts','spacings',...
      'Fvals','xorig','localIters','nZigZagsAll','term1time','term2meanRBC',...
      'term3rigid','term4nZZ')
end % j = 1 : nParEvals
%------------------------------------------------------------------------
tIter = toc(tIter);
% save the log info
if nParEvals == 1
  message = ['Only one simulation is performed ' ...
      '(either for mean or resampled params.)']; ...
  o.saveLog(o.logFile,message,1);
  message = ['It took ' num2str(tIter,'%2.2e') ' seconds.'];     
  o.saveLog(o.logFile,message,1);
else
  message = [num2str(nParEvals) ' simulations are performed in parallel.'];
  o.saveLog(o.logFile,message,1);
  message = ['It took ' num2str(tIter,'%2.2e') ' seconds'];
  o.saveLog(o.logFile,message,1);
end

end %end DLD_objFunc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveLog(o,fileName,message,iCont)
    
if iCont % continue writing
  fid = fopen(fileName,'a');
  fprintf(fid,'%s\n',message);
  fclose(fid);
else % open the file for the first time
  fid = fopen(fileName,'w');
  fprintf(fid,'%s\n',message);
  fclose(fid);
end
    
end % end saveLog
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end class
