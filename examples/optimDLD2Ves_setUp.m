classdef optimDLD2Ves_setUp
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
fileParams % file to save spline control points
saveFolder % folder in which we save

period     % period
initFlux   % initial Dy*Ufar, scale Ufar as Dy changes

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
function o = optimDLD2Ves_setUp(prams,options)
    
o.runName = prams.runName;
o.period = prams.period;

o.Next = prams.Next;
o.Nint = prams.Nint;
o.Nves = prams.Nves;
o.Ufar = prams.Ufar;
o.VCs = prams.VCs;
o.kappas = prams.kappas;

o.initFlux = prams.Ufar*prams.Dy;

o.useFMM = options.useFMM;
o.useHODLR = options.useHODLR;
o.HODLRforW2W = options.HODLRforW2W;
o.repulsion = options.repulsion;

o.fileData = prams.fileData;
o.logFile = prams.logFile;
o.saveFolder = prams.saveFolder;
o.fileParams = prams.fileParams;

o.wbox = prams.wbox;
o.spaLowLim = prams.spaLowLim;

% initialize the files
o.saveData(o.fileData,[o.Nint;o.wbox;o.initFlux],0);

message = ['Optimization for post shape and wbox is ' num2str(o.wbox) '.'];
paramSize = 16;
o.saveLog(o.logFile,message,0);
o.saveData(o.fileParams,paramSize,0);

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
  Dx = o.wbox-Dpostx; spacings(1,j) = Dx;
  Dy = o.wbox-Dposty; spacings(2,j) = Dy;
  
  % Also adjust Ufar based on given Dy so that we still have the same flux,
  % and capillary number
  Ufar = o.initFlux/Dy;
  
  % check if Dx and Dy are greater than the lower limit for spacing
  if Dx < o.spaLowLim || Dy < o.spaLowLim
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
  job{j} = batch(parprofile,@submitDLDrun2Ves,0,{flgSmooth(j),flgIntersect(j),...
      flgSpacing(j),runNames{j},o.period,Xpost,Dx,Dy,o.VCs,o.kappas,Ufar,...
      o.Nint,o.Nves,o.Next,o.useFMM,o.useHODLR,o.HODLRforW2W,...
      o.repulsion,o.saveFolder});

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
    % check if file exists. If not, it means that the post shape is fine
    % to run a DLD simulation but somehow we could not perform the
    % simulation, and we do not have a file. So
    % just ignore this simulation and keep going.

    if exist([o.saveFolder runNames{j} '_Data.bin'],'file') == 2
      [posx,posy,~,~,~,~,intWallx,intWally,...
          time,~,~] = loadFile([o.saveFolder runNames{j} '_Data.bin']);
    
      % get the vesicles' trajectories (we have 2 vesicles in a sim.)
      centx = zeros(numel(time),2); centy = zeros(numel(time),2);
      for tt = 1 : numel(time)
        for k = 1 : 2
          centx(tt,k) = mean(interpft(posx(:,k,tt),96));
          centy(tt,k) = mean(interpft(posy(:,k,tt),96));
        end % for k
      end % for tt
      
      
      if numel(time)>1
        % get number of putbacks and lateral displacement between initial
        % and the last states
        nPutBacks = zeros(2,1);
%         delYs = zeros(2,1);
        finalYs = zeros(2,1);
        for k = 1 : 2
          nPutBacks(k) = numel(find(diff(centx(:,k))<=-Dpostxs(j)));
%           delYs(k) = centy(end,k)-centy(1,k);
          % final lat. pos. w.r.t. to the top of the pillar below
          finalYs(k) = centy(end,k)-(centy(1,k)-0.5*spacings(2,j));
        end
      else
        % Vesicles collide with pillars when they are initialized
        nPutBacks = ones(2,1);
%         delYs = [NaN; NaN];
        finalYs = [NaN; NaN];
      end
      
      % If centy(end,1)<0, centy(end,2)<-Dy-Dposty(== wbox), they both zig-zagged
      % If the above is not correct, and nPutBack ==
      % ceil(prams.periodN/2)+1, then the cell displaces.
      % If none of the above is the case, then the simulation has failed.
      % So check these conditions to find out if the simulation ran through
      
      iZZ = zeros(2,1); % flag showing that cell zig-zagged
      
      % FIRST CELL WHICH IS SUPPOSED TO DISPLACE
      if centy(end,1) < 1.1*(Dpostys(j)/2+spacings(2,j)/2) 
        iZZ(1) = 1;
      else
        if nPutBacks(1) >= o.period
          iZZ(1) = 0;
        else
          iZZ(1) = NaN;
        end
      end
      
      % SECOND CELL WHICH IS SUPPOSED TO ZIG-ZAG
      if centy(end,2) < 0.9*(-Dpostys(j)/2-spacings(2,j)/2) 
        iZZ(2) = 1;
      else
        if nPutBacks(2) >= o.period
          iZZ(2) = 0;
        else
          iZZ(2) = NaN;
        end
      end
      
      % IF WE CAN DETERMINE IF THE CEL ZIG-ZAGGED OR NOT, THEN WE CAN
      % COMPUTE THE OBJECTIVE FUNCTION
      
      if ~any(isnan(iZZ)) && ~any(isnan(finalYs))
        % normalize with the initial delY w.r.t. to the post's center  
%         Fvals(j) = -delYs(1)/(0.5*spacings(2,j)) ...
%             + 1/(1+nPutBacks(2))*delYs(2)/...
%             (0.5*spacings(2,j)+0.5*Dpostys(j));

        if iZZ(1) == iZZ(2) % means that no separation
          % Fvals is not greater than 1 for separation cases, so 10 is
          % large enough to use as a penalty constant. The other term
          % evaluates how different the cell dynamics are. The smaller the
          % term is, the more different the dynamics are and the closer to
          % get separated the cells are.
          if iZZ(1) == 1 % both zig-zag
            Fvals(j) = 10 + (-abs(diff(nPutBacks))/(o.period+2)); 
          else % both displacement
            Fvals(j) = 10 + (-abs(diff(finalYs))/spacings(2,j));
          end
        else % means that separation occurs
          % minimize nPutBacks for ZZ cell and maximize final lat. pos. 
          % D cell, we do not care which one is doing what as long as
          % separation occurs
          Fvals(j) = 0.2*nPutBacks(iZZ==1)/(o.period+2)-finalYs(iZZ==0)/spacings(2,j);
        end
           
      else
        Fvals(j) = 1E+4;
        message = [num2str(localIters(j)) ...
            'th iteration could not be completed, ' ...
            'so penalize this shape'];
          o.saveLog(o.logFile,message,1)
      end

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
  o.saveData(o.fileData,[Xposts(:,j);spacings(:,j);Fvals(j)],1);
  % also save the b-spline control points
  o.saveData(o.fileParams,xorig(:,j),1); % so that we can rerun
  
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
function saveData(o,fileName,X,iCont)
    
if iCont % continue writing
  fid = fopen(fileName,'a');
  fwrite(fid,X,'double');
  fclose(fid);
else % open the file for the first time
  fid = fopen(fileName,'w');
  fwrite(fid,X,'double');
  fclose(fid);
end
    
end % end saveData
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
