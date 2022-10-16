function MLARM_DLD(X,XwallsInt,XwallsExt,prams,options)

disp(['Solving with DNNs?:' num2str(~options.exactSolve)'])
disp(['Background flow: ' prams.bgFlow])
disp(['Flow scaling: ' num2str(prams.speed)])
disp(['Time horizon: ' num2str(prams.Th)])
disp(['Time step size: ' num2str(prams.dt)])
disp(['Bending stiffness: ' num2str(prams.kappa)])
disp(['Num. points per vesicle: ' num2str(prams.N)])
disp(['Num. of vesicles: ' num2str(prams.nv)])
disp(['Num. of pillars: ' num2str(prams.nvbdInt)])
disp(['Num. points per pillar: ' num2str(prams.NbdInt)])
disp(['Num. points on exterior wall: ' num2str(prams.NbdExt)])
oc = curve;

% BUILD DNN CLASS
% -------------------------------------------------------------------------
dnnOpts.confined = true; % confined flow?
dnnOpts.Nnet = 256; % num. points per vesicle at which nets trained
dnnOpts.nVelModes = 24; % # of modes for M*vinf's network
dnnOpts.nTenModes = 32; % number of modes for inv(DivGT)*(DivGB)x
dnnOpts.nComp = 16; % num. of components for networks except relaxation problem
dnnOpts.nCompRelax = options.nCompRelax;
% Set the remaining MLARM parameters

prams = initMLARM(prams);

dnn = dnnClass(X,[],XwallsInt,XwallsExt,prams,dnnOpts);
tt = dnn.tt; dnn.oc = oc; % time step and curve classes
% load PCA basis vectors and column mean (evects, colMeans) for Nnet = 256
load necessaryMatFiles/pcaBasisNewest.mat
dnn.colMeans = colMeans; dnn.evects = evects; 

% If there are circular particles, precompute matrices related to them
circularPartExists = false;
if prams.nCircular > 0 % if so, then they are in the first nCircular columns of X
  circularPartExists = true;  
  [Xstand,~,~,~,~] = dnn.standardizationStep(X(:,1));
  vesicle = capsules(Xstand,[],[],1,1,prams.kappa);
  vesicle.setUpRate([]);
  opCirc = poten(dnnOpts.Nnet);
  GCirc = opCirc.stokesSLmatrix(vesicle);
  [BenCirc,TenCirc,DivCirc] = vesicle.computeDerivs;
  dnn.invTenMatCirc = inv(DivCirc*GCirc*TenCirc)*DivCirc;
  dnn.selfBendMatCirc = inv(DivCirc*GCirc*TenCirc)*(DivCirc*GCirc*(-BenCirc))*prams.kappa;
  dnn.Mcirc = GCirc*TenCirc*((DivCirc*GCirc*TenCirc)\eye(dnnOpts.Nnet))*DivCirc;
  dnn.relaxMatCirc = (eye(2*dnnOpts.Nnet)-prams.kappa*prams.dt*...
      (-GCirc*BenCirc+dnn.Mcirc*GCirc*BenCirc))\eye(2*dnnOpts.Nnet);
end

% INITIALLY TAKE SMALL TIME STEPS WITH IMPLICIT SOLVER
% ------------------------------------------------------------------------
% necessary to find correct initial tension, density, rotlet and stokeslet
tt.dt = min(prams.dt,1E-5); sig = zeros(prams.N,prams.nv); 
etaInt = zeros(2*prams.NbdInt,prams.nvbdInt); 
etaExt = zeros(2*prams.NbdExt,1); RS = zeros(3,prams.nvbdInt+1);

Xfrozen = []; frozen_ids = []; rem_ids = [1:prams.nv]'; nputLater = 0;
if prams.nCircular > 0
  nPutBacks = zeros(prams.nCircular,1);
else
  nPutBacks = zeros(prams.nv,1);
end

initCys = zeros(prams.nv,1);
for k = 1 : prams.nv
  initCys(k) = mean(X(end/2+1:end,k));
end
finalCys = zeros(prams.nv,1);


for iter = 1 : 2
  vesicle = capsules(X,[],[],prams.kappa,ones(prams.nv,1),1); vesicle.setUpRate();
  [X,sig,etaInt,etaExt,RS] = tt.timeStepSimpleDiffDisc(X,sig,etaInt,etaExt,...
      RS,ones(prams.nv,1),dnn.wallsInt,dnn.wallsExt,vesicle);
end
Xstore{1} = X;
tt.dt = prams.dt;
% ------------------------------------------------------------------------

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
time = (0:prams.dt:prams.Th)'; ntime = numel(time);
errAreaLength = zeros(ntime,1);

% Save in a file:
save(options.fileName,'Xstore','XwallsInt','XwallsExt','time',...
    'it','errAreaLength','-append')

% TIME STEPPING
for it = 2 : ntime
  disp('********************************************') 
  disp([num2str(it) 'th (/' num2str(ntime) ...
    ') time step, time: ' num2str(time(it))])

  
  if options.exactSolve     
    % EXACT SOLVE  
    disp('Taking an exact time step...'); tStart = tic;  
    vesicle = capsules(X,[],[],prams.kappa,ones(prams.nv,1),1);
    vesicle.setUpRate();
    [X,sig,etaInt,etaExt,RS] = tt.timeStepSimpleDiffDisc(X,sig,etaInt,etaExt,...
      RS,ones(prams.nv,1),dnn.wallsInt,dnn.wallsExt,vesicle);
  else 
    % MLARM SOLVE 
    disp('Taking a step with DNNs...');  tStart = tic;
    [X,sig,~,etaInt,etaExt,RS] = dnn.DNNsolveMixture(X,sig,[],etaInt,...
        etaExt,RS);
  end % if options.exactSolve
  
  
  % JIGGLING
  if options.jiggle % jiggle vesicles pointwise if there is a near-collision
    disp('Handling collisions...')  
    X = oc.fixCollisionsZinc(X,XwallsInt);
    % Equally distribute points in arc-length b/c jiggling messes up points
    disp('Equally distributing points on arc-length...') 
    Xiter = X;
    for iter = 1 : 5
      [Xiter,~,~] = oc.redistributeArcLength(Xiter);
    end
    % Fix misalignment in center and angle due to reparametrization
    X = oc.alignCenterAngle(X,Xiter);
  end
  
  % AREA-LENGTH CORRECTION
  disp('Correcting area-length errors...')
  [Xcorr,ifail] = oc.correctAreaAndLength2(X,prams.area0,prams.len0);
  if ifail
    disp('AREA-LENGTH CANNOT BE CORRECTED!!!')
  end
  X = oc.alignCenterAngle(X,Xcorr); 

  % check if shape is intersecting
  for k = 1 : prams.nv
    [xIntersect,~,~] = oc.selfintersect(X(:,k));
    if ~isempty(xIntersect)
      disp('NEW VESICLE SHAPE IS SELF-INTERSECTING!!!')
    end
  end
  
  % Equally distribute points in arc-length
  disp('Equally distributing points on arc-length...')
  Xiter = X;
  for iter = 1 : 5
    [Xiter,~,~] = oc.redistributeArcLength(Xiter);
  end
  
  % Fix misalignment in center and angle due to reparametrization
  X = oc.alignCenterAngle(X,Xiter);
  
  % Compute error in area and length
  [~,area,len] = oc.geomProp(X);
  errArea = max(abs(area-prams.area0)./prams.area0); 
  errLen = max(abs(len-prams.len0)./prams.len0);
  errAreaLength(it) = max(errArea,errLen);
  
  % Freeze and Stream
  [X,sig,Xfrozen,rem_ids,frozen_ids,prams,nPutBacks,nputLater,nZigZags,...
      initCys,finalCys] = freezeAndStream(oc,X,sig,Xfrozen,prams,...
    rem_ids,frozen_ids,nPutBacks,nputLater,...
    nZigZags,initCys,finalCys);

  % save X into cell structure (as it might be changing size every step)
  Xstore{end+1} = X;
  
  if rem(it,100) == 0
    disp('Saving data...')  
    save(options.fileName,'Xstore','Xfrozen','frozen_ids','rem_ids',...
        'initCys','finalCys','nPutBacks','nputLater','nZigZags','it',...
        'errAreaLength','-append')
  end
  
  disp(['Time step took ' num2str(toc(tStart)) ' seconds.'])
  disp('********************************************') 
  disp(' ')
  
  % Check if we need to terminate the simulation
  if numel(rem_ids) == 0
    disp('All particles have zig-zagged. COMPLETED!')
    save(options.fileName,'Xstore','Xfrozen','frozen_ids','rem_ids',...
        'initCys','finalCys','nPutBacks','nputLater','nZigZags','it',...
        'errAreaLength','-append')
    break;
  end
  
  if min(nPutBacks(rem_ids)) >= prams.putBackLimit
    disp('Vesicles or circular particles have been put back at least the given limit. COMPLETED!')
    save(options.fileName,'Xstore','Xfrozen','frozen_ids','rem_ids',...
        'initCys','finalCys','nPutBacks','nputLater','nZigZags','it',...
        'errAreaLength','-append')
    break;
  end
  
  if circularPartExists && prams.nCircular == 0
    disp('All circular particles are frozen. COMPLETED!')      
    save(options.fileName,'Xstore','Xfrozen','frozen_ids','rem_ids',...
        'initCys','finalCys','nPutBacks','nputLater','nZigZags','it',...
        'errAreaLength','-append')  
  end
  
  if options.plotOnFly
  figure(1);clf;
  vecx = [XwallsExt(1:end/2);XwallsExt(1)];
  vecy = [XwallsExt(end/2+1:end);XwallsExt(end/2+1)];
  plot(vecx,vecy,'k')
  hold on;
  
  vecx = [XwallsInt(1:end/2,:);XwallsInt(1,:)];
  vecy = [XwallsInt(end/2+1:end,:);XwallsInt(end/2+1,:)];
  plot(vecx,vecy,'k')
  
  vecx = [interpft(X(1:end/2,:),128);X(1,:)];
  vecy = [interpft(X(end/2+1:end,:),128);X(end/2+1,:)];
  plot(vecx,vecy,'r','linewidth',2)
  axis equal
  pause(0.1)
  end
  
end % for it = 2 : ntime
 
end % end MLARM_DLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xhist,sigStore,Xfrozen,rem_ids,frozen_ids,prams,nPutBacks,...
    nputLater,nZigZags,initCys,finalCys] = freezeAndStream(oc,X,sigma,...
    Xfrozen,prams,rem_ids,frozen_ids,nPutBacks,nputLater,...
    nZigZags,initCys,finalCys)

% THIS ASSUMES ONLY ONE CIRCULAR PARTICLE IN THE FLOW

% in this one we randomly put vesicle in the entrance if there is a vesicle 
% that reaches to the end of the device. If cannot put a new vesicle, then,
% freezes. In later time steps, it tries to put new vesicles.

[~,~,lenscale] = oc.geomProp(X); lenscale = max(lenscale);

% First check if the first vesicle (RA = 0.9) needs to be translated back
Xlarge = zeros(size(X));
for k = 1 : prams.nv
  cx = mean(X(1:end/2,k)); cy = mean(X(end/2+1:end,k));
  Xlarge(:,k) = [1.1*(X(1:end/2,k)-cx)+cx;1.1*(X(end/2+1:end,k)-cy)+cy];
end
vesicle = capsules(Xlarge,[],[],[],[],1);
vesicle.setUpRate();
NearV2V = vesicle.getZone([],1);
zoneV2V = NearV2V.zone;

pillarGrad = prams.epsilon*(prams.Dposty+prams.Dy)/(prams.Dpostx+prams.Dx);

% if vesicle is beyond xfreeze, then move it back to 1st gap

% ids of the vesicles to be frozen and to be remained
will_freeze=[]; will_remain=[]; 

for k = 1 : prams.nCircular
  centxn = mean(X(1:end/2,k));  
  if centxn>=prams.xfreeze
    centyn = mean(X(end/2+1:end,k));
    centxNew = centxn-prams.moveBackX;
    centyNew = centyn-pillarGrad*prams.moveBackX;
    K = [(1:k-1) (k+1:prams.nv)];
    centxK = mean(X(1:end/2,K),1); centyK = mean(X(end/2+1:end,K),1);
    vesIDsInZone = find(centyK >= centyNew-0.5*prams.Dy & ...
        centyK <= centyNew+0.5*prams.Dy & centxK >= ...
        centxNew-prams.Dx/2-prams.Dpostx/2 & ...
        centxK <= centxNew+prams.Dx/2+prams.Dpostx/2,1);
    % if there is no vesicle in the inlet
    if isempty(vesIDsInZone) && min(X(1:end/2,k)-centxn+centxNew)>prams.xrange
      X(:,k) = [X(1:end/2,k)-centxn+centxNew;X(end/2+1:end,k)-centyn+centyNew]; 
      nPutBacks(rem_ids(k)) = nPutBacks(rem_ids(k))+1; 
      % check if there is any vesicle in near-zone of k, if so, move it back too
      nearZoneIDs = K(sqrt((centxn-centxK).^2+(centyn-centyK).^2)<=3*lenscale);
      for k2 = nearZoneIDs
        J = find(zoneV2V{k2}(:,k) == 1);
        if numel(J)~=0
          centyn = mean(X(end/2+1:end,k2));
          centxn = mean(X(1:end/2,k2));  
          centxNew = centxn-prams.moveBackX;
          centyNew = centyn-pillarGrad*prams.moveBackX;
          if min(X(1:end/2,k2)-centxn+centxNew)>prams.xrange
            X(:,k2) = [X(1:end/2,k2)-centxn+centxNew;X(end/2+1:end,k2)-centyn+centyNew];
            nPutBacks(rem_ids(k2)) = nPutBacks(rem_ids(k2)) +1;
          else % if we do not move k2, freeze it, also freeze those near to it
            will_freeze = [will_freeze;k2];
            nputLater = nputLater + 1;  
          end % if we move k2 near to RA = 0.9
        end % numel(J)~=0
      end % k2 = K
    else
      will_freeze = [will_freeze;K(vesIDsInZone)];
      nputLater = nputLater + numel(vesIDsInZone);
    end % if vesExistsInlet 
  end % if centxn >= prams.xfreeze
end % for k = 1 : prams.nCircular

% regions where a new vesicle can be introduced
xrange = prams.xrange + [0 prams.Dx+prams.Dpostx];
xrangePlace = prams.xrange + (prams.Dx+prams.Dpostx)/2;
yranges(1,:) = prams.Dy/2+prams.Dposty+[0.05*prams.Dy 0.95*prams.Dy]; % first lane
yranges(2,:) = -1.5*prams.Dy-prams.Dposty+[0.05*prams.Dy 0.95*prams.Dy]; % bottom lane
yranges(3,:) = [-0.45 0.45]*prams.Dy; % middle lane

% find vesicles' centers
centxS = mean(X(1:end/2,:),1)';
centyS = mean(X(end/2+1:end,:),1)';

% if we frozen some before, and need to seed again, randomly initialize vesicles
countPutLater = 0;
doItOneMore = true;
while (countPutLater<nputLater && doItOneMore)
  idOfLast = max([rem_ids;frozen_ids]); % id of last vesicle  
  refVesId = randi([prams.nCircular+1,prams.nv]);
  Xref = X(:,refVesId);
  Xref = [Xref(1:end/2)-mean(Xref(1:end/2));Xref(end/2+1:end)-mean(Xref(end/2+1:end))];
  IA = oc.getIncAngle(Xref);
  xref = Xref(1:end/2)*cos(-IA)-Xref(end/2+1:end)*sin(-IA);
  yref = Xref(1:end/2)*sin(-IA)+Xref(end/2+1:end)*cos(-IA);
  Xref = [xref;yref];
  whichLane = randperm(2,2); % only in the bottom or top lanes
  vesPlaced = false; iter = 1;
  while (~vesPlaced && iter <= 2)
    % check if there is any vesicle in entrance of that lane
    vesExistsInlet = isempty(find(centyS >= yranges(whichLane(iter),1) & ...
          centyS <= yranges(whichLane(iter),2) & centxS >= xrange(1) & centxS <= xrange(2),1));
    if vesExistsInlet % then there is no vesicle in the inlet
      centxNew = xrangePlace;
      centyNew = (yranges(whichLane(iter),2)+yranges(whichLane(iter),1))*0.5;
      Xnew = [Xref(1:end/2)+centxNew; Xref(end/2+1:end)+centyNew];  
      vesPlaced = true; 
      countPutLater = countPutLater + 1;
      % update related information
      X(:,end+1) = Xnew;
      rem_ids = [rem_ids;idOfLast+1];
      prams.area0(end+1) = prams.area0(refVesId); 
      prams.len0(end+1) = prams.len0(refVesId);
      prams.nv = prams.nv + 1; 
      prams.viscCont(end+1) = prams.viscCont(refVesId);
      sigma(:,end+1) = sigma(:,refVesId);
      centxS(end+1) = centxNew; centyS(end+1) = centyNew;
      initCys(end+1) = centyNew; finalCys(end+1) = centyNew; 
      nPutBacks(idOfLast+1) = 0;
    else
      iter = iter + 1;
    end
  end % end while      
  doItOneMore = vesPlaced; % if vesicle is placed, then do it one more time, if we need to do so
end % while(countPutLater<nputLater && doItOneMore)
% update the number of vesicles to be put later
nputLater = nputLater - countPutLater;

% Get the near zone again
for k = 1 : prams.nv
  cx = mean(X(1:end/2,k)); cy = mean(X(end/2+1:end,k));
  Xlarge(:,k) = [1.1*(X(1:end/2,k)-cx)+cx;1.1*(X(end/2+1:end,k)-cy)+cy];
end
vesicle = capsules(Xlarge,[],[],[],[],1);
vesicle.setUpRate();
NearV2V = vesicle.getZone([],1);
zoneV2V = NearV2V.zone;

% if vesicle is beyond xfreeze, then move it back to
for k = prams.nCircular+1 : prams.nv
  centxn = mean(X(1:end/2,k));  
  if centxn>=prams.xfreeze
    K = [(prams.nCircular+1:k-1) (k+1:prams.nv)];  
    centyn = mean(X(end/2+1:end,k));
    centxNew = centxn-prams.moveBackX;
    centyNew = centyn-pillarGrad*prams.moveBackX;
    centxK = mean(X(1:end/2,K)); centyK = mean(X(end/2+1:end,K));
    
    vesIDsInZone = find(centyK >= centyNew-0.5*prams.Dy & ...
        centyK <= centyNew+0.5*prams.Dy & centxK >= ...
        centxNew-prams.Dx/2-prams.Dpostx/2 & ...
        centxK <= centxNew+prams.Dx/2+prams.Dpostx/2,1);
    
    % if there is no vesicle in the inlet  
    if isempty(vesIDsInZone) && min(X(1:end/2,k)-centxn+centxNew)>prams.xrange
      X(:,k) = [X(1:end/2,k)-centxn+centxNew;X(end/2+1:end,k)-centyn+centyNew];
      nPutBacks(rem_ids(k)) = nPutBacks(rem_ids(k)) + 1;
      % check if there is any vesicle in near-zone of k, if so, move it back,
      % too
      nearZoneIDs = K(sqrt((centxn-centxK).^2+(centyn-centyK).^2)<=3*lenscale);
      for k2 = nearZoneIDs
        J = find(zoneV2V{k2}(:,k) == 1);
        if numel(J)~=0
          centyn = mean(X(end/2+1:end,k2));
          centxn = mean(X(1:end/2,k2));  
          centxNew = centxn-prams.moveBackX;
          centyNew = centyn-pillarGrad*prams.moveBackX;
          if min(X(1:end/2,k2)-centxn+centxNew)>prams.xrange
            X(:,k2) = [X(1:end/2,k2)-centxn+centxNew;X(end/2+1:end,k2)-centyn+centyNew];  
            nPutBacks(rem_ids(k2)) = nPutBacks(rem_ids(k2)) + 1;
          else
            nputLater = nputLater + 1;
            will_freeze = [will_freeze; k2];
          end
        end
      end
    else % if there is a vesicle where the new one needs to be placed 
      will_freeze = [will_freeze; k];
      nputLater = nputLater + 1;
    end % if vesExistsInlet
  end % if centxn >= prams.xfreeze 
end % for k = 2 : prams.nv

% use different threshold for others
yFreeze = -1.5*prams.Dy-2*prams.Dposty;
for iv = 1:prams.nv
  if mean(X(end/2+1:end,iv)) <= yFreeze
    % tag the ones we need to freeze
    will_freeze = [will_freeze;iv];
    nZigZags = nZigZags + 1;
    message = ['Vesicle ' num2str(rem_ids(iv),'%d') ' is frozen'];
    disp(message)
  end % if mean    
end % for iv

% save the frozen ones and remained ones
if~isempty(will_freeze)
  will_freeze = unique(will_freeze);
  for k = 1 : prams.nv
    if ~any(will_freeze == k)
      will_remain = [will_remain; k];
    end
  end
  % then there are vesicles to be frozen  
  Xfrozen = [Xfrozen X(:,will_freeze)];
  frozen_ids = [frozen_ids;rem_ids(will_freeze)]; 
else
  will_remain = [1:numel(X(1,:))]'; 
end

for k = 1 : prams.nv
  finalCys(rem_ids(k)) = mean(X(end/2+1:end,k));
end
rem_ids = rem_ids(will_remain);
Xhist = X(:,will_remain);
sigStore = sigma(:,will_remain);
prams.viscCont = prams.viscCont(will_remain);
prams.nv = numel(rem_ids);
prams.area0 = prams.area0(will_remain);
prams.len0 = prams.len0(will_remain);
nPutBacks = nPutBacks(will_remain);
prams.nv = numel(rem_ids);


for k = 1 : prams.nCircular
  if ~any(rem_ids == k)
    % then, kth circular particle is frozen      
    prams.nCircular = prams.nCircular - 1;  
  end
end

end % freezeAndStream