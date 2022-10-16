function Xhist = MLARM(X,Xwalls,prams,options)

disp(['Solving with DNNs?:' num2str(~options.exactSolve)'])
disp(['Background flow: ' prams.bgFlow])
disp(['Flow scaling: ' num2str(prams.speed)])
disp(['Time horizon: ' num2str(prams.Th)])
disp(['Time step size: ' num2str(prams.dt)])
disp(['Bending stiffness: ' num2str(prams.kappa)])
disp(['Num. points per vesicle: ' num2str(prams.N)])
disp(['Num. of vesicles: ' num2str(prams.nv)])
if ~isempty(Xwalls)
disp(['Num. points per wall: ' num2str(prams.Nbd)])
disp(['Num. of walls: ' num2str(prams.nvbd)])
end

oc = curve;

% BUILD DNN CLASS
% -------------------------------------------------------------------------
dnnOpts.confined = ~isempty(Xwalls); % confined flow?
dnnOpts.Nnet = 256; % num. points per vesicle at which nets trained
dnnOpts.nVelModes = 24; % # of modes for M*vinf's network
dnnOpts.nTenModes = 32; % number of modes for inv(DivGT)*(DivGB)x
dnnOpts.nComp = 16; % num. of components for networks except relaxation problem
dnnOpts.nCompRelax = options.nCompRelax;
% Set the remaining MLARM parameters

prams = initMLARM(prams);

dnn = dnnClass(X,Xwalls,[],[],prams,dnnOpts);
tt = dnn.tt; dnn.oc = oc; % time step and curve classes
% load PCA basis vectors and column mean (evects, colMeans) for Nnet = 256
load necessaryMatFiles/pcaBasisNewest.mat
dnn.colMeans = colMeans; dnn.evects = evects; 


% INITIALLY TAKE SMALL TIME STEPS WITH IMPLICIT SOLVER
% ------------------------------------------------------------------------
% necessary to find correct initial tension, density, rotlet and stokeslet
tt.dt = min(prams.dt,1E-5); sig = zeros(prams.N,prams.nv); 
eta = zeros(size(Xwalls)); RS = zeros(3,prams.nvbd);

for iter = 1 : 2
  vesicle = capsules(X,[],[],prams.kappa,ones(prams.nv,1),1); vesicle.setUpRate();
  [X,sig,eta,RS] = tt.timeStepSimple(X,sig,eta,RS,ones(prams.nv,1),dnn.walls,vesicle);
end
tt.dt = prams.dt;
% ------------------------------------------------------------------------

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
time = (0:prams.dt:prams.Th)'; ntime = numel(time);
Xhist = zeros(2*prams.N,prams.nv,ntime); Xhist(:,:,1) = X; 
sigStore = zeros(prams.N,prams.nv,ntime); sigStore(:,:,1) = sig;
errAreaLength = zeros(ntime,1);
etaStore = zeros(2*prams.Nbd,prams.nvbd,ntime); 
RSstore = zeros(3,prams.nvbd,ntime);
etaStore(:,:,1) = eta; RSstore(:,:,1) = RS;

% Save in a file:
save(options.fileName,'Xhist','Xwalls','time','sigStore',...
        'etaStore','errAreaLength')

% TIME STEPPING
for it = 2 : ntime
  disp('********************************************') 
  disp([num2str(it) 'th (/' num2str(ntime) ...
    ') time step, time: ' num2str(time(it))])

  
  if options.exactSolve     
    % EXACT SOLVE  
    disp('Taking an exact time step...'); tStart = tic;  
    vesicle = capsules(Xhist(:,:,it-1),[],[],prams.kappa,ones(prams.nv,1),1);
    vesicle.setUpRate();
    [Xhist(:,:,it),sigStore(:,:,it),etaStore(:,:,it),RSstore(:,:,it),...
        ~,~] = tt.timeStepSimple(Xhist(:,:,it-1),...
        sigStore(:,:,it-1),etaStore(:,:,it-1),RSstore(:,:,it-1),...
        ones(prams.nv,1),dnn.walls,vesicle);
  else 
    % MLARM SOLVE 
    disp('Taking a step with DNNs...');  tStart = tic;
    [Xhist(:,:,it),sigStore(:,:,it),etaStore(:,:,it),~,~,...
        RSstore(:,:,it)] = dnn.DNNsolveAltern(Xhist(:,:,it-1),...
        sigStore(:,:,it-1),etaStore(:,:,it-1),[],[],RSstore(:,:,it-1));
  end % if options.exactSolve
  
  
  % JIGGLING
  if options.jiggle % jiggle vesicles pointwise if there is a near-collision
    disp('Handling collisions...')  
    Xhist(:,:,it) = oc.fixCollisionsZinc(Xhist(:,:,it),Xwalls);
    % Equally distribute points in arc-length b/c jiggling messes up points
    disp('Equally distributing points on arc-length...') 
    Xiter = Xhist(:,:,it);
    for iter = 1 : 5
      [Xiter,~,~] = oc.redistributeArcLength(Xiter);
    end
    % Fix misalignment in center and angle due to reparametrization
    Xhist(:,:,it) = oc.alignCenterAngle(Xhist(:,:,it),Xiter);
  end
  
  % AREA-LENGTH CORRECTION
  disp('Correcting area-length errors...')
  [Xnew,ifail] = oc.correctAreaAndLength2(Xhist(:,:,it),prams.area0,prams.len0);
  if ifail
    disp('AREA-LENGTH CANNOT BE CORRECTED!!!')
  end
  Xhist(:,:,it) = oc.alignCenterAngle(Xhist(:,:,it),Xnew); 

  % check if shape is intersecting
  for k = 1 : prams.nv
    [xIntersect,~,~] = oc.selfintersect(Xhist(:,k,it));
    if ~isempty(xIntersect)
      disp('NEW VESICLE SHAPE IS SELF-INTERSECTING!!!')
    end
  end
  
  % Equally distribute points in arc-length
  disp('Equally distributing points on arc-length...')
  Xiter = Xhist(:,:,it);
  for iter = 1 : 5
    [Xiter,~,~] = oc.redistributeArcLength(Xiter);
  end
  
  % Fix misalignment in center and angle due to reparametrization
  Xhist(:,:,it) = oc.alignCenterAngle(Xhist(:,:,it),Xiter);
  
  % Compute error in area and length
  [~,area,len] = oc.geomProp(Xhist(:,:,it));
  errArea = max(abs(area-prams.area0)./prams.area0); 
  errLen = max(abs(len-prams.len0)./prams.len0);
  errAreaLength(it) = max(errArea,errLen);
  
  if rem(it,100) == 0
    disp('Saving data...')  
    save(options.fileName,'Xhist','it','sigStore','RSstore',...
        'etaStore','errAreaLength','-append')
  end
  disp(['Time step took ' num2str(toc(tStart)) ' seconds.'])
  disp('********************************************') 
  disp(' ')
  
  if options.plotOnFly
  figure(1);clf;
  vecx = [Xwalls(1:end/2,:);Xwalls(1,:)];
  vecy = [Xwalls(end/2+1:end,:);Xwalls(end/2+1,:)];
  plot(vecx,vecy,'k')
  hold on;
  vecx = [interpft(Xhist(1:end/2,:,it),128);Xhist(1,:,it)];
  vecy = [interpft(Xhist(end/2+1:end,:,it),128);Xhist(end/2+1,:,it)];
  plot(vecx,vecy,'r','linewidth',2)
  axis equal
  pause(0.1)
  end
  
end % for it = 2 : ntime
save(options.fileName,'Xhist','it','sigStore','RSstore',...
        'etaStore','errAreaLength','-append') 
end % end MLARM
