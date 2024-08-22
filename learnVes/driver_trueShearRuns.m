function driver_trueShearRuns(speed, timeHorizon, dtGiven)
addpath ../src/
disp('Single vesicle with background fluid, Exact Solve')

% FLAGS
%-------------------------------------------------------------------------
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------

N =  128; % num. points
dt = dtGiven;
Th = timeHorizon;
oc = curve;


%-------------------------------------------------------------------------

% VESICLE:
% -------------------------------------------------------------------------
cx = [-0.3; 0];
cy = [0.040; 0];

IA = [0; pi/2];


X0 = oc.initConfig(N,'ellipse');


[~,~,len] = oc.geomProp(X0);
X0 = X0./len;
X = zeros(2*N,2);
for k = 1 : 2
X(1:N,k) = cos(IA(k)) * X0(1:N) - ...
      sin(IA(k)) * X0(N+1:2*N) + cx(k);
X(N+1:2*N,k) = sin(IA(k)) * X0(1:N)  + ...
      cos(IA(k)) * X0(N+1:2*N) + cy(k);
end
[~,area0,len0] = oc.geomProp(X);
nv = 2;

prams.bgFlow = 'shear'; % 'shear','tayGreen','relax','parabolic'
prams.speed = speed; % 500-3000 for shear, 70 for rotation, 100-400 for parabolic
prams.Th = Th;
prams.N = N; % num. points for true solve in DNN scheme
prams.Nfmm = N;
prams.nv = 2; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.dt = dt; % time step size
prams.dtRelax = prams.dt;
prams.Nbd = 0;
prams.nvbd = 0;
prams.interpOrder = 1;
prams.chanWidth = 0;
dnn = dnnToolsManyVesFree(X,prams);
% -------------------------------------------------------------------------


% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
% folderName = '/work2/03353/gokberk/frontera/truePoisRuns/';
fileName = ['./output/N128again_shearTrueRuns_dt' num2str(dt) '_speed' num2str(speed) '.bin'];

fid = fopen(fileName,'w');
output = [N;nv];
fwrite(fid,output,'double');

x = X(1:end/2,:); y = X(end/2+1:end,:);
output = [x(:); y(:)];
fwrite(fid,output,'double');
fclose(fid);

timeTrue = [0];
XhistTrue = X; sigStore = zeros(N,nv);

writeData(fileName,XhistTrue,sigStore,timeTrue(end),0,0);

% TIME STEPPING
it = 1;
while timeTrue(end) < Th
  disp('********************************************') 
  disp([num2str(it) 'th time step, time: ' num2str(timeTrue(it))])
  
  tStart = tic;
  disp('---------------------------')
  
  % SOLVE WITHOUT OPERATOR SPLITTING  

  [Xnew,sigStore] = dnn.DNNlikeExactSolve(XhistTrue,sigStore);
  
  
  % AREA-LENGTH CORRECTION
  [Xnew2,ifail] = oc.correctAreaAndLength2(Xnew,area0,len0);
  if ifail
    disp('Error in AL cannot be corrected!!!')
  else
    Xnew = oc.alignCenterAngle(Xnew, Xnew2);
  end
  
  [Xiter,~] = oc.reparametrize(Xnew,[],6,20);
  XhistTrue = oc.alignCenterAngle(Xnew,Xiter);

  it = it + 1;
  timeTrue(it) = timeTrue(it-1) + dt;  


  disp(['took ' num2str(toc(tStart)) ' seconds.'])
  disp('---------------------------')    
  
  % Compute error in area and length
  [~,area,len] = oc.geomProp(XhistTrue);
  errALTrue = max(abs(area-area0)/area0,abs(len-len0)/len0);
  disp(['Error in area and length: ' num2str(errALTrue)])   
  disp('********************************************') 
  disp(' ')
  
  

  if rem(it,10) == 0
    writeData(fileName,XhistTrue,sigStore,timeTrue(end),0,0);  
    % figure(1);clf;
    % hold on;
    % x = [XhistTrue(1:end/2,:); XhistTrue(1,:)];
    % y = [XhistTrue(1+end/2:end,:); XhistTrue(end/2+1,:)];
    % plot(x,y,'r','linewidth',2)
    % hold on
    % plot(XhistTrue(1,:), XhistTrue(end/2+1,:),'o','markerfacecolor','r','markersize',8)
    % xlim([-1 1])
    % ylim([-1 1])
    % axis equal
    % pause(0.1)
  end
  
end % while
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(filename,X,sigma,time,ncountNN,ncountExact)
x = X(1:end/2,:);
y = X(end/2+1:end,:);
output = [time;ncountNN;ncountExact;x(:);y(:);sigma(:)];

fid = fopen(filename,'a');
fwrite(fid,output,'double');
fclose(fid);


end




