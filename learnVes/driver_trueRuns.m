function driver_trueRuns(chanWidth, speed, timeHorizon,dtGiven)
addpath ../src/
disp('Single vesicle with background fluid, Exact Solve')

% FLAGS
%-------------------------------------------------------------------------
bgFlow = 'parabolic'; % 'shear','tayGreen','relax','parabolic','rotation'
kappa = 1;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------

N =  32; % num. points
dt = dtGiven;
Th = timeHorizon;
oc = curve;
op = poten(N);
dnn = dnnTools;

%-------------------------------------------------------------------------

% VESICLE:
% -------------------------------------------------------------------------
if 0
load('./necessaryMatFiles/X100KinitShapes.mat')
X0 = Xstore(:,88);
X0 = [interpft(X0(1:end/2),N);interpft(X0(end/2+1:end),N)];
else
X0 = oc.initConfig(N,'ellipse');
end
% 
[~,~,len] = oc.geomProp(X0);
X0 = X0./len;
IA = 0;
cent = [0; chanWidth/2];
X = zeros(size(X0));
X(1:N) = cos(IA) * X0(1:N) - ...
      sin(IA) * X0(N+1:2*N) + cent(1);
X(N+1:2*N) = sin(IA) * X0(1:N) +  ...
      cos(IA) * X0(N+1:2*N) + cent(2);
[~,area0,len0] = oc.geomProp(X);


% 
% load equilX_AllFour_Speed400_Widthp17213
% [~,area0,len0] = oc.geomProp(X);
% -------------------------------------------------------------------------

% BACKGROUND VELOCITY
%-------------------------------------------------------------------------
vinf = @(X) [speed*(1-(X(end/2+1:end,:)/chanWidth).^2);...
      zeros(size(X(1:end/2,:)))];
% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
% folderName = '/work2/03353/gokberk/frontera/truePoisRuns/';
fileName = ['./output/32modes_poisTrueRuns_dt' num2str(dt) '_speed' num2str(speed) '_width' num2str(chanWidth) '.bin'];

fid = fopen(fileName,'w');
output = [N;1];
fwrite(fid,output,'double');

x = X(1:end/2,:); y = X(end/2+1:end,:);
output = [x(:); y(:)];
fwrite(fid,output,'double');
fclose(fid);

timeTrue = [0];
XhistTrue = X; sigStore = zeros(N,1);

writeData(fileName,XhistTrue,sigStore,timeTrue(end),0,0);

% TIME STEPPING
it = 1;
cx = []; cy = [];
while timeTrue(end) < Th
  disp('********************************************') 
  disp([num2str(it) 'th time step, time: ' num2str(timeTrue(it))])
  
  tStart = tic;
  disp('---------------------------')
  
  % SOLVE WITHOUT OPERATOR SPLITTING  
  disp('Solving without operator splitting...') 

  vesicle = capsules(XhistTrue,[],[],kappa,1,0);

  Xnew = dnn.relaxExactSolve(vesicle,vinf(XhistTrue),dt,XhistTrue,op);
  
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

  cx = [cx; mean(XhistTrue(1:end/2))];
  cy = [cy; mean(XhistTrue(end/2+1:end))];


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




