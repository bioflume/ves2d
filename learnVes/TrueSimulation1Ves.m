function TrueSimulation1Ves(speed)
addpath ../src/
disp('Single vesicle with background fluid, Exact Solve')

% FLAGS
%-------------------------------------------------------------------------
bgFlow = 'parabolic'; % 'shear','tayGreen','relax','parabolic','rotation'
kappa = 1;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
if speed == 100; Th = 20; dt = 1e-3; end; % Ca = 2.5
if speed == 200; Th = 10; dt = 1e-3; end; % Ca = 5
if speed == 400; Th = 5; dt = 1e-3; end; % Ca = 10

if speed == 8000; Th = 0.25; dt = 1e-4; end;
if speed == 12000; Th = 0.15; dt = 1e-5; end;
if speed == 16000; Th = 0.15; dt = 1e-5; end;
if speed == 24000; Th = 0.05; dt = 1e-5; end;
if speed == 30000; Th = 0.05; dt = 1e-5; end;

N =  128; % num. points
% dt = 1e-4/(speed/8000); % time step size
oc = curve;
op = poten(N);
dnn = dnnTools;
%-------------------------------------------------------------------------

% VESICLE:
% -------------------------------------------------------------------------
X0 = oc.initConfig(N,'ellipse');
[~,~,len] = oc.geomProp(X0);
X0 = X0./len;
IA = pi/2;
cent = [0; 0.065];
X = zeros(size(X0));
X(1:N) = cos(IA) * X0(1:N) - ...
      sin(IA) * X0(N+1:2*N) + cent(1);
X(N+1:2*N) = sin(IA) * X0(1:N) +  ...
      cos(IA) * X0(N+1:2*N) + cent(2);

[~,area0,len0] = oc.geomProp(X);
% -------------------------------------------------------------------------

% BACKGROUND VELOCITY
%-------------------------------------------------------------------------
vinf = dnn.setBgFlow(bgFlow,speed);

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
% folderName = '/work2/03353/gokberk/frontera/truePoisRuns/';
folderName = './output/';
fileName = [folderName 'speed' num2str(speed) '.bin'];

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

  vesicle = capsules(XhistTrue,[],[],kappa,1,1);
  vesicle.setUpRate();
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
  
  % figure(1);clf;
  % hold on;
  % x = [XhistTrue(1:end/2); XhistTrue(1)];
  % y = [XhistTrue(1+end/2:end); XhistTrue(end/2+1)];
  % plot(x,y,'r','linewidth',2)
  % plot(XhistTrue(1), XhistTrue(end/2+1),'o','markerfacecolor','r','markersize',8)
  % axis equal
  % pause(0.1)

  if rem(it,10) == 0
    writeData(fileName,XhistTrue,sigStore,timeTrue(end),0,0);  
    % figure(2); clf;
    % plot(cx, cy, 'linewidth',2)
    % axis square
    % grid
    % xlabel('center in x')
    % ylabel('center in y')
    % title(['Time: ' num2str(timeTrue(it))])
  end
  
end % while
end % end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(filename,X,sigma,time,ncountNN,ncountExact)
x = X(1:end/2,:);
y = X(end/2+1:end,:);
output = [time;ncountNN;ncountExact;x(:);y(:);sigma(:)];

fid = fopen(filename,'a');
fwrite(fid,output,'double');
fclose(fid);


end




