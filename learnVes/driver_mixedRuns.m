function driver_mixedRuns(chanWidth,speed,Th)
addpath ../src/
addpath ../examples/
addpath ./shannets/
addpath ./shannets/ves_fft_models/

pathofDocument = fileparts(which('Net_ves_relax_midfat.py'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pathofDocument = fileparts(which('Net_ves_adv_fft.py'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pathofDocument = fileparts(which('ves_fft_mode2.pth'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pe = pyenv('Version', '/Users/gokberk/opt/anaconda3/envs/mattorch/bin/python');


% FLAGS
%-------------------------------------------------------------------------
interpOrder = 5; % must be an odd number
prams.bgFlow = 'parabolic'; % 'shear','tayGreen','relax','parabolic'
prams.speed = speed; 
prams.chanWidth = chanWidth;
iplot = 0;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
errTol = 1e-2;
maxDt = 1e-5;
prams.Th = Th;%2.5; % time horizon
prams.N = 128; % num. points for true solve in DNN scheme
prams.nv = 1; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.interpOrder = interpOrder;
prams.dt = 1E-5; % time step size
prams.dtRelax = prams.dt;
prams.Nbd = 0;
prams.nvbd = 0;

oc = curve;
Th = prams.Th; N = prams.N; nv = prams.nv; dt = prams.dt; 
bgFlow = prams.bgFlow; speed = prams.speed;

% net parameters
Nnet = 256; % num. points
nCompRelax = 32; % number of PCA components for relaxation problem
% # of modes for M*vinf's network, inv(DivGT)*Div*vinf uses the same;
% load PCA basis vectors and column mean (evects, colMeans) for Nnet = 256
% load necessaryMatFiles/pcaCoeffsBasis1step.mat
load necessaryMatFiles/pcaBasisNewest.mat
%-------------------------------------------------------------------------
disp(['Flow: ' prams.bgFlow ', N = ' num2str(N) ', nv = ' num2str(nv) ...
    ', Th = ' num2str(Th)])
%-------------------------------------------------------------------------

% VESICLES and WALLS:
% -------------------------------------------------------------------------
X0 = oc.initConfig(N,'ellipse');
[~,~,len] = oc.geomProp(X0);
X0 = X0./len;
IA = 0;
cent = [0; prams.chanWidth/2];
X = zeros(size(X0));
X(1:N) = cos(IA) * X0(1:N) - ...
      sin(IA) * X0(N+1:2*N) + cent(1);
X(N+1:2*N) = sin(IA) * X0(1:N) +  ...
      cos(IA) * X0(N+1:2*N) + cent(2);
[~,area0,len0] = oc.geomProp(X);

% -------------------------------------------------------------------------

solveType = 'DNN';
fileName = ['./output/mixedNets_DiffNet_poisRuns_speed' num2str(prams.speed) '_width' num2str(prams.chanWidth) '.bin'];
fid = fopen(fileName,'w');
output = [N;nv];
fwrite(fid,output,'double');
x = X(1:end/2,:); y = X(end/2+1:end,:);
output = [x(:); y(:)];
fwrite(fid,output,'double');
fclose(fid);

% BUILD DNN CLASS
% -------------------------------------------------------------------------
dnn = dnnToolsSingleVes(X,prams);
tt = dnn.tt; dnn.oc = oc; dnn.variableKbDt = false;
load ./shannets/ves_fft_in_param.mat
load ./shannets/ves_fft_out_param.mat
dnn.torchAdvInNorm = in_param;
dnn.torchAdvOutNorm = out_param;

% -------------------------------------------------------------------------

% LOAD NETWORKS (SEE THE END OF THE CODE FOR THESE FUNCTIONS)
% -------------------------------------------------------------------------
% LOAD THE COM NET
load ./netRelaxDrift_dt1e-05_8modes_Network.mat
% this loads the network 
dnn.driftNet = net;
dnn.driftNet_muChan = muChan1;
dnn.driftNet_sdevChan = sdevChan1;
dnn.driftNet_muOutput = muOutput;
dnn.driftNet_sdevOutput = stdOutput;

dnn.nCompRelax = nCompRelax;
% LOAD PCA Network for Relaxation Problem
dnn = loadAllPCAnets4Relax(dnn);
% LOAD FFT Based Network for M*Vinf 
% save PCA matrices 
dnn.colMeans = colMeans; dnn.evects = evects; 
% -------------------------------------------------------------------------

% INITIALLY TAKE SMALL TIME STEPS WITH IMPLICIT SOLVER
% ------------------------------------------------------------------------
% necessary to find correct initial tension, density, rotlet and stokeslet
tt.dt = maxDt; sig = zeros(N,nv); 
% ------------------------------------------------------------------------

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
time = [0];
Xhist = X; sigStore = sig; 
errALPred = 0;
ncountCNN = 0;
ncountExct = 0;
driftyNet = [];
driftyAdv = [];
% ------------------------------------------------------------------------
writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);

% TIME STEPPING
it = 1;
cx = []; cy = [];
exactCnter = 0; iExact = false;
while time(end) < prams.Th
  disp('********************************************') 
  disp([num2str(it) 'th time step, time: ' num2str(time(it))])
  
  shouldSolveExact = 0;
  
  disp('Taking a step with DNNs...');  tStart = tic;    
  [Xnew, dyNet, dyAdv] = dnn.DNNsolveMixedNets(Xhist,area0,len0,iExact);
  
  exactCnter = exactCnter + 1;
  % if rem(exactCnter+1,exactFreq) == 0
  %   iExact = true;
  %   exactCnter = 0;
  % else
  %   iExact = false;
  % end

  driftyNet = [driftyNet;dyNet];
  driftyAdv = [driftyAdv;dyAdv];

  % figure(2);clf;
  % plot(driftyNet,'linewidth',2)
  % hold on
  % plot(driftyAdv,'linewidth',2)
  % legend('Relax','Adv')
  % grid
  % axis square
  % pause(0.1)
  
  [xIntersect,~,~] = oc.selfintersect(Xnew);
  if ~isempty(xIntersect); disp('New vesicle shape is self-intersecting!!!'); end;
  
  [~,area,len] = oc.geomProp(Xnew);
  errArea = max(abs(area-area0)./area0); errLen = max(abs(len-len0)./len0);
  
  it = it + 1;
  Xhist = Xnew;
  errALPred(it) = max(errArea,errLen);
  time(it) = time(it-1) + prams.dt;  
  

  cx = [cx; mean(Xhist(1:end/2))];
  cy = [cy; mean(Xhist(end/2+1:end))];
  
  
  
  disp(['Error in area and length: ' num2str(max(errArea, errLen))])   
  disp('********************************************') 
  disp(' ')
  
  if rem(it,10) == 0
    writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);  
  end

  disp(['took ' num2str(toc(tStart)) ' seconds.'])
  
  if iplot
  figure(1);clf;
  hold on;
  x = [Xhist(1:end/2); Xhist(1)];
  y = [Xhist(1+end/2:end); Xhist(end/2+1)];
  plot(x,y,'r','linewidth',2)
  plot(Xhist(1), Xhist(end/2+1),'o','markerfacecolor','r','markersize',8)
  xlim([-1 1])
  ylim([-1 1])
  axis equal
  pause(0.1)
  end
end

% Save data to a mat-file:
writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dnn = loadAllPCAnets4Relax(dnn)
interpOrder = dnn.interpOrder;
% list of networks for different Kb*Dt values
dnn.KbDts = 2.^(0:10)' * 1E-6;
files = [];
fileName = './networks2/fcPCArelaxN256Dt';
for ifl = 1 : numel(dnn.KbDts)
  files{ifl} = [fileName num2str(dnn.KbDts(ifl)) 'Kb1nModes'];
end

% current flow's Kb*Dt
flowKbDt = dnn.kappa * dnn.dt;

[~, idx] = min(abs(dnn.KbDts-flowKbDt));
if idx == 1 || idx == numel(dnn.KbDts)
  disp('Extrapolation needed for the given bending stiffness and Dt, stop!')
  pause
else
  % Choose 5 networks and we do 5th order Lagrange Interpolation 
  if idx == numel(dnn.KbDts)-1
    whichNets = (idx-(interpOrder-2):idx+1);
  elseif idx == 2
    whichNets = (idx-1:idx+(interpOrder-2));
  else
    whichNets = (idx-(interpOrder-1)/2:idx+(interpOrder-1)/2);
  end
%   whichNets = (1:numel(dnn.KbDts))';
end
% Load the networks needed for Lagrange interpolation
for k = 1 : numel(whichNets)
  if 1%dnn.KbDts(k) >= 5E-4
  load([files{whichNets(k)} '1to16_fcXlarge_tstepFCNet_flow.mat'])
  else
  load([files{whichNets(k)} '1to16_fcMedium_tstepFCNet_flow.mat'])    
  end
  dnn.bendNets{k,1} = net; 
  dnn.muChan_bend(k,1) = muChan1; 
  dnn.sdevChan_bend(k,1) = sdevChan1; 
  dnn.scale_bend(k,1) = scale; 
  dnn.offset_bend(k,1) = offset;
  if dnn.nCompRelax > 16
  if 1%dnn.KbDts(k) >= 5E-4
  load([files{whichNets(k)} '17to32_fcXlarge_tstepFCNet_flow.mat'])
  else
  load([files{whichNets(k)} '17to32_fcMedium_tstepFCNet_flow.mat']) 
  end
  dnn.bendNets{k,2} = net; 
  dnn.muChan_bend(k,2) = muChan1; 
  dnn.sdevChan_bend(k,2) = sdevChan1; 
  dnn.scale_bend(k,2) = scale; 
  dnn.offset_bend(k,2) = offset;
  end  
end
dnn.KbDts = dnn.KbDts(whichNets);

end % loadAllPCAnets4Relax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(filename,X,sigma,time,ncountNN,ncountExact)
x = X(1:end/2,:);
y = X(end/2+1:end,:);
output = [time;ncountNN;ncountExact;x(:);y(:);sigma(:)];

fid = fopen(filename,'a');
fwrite(fid,output,'double');
fclose(fid);


end
