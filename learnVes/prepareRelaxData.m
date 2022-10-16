clear; clc;
addpath ../src/

% Which set is being processed:
iset = 4
kappa = 1;
oc = curve;
% set of vesicle shapes
load /workspace/gokberk/X100KinitShapes.mat
% Xstore

istandard = false; % if the shapes are already standardized

% number of vesicles
nves = size(Xstore,2);
% num. points per vesicle
N = 256;

% time step size
dt = 2.5E-3;

% store aligned shapes
XstandStore = [];

% store right hand sides
XnewStandStore = [];

% poten operator
op = poten(N);
dnn = dnnTools;

% Separate into four sets, so we can quickly process
% number of samples in the sets
nSamples = [floor(nves/4); floor(nves/4); floor(nves/4); nves-3*floor(nves/4)];


% background velocity, 0 for relaxation
vinf = @(X) zeros(size(X));

idx = 0;
for k = sum(nSamples(1:iset-1))+1 : sum(nSamples(1:iset))

  % change the resolution
  X = [interpft(Xstore(1:end/2,k),N);...
        interpft(Xstore(end/2+1:end,k),N)]; 
  
  if ~istandard
    [X,scaling,rotate,trans,sortIdx] = dnn.standardizationStep(X,oc);
  end

  [~,area,len] = oc.geomProp(X);
  [Xnew,errAL] = dnn.exactlySolve(X,vinf,dt,area,len,oc,op,kappa);
  [xinter,~,~] = oc.selfintersect(Xnew);

  if errAL <= 1E-1 && isempty(xinter) 
  % if error in area-length is smaller than 1E-1 and not self-intersecting,
  % add the shape to the set
   idx = idx + 1;
   XstandStore(:,idx) = X;
   XnewStandStore(:,idx) = Xnew;
  end

  disp('******************************************************************')
  disp([num2str(k-sum(nSamples(1:iset-1))) 'th vesicle of out of ' num2str(nSamples(iset)) ' vesicles is done.'])
  disp(['There are ' num2str(idx) ' samples.'])   
end
nInstances = idx;

fileName = ['/workspace/gokberk/relax1step/n256Dt2p5E3Kb1Relax100K_part' num2str(iset) '.mat']; 
save(fileName,'nves','N','XstandStore','nInstances','XnewStandStore','dt')


