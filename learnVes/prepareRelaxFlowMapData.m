function prepareRelaxFlowMapData(dtEffect, iset, npar)
addpath ../src/
%dtEffect = 1e-6; % effective time step size, i.e., reach to this one using the small, fixed time step size

% Which set is being processed:
%iset = 1;
kappa = 1;
oc = curve;
% set of vesicle shapes
%load /workspace/gokberk/X100KinitShapes.mat
% load /mnt/ceph/users/gkabacaoglu/SVTGRuns/workspace/X100KinitShapes.mat
load ./necessaryMatFiles/X100KinitShapes.mat
% Xstore

istandard = false; % if the shapes are already standardized

% number of vesicles
nves = size(Xstore,2);
% num. points per vesicle
N = 128;

% time step size
dt = 1E-4;

% store aligned shapes
XstandStore = [];

% store right hand sides
XnewStandStore = [];

% poten operator
op = poten(N);
dnn = dnnTools;

% Separate into four sets, so we can quickly process
% number of samples in the sets for 1:5
%nSamples = [floor(nves/4); floor(nves/4); floor(nves/4); nves-3*floor(nves/4)];

% for 6:16
nSamples = ones(npar,1)*floor(nves/npar);
nSamples(end) = nSamples(end) + nves-sum(nSamples);

% background velocity, 0 for relaxation
vinf = @(X) zeros(size(X));

idx = 0;
for k = sum(nSamples(1:iset-1))+1 : sum(nSamples(1:iset))
  % change the resolution
  Xinit = [interpft(Xstore(1:end/2,k),N); interpft(Xstore(end/2+1:end,k),N)]; 
  % Standardize if not done yet
  if ~istandard
    [Xinit,scaling,rotate,trans,sortIdx] = dnn.standardizationStep(Xinit,oc);
  end
  [~,area,len] = oc.geomProp(Xinit);

  it = 0; X = Xinit;
  while it * dt < dtEffect
  [Xnew,errAL] = dnn.exactlySolve(X,vinf,dt,area,len,oc,op,kappa);
  X = Xnew;
  it = it + 1;
  end
  [xinter,~,~] = oc.selfintersect(Xnew);
  if errAL <= 1E-1 && isempty(xinter) 
  % if error in area-length is smaller than 1E-1 and not self-intersecting,
  % add the shape to the set
   idx = idx + 1;
   XstandStore(:,idx) = Xinit;
   XnewStandStore(:,idx) = Xnew;
  end

  disp('******************************************************************')
  disp([num2str(k-sum(nSamples(1:iset-1))) 'th vesicle of out of ' num2str(nSamples(iset)) ' vesicles is done.'])
  disp(['There are ' num2str(idx) ' samples.'])   

  if rem(idx,10) == 0 
    nInstances = idx; 
    fileName = ['/work2/03353/gokberk/frontera/relaxData/n256Dt' num2str(dtEffect) 'Relax100K_part' num2str(iset) '.mat'];
%     fileName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/workspace/n256Dt' num2str(dtEffect) 'Relax100K_part' num2str(iset) '.mat'];
    %fileName = ['/workspace/gokberk/relax1step/n256Dt' num2str(dtEffect) 'Relax100K_part' num2str(iset) '.mat']; 
    save(fileName,'nves','N','XstandStore','nInstances','XnewStandStore','dt','idx')
  end
end
nInstances = idx;

fileName = ['./relaxData/n256Dt' num2str(dtEffect) 'Relax100K_part' num2str(iset) '.mat'];
% fileName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/workspace/n256Dt' num2str(dtEffect) 'Relax100K_part' num2str(iset) '.mat'];
%fileName = ['/workspace/gokberk/relax1step/n256Dt' num2str(dtEffect) 'Relax100K_part' num2str(iset) '.mat']; 
save(fileName,'nves','N','XstandStore','nInstances','XnewStandStore','dt')
end


