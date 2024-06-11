function prepareRelaxFlowMapData(dtEffect, iset, npar)
% addpath ../src/
%dtEffect = 1e-6; % effective time step size, i.e., reach to this one using the small, fixed time step size

% Which set is being processed:
%iset = 1;
kappa = 1;
oc = curve;

load ./X625K_mirrdConfigs.mat
% Xstore


% number of vesicles
nves = size(Xstore,2);
% num. points per vesicle
N = 128;

% time step size
dt = 1E-5;

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
  Xinit = Xstore(:,k);
  if numel(Xinit)/2 ~= 128
    Xinit = [interpft(Xinit(1:end/2),128); interpft(Xinit(end/2+1:end),128)]; 
  end

  % Relax the shape
  [~,area,len] = oc.geomProp(Xinit);
  it = 0; X = Xinit;
  while it * dt < dtEffect
  [Xnew,errAL] = dnn.exactlySolve(X,vinf,dt,area,len,oc,op,kappa);
  X = Xnew;
  it = it + 1;
  end
  
  % Standardize input (Xinit) and output (Xnew)
  % Input
  % Standardize if not done yet
  [Xstand, scaling, rotate, trans, sortIdx] = standardizeVesicle(Xinit,oc,[],[],[],[]);
  
  
  % standardize Xnew with Xinit standardization input
  XnewStand = standardizeVesicle(Xnew,oc,trans,rotate,scaling,sortIdx);

  [xinter,~,~] = oc.selfintersect(Xnew);
  if errAL <= 1E-1 && isempty(xinter) 
   % if error in area-length is smaller than 1E-1 and not self-intersecting,
   % add the shape to the set
   idx = idx + 1;
   XstandStore(:,idx) = Xstand;
   XnewStandStore(:,idx) = XnewStand;
  end
  centerSt = oc.getPhysicalCenter(Xstand);
  centerNewSt = oc.getPhysicalCenter(XnewStand);

  centerInt = oc.getPhysicalCenter(Xinit);
  centerSol = oc.getPhysicalCenter(Xnew);

  figure(1);clf;
  plot(Xstand(1:end/2), Xstand(end/2+1:end),'k-o','linewidth',2)
  hold on
  plot(XnewStand(1:end/2), XnewStand(end/2+1:end),'r-o','linewidth',2)
  plot(centerSt(1),centerSt(2),'ko','markersize',10,'markerfacecolor','k')
  plot(centerNewSt(1),centerNewSt(2),'ro','markersize',10,'markerfacecolor','r')
  axis equal
  title('Normalized')
  % 
  % 
  % figure(2);clf;
  % plot(Xinit(1:end/2), Xinit(end/2+1:end),'k-o','linewidth',2)
  % hold on
  % plot(Xnew(1:end/2), Xnew(end/2+1:end),'r-o','linewidth',2)
  % plot(centerInt(1),centerInt(2),'ko','markersize',10,'markerfacecolor','k')
  % plot(centerSol(1),centerSol(2),'ro','markersize',10,'markerfacecolor','r')
  % axis equal
  % title('Real')
  pause
  save someData XstandStore XnewStandStore
  disp('******************************************************************')
  disp([num2str(k-sum(nSamples(1:iset-1))) 'th vesicle of out of ' num2str(nSamples(iset)) ' vesicles is done.'])
  disp(['There are ' num2str(idx) ' samples.'])   

%   if rem(idx,100) == 0 
%     nInstances = idx; 
%     fileName = ['./output/relaxNewData/n128Dt' num2str(dtEffect) 'Relax_part' num2str(iset) '.mat'];
% %     fileName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/workspace/n256Dt' num2str(dtEffect) 'Relax100K_part' num2str(iset) '.mat'];
%     %fileName = ['/workspace/gokberk/relax1step/n256Dt' num2str(dtEffect) 'Relax100K_part' num2str(iset) '.mat']; 
%     save(fileName,'nves','N','XstandStore','nInstances','XnewStandStore','dt','idx')
%   end
end
nInstances = idx;

% fileName = ['./output/relaxNewData/n128Dt' num2str(dtEffect) 'Relax_part' num2str(iset) '.mat'];
% % fileName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/workspace/n256Dt' num2str(dtEffect) 'Relax100K_part' num2str(iset) '.mat'];
% %fileName = ['/workspace/gokberk/relax1step/n256Dt' num2str(dtEffect) 'Relax100K_part' num2str(iset) '.mat']; 
% save(fileName,'nves','N','XstandStore','nInstances','XnewStandStore','dt')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xstand, scaling, rotate, trans, sortIdx] = standardizeVesicle(Xinit,oc,trans,rotate,scaling,sortIdx)
N = numel(Xinit)/2;

if isempty(trans)
% equally distribute in arc-length and find the first point

for iter = 1 : 10
  [Xinit,~,~] = oc.redistributeArcLength(Xinit);
end

[~,~,length] = oc.geomProp(Xinit);
scaling = 1/length;
% scale down 
Xinit = 1/length * Xinit;

% this is the reference vesicle - no standardization known about it
center = oc.getPhysicalCenter(Xinit); % get the center
V = oc.getPrincAxesGivenCentroid(Xinit,center); % get the principal axis

% find rotation angle
w = [0;1]; % y-axis
theta = atan2(w(2)*V(1)-w(1)*V(2), w(1)*V(1)+w(2)*V(2));

% Rotate the vesicle (with centering/decentering)
Xrot = zeros(size(Xinit));
x = Xinit(1:end/2); y = Xinit(end/2+1:end);
Xrot(1:end/2) = (x-mean(x))*cos(theta) - (y-mean(y))*sin(theta) + mean(x);
Xrot(end/2+1:end) = (x-mean(x))*sin(theta) + (y-mean(y))*cos(theta) + mean(y);

center = oc.getPhysicalCenter(Xrot);
Xrot(1:end/2) = Xrot(1:end/2)-center(1);
Xrot(end/2+1:end) = Xrot(end/2+1:end)-center(2);

rotate = theta;
trans = -center;

% Find ordering
centerRot = oc.getPhysicalCenter(Xrot); % get the center
firstQuad = find(Xrot(1:end/2)>=0 & Xrot(end/2+1:end)>=0);
theta = atan2(Xrot(end/2+1:end),Xrot(1:end/2));
[~,idx]= min(theta(firstQuad));
sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];

% now order the points
Xstand = [Xrot(sortIdx);Xrot(sortIdx+N)];
else

% this is the vesicle standardized based on another vesicle
Xrot = zeros(size(Xinit));
x = Xinit(1:end/2); y = Xinit(end/2+1:end);
Xrot(1:end/2) = (x-mean(x))*cos(rotate) - (y-mean(y))*sin(rotate) + mean(x) + trans(1);
Xrot(end/2+1:end) = (x-mean(x))*sin(rotate) + (y-mean(y))*cos(rotate) + mean(y) + trans(2);
% scale down 
Xrot = scaling * Xrot;
% equally distribute in arc-length and find the first point
for iter = 1 : 10
  [Xrot,~,~] = oc.redistributeArcLength(Xrot);
end

% Sort the points
Xstand = [Xrot(sortIdx);Xrot(sortIdx+N)];

end





end
