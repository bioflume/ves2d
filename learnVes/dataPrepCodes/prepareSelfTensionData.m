function prepareSelfTensionData(iset,npar)
%load /workspace/gokberk/relax1step/N96DataSets/n96Relax100KAllDataSet.mat
%load /workspace/gokberk/relax1step/only1timeStepPCAData.mat
%load /workspace/gokberk/relax1step/n256Dt1E4Kb1RelaxAllDataSet.mat
%clear XnewStandStore;

load ./workingOnDataSet/X156Kconfigs.mat

% Xstore Loaded
% clear sum of the data to have some space

%clear XnewCoeffs; clear XnewRec; clear XnewStandStore;
%clear XoldCoeffs; clear Xrec; clear colMeans;
%clear errInNew; clear evects;

addpath ../src/
oc = curve;
% num. points
N = 128;
op = poten(N);

% store aligned shapes
XstandStore = [];
nInstances = size(Xstore,2);
% 


nSamples = ones(npar,1)*floor(nInstances/npar);
nSamples(end) = nInstances-(npar-1)*nSamples(1);
dnn = dnnTools;

selfTenStore = [];

idx = 1;
for ives = sum(nSamples(1:iset-1))+1:sum(nSamples(1:iset))
  disp(['Vesicle #' num2str(idx) ' out of ' num2str(nSamples(iset)) ' being processed...'])
  tstart = tic;
  % change the resolution
  Xinit = [interpft(Xstore(1:end/2,ives),N); interpft(Xstore(end/2+1:end,ives),N)]; 
  
   % Build vesicle
  vesicle = capsules(Xinit,[],[],1,1,0);

  % derivatives
  [Ben,Ten,Div] = vesicle.computeDerivs;
  
  % SLP
  G = op.stokesSLmatrix(vesicle);
  
  % Build M and Multiply with Basis (predict Z11 on 48 points)
  M = ((Div*G*Ten)\eye(vesicle.N))*Div;
  
  % Compute the action of M on G*(-Ben)
  tension = M*G*(-Ben)*Xinit; 
  
  % standardize
  [Xstand,scaling,rotate,trans,sortIdx] = dnn.standardizationStep(Xinit,oc);

  XstandStore(:,idx) = Xinit;

  selfTenStore(:,idx) = scaling*tension(sortIdx);
  
  idx = idx + 1;
  tend = toc(tstart);
  disp(['took ' num2str(tend) ' seconds'])
end


fileName = ['E:\selfTensionData\selfTensionData_' num2str(iset) '.mat']; 
nsampInSet = nSamples(iset);
save(fileName,'nInstances','nsampInSet','XstandStore','selfTenStore','-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xstand, scaling, rotate, rotCent, trans, sortIdx] = standardizeVesicle(Xinit,oc,trans,rotate,scaling,sortIdx,rotCent)
N = numel(Xinit)/2;

% equally distribute in arc-length and find the first point
for iter = 1 : 10
  [Xinit,~,~] = oc.redistributeArcLength(Xinit);
end


if isempty(trans)
[~,~,length] = oc.geomProp(Xinit);
scaling = 1/length;

% this is the reference vesicle - no standardization known about it
center = oc.getPhysicalCenter(Xinit); % get the center
V = oc.getPrincAxesGivenCentroid(Xinit,center); % get the principal axis
rotCent = center;
% find rotation angle
w = [0;1]; % y-axis
theta = atan2(w(2)*V(1)-w(1)*V(2), w(1)*V(1)+w(2)*V(2));

% Rotate the vesicle (with centering/decentering)
Xrot = zeros(size(Xinit));
x = Xinit(1:end/2); y = Xinit(end/2+1:end);
Xrot(1:end/2) = (x-center(1))*cos(theta) - (y-center(2))*sin(theta) + center(1) ;
Xrot(end/2+1:end) = (x-center(1))*sin(theta) + (y-center(2))*cos(theta) + center(2) ;

center = oc.getPhysicalCenter(Xrot); % get the center
Xrot(1:end/2) = Xrot(1:end/2)-center(1);
Xrot(end/2+1:end) = Xrot(end/2+1:end)-center(2);

rotate = theta;
trans = -center;

% scale down 
Xrot = 1/length * Xrot;

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
Xrot(1:end/2) = (x-rotCent(1))*cos(rotate) - (y-rotCent(2))*sin(rotate) + rotCent(1) + trans(1);
Xrot(end/2+1:end) = (x-rotCent(1))*sin(rotate) + (y-rotCent(2))*cos(rotate) + rotCent(2) + trans(2);
% scale down 
Xrot = scaling * Xrot;

% Sort the points
Xstand = [Xrot(sortIdx);Xrot(sortIdx+N)];

end


