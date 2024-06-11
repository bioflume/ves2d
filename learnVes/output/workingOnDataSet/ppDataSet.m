% load newStand625K_n128Dt1e-05RelaxDataSet.mat

cx = zeros(nInstances,1); 
cy = cx;
cxN = cx;
cyN = cy;

cxP = cx;
cyP = cx;
cxPN = cx;
cyPN = cy;

oc = curve;

for k = 1 : nInstances
  cx(k,1) = mean(XstandStore(1:end/2,k));
  cy(k,1) = mean(XstandStore(end/2+1:end,k));

  cxN(k,1) = mean(XnewStandStore(1:end/2,k));
  cyN(k,1) = mean(XnewStandStore(end/2+1:end,k));

  center = oc.getPhysicalCenterShan(XstandStore(:,k));
  centerN = oc.getPhysicalCenterShan(XnewStandStore(:,k));

  cxP(k,1) = center(1);
  cyP(k,1) = center(2);

  cxPN(k,1) = centerN(1);
  cyPN(k,1) = centerN(2);
  
  Xinit = XstandStore(:,k); Xnew = XnewStandStore(:,k);

  [XstandStore(:,k),scaling,rotate,trans,sortIdx] = standardizeVesicle(Xinit,oc,[],[],[],[]);
  XnewStandStore(:,k) = standardizeVesicle(Xnew,oc,trans,rotate,scaling,sortIdx);

  % XstandStore(:,k) = [XstandStore(1:end/2,k)-cxP(k); XstandStore(end/2+1:end,k)-cyP(k)];
  % XnewStandStore(:,k) = [XnewStandStore(1:end/2,k)-cxP(k); XnewStandStore(end/2+1:end,k)-cyP(k)];

  center = oc.getPhysicalCenterShan(XstandStore(:,k));
  centerN = oc.getPhysicalCenterShan(XnewStandStore(:,k));

  cxP(k,1) = center(1);
  cyP(k,1) = center(2);

  cxPN(k,1) = centerN(1);
  cyPN(k,1) = centerN(2);
  
end

dx = cxN-cx;
dy = cyN-cy;
dxP = cxPN-cxP;
dyP = cyPN-cyP;

save newStand625K_n128Dt1e-05RelaxDataSet_IT3.mat XnewStandStore XstandStore nInstances N kappa dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xstand, scaling, rotate, trans, sortIdx] = standardizeVesicle(Xinit,oc,trans,rotate,scaling,sortIdx)
N = numel(Xinit)/2;

if isempty(trans)
% equally distribute in arc-length and find the first point

% for iter = 1 : 10
%   [Xinit,~,~] = oc.redistributeArcLength(Xinit);
% end

% this is the reference vesicle - no standardization known about it
center = oc.getPhysicalCenterShan(Xinit); % get the center
V = oc.getPrincAxesGivenCentroid(Xinit,center); % get the principal axis

% find rotation angle
w = [0;1]; % y-axis
theta = atan2(w(2)*V(1)-w(1)*V(2), w(1)*V(1)+w(2)*V(2));

% Rotate the vesicle (with centering/decentering)
Xrot = zeros(size(Xinit));
x = Xinit(1:end/2); y = Xinit(end/2+1:end);
Xrot(1:end/2) = (x-mean(x))*cos(theta) - (y-mean(y))*sin(theta) + mean(x);
Xrot(end/2+1:end) = (x-mean(x))*sin(theta) + (y-mean(y))*cos(theta) + mean(y);

center = oc.getPhysicalCenterShan(Xrot);
Xrot(1:end/2) = Xrot(1:end/2)-center(1);
Xrot(end/2+1:end) = Xrot(end/2+1:end)-center(2);

rotate = theta;
trans = -center;

% Find ordering
firstQuad = find(Xrot(1:end/2)>=0 & Xrot(end/2+1:end)>=0);
theta = atan2(Xrot(end/2+1:end),Xrot(1:end/2));
[~,idx]= min(theta(firstQuad));
sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];

% now order the points
Xstand = [Xrot(sortIdx);Xrot(sortIdx+N)];

[~,~,length] = oc.geomProp(Xstand);
scaling = 1/length;
% scale down 
Xstand = 1/length * Xstand;
else

% for iter = 1 : 10
%   [Xinit,~,~] = oc.redistributeArcLength(Xinit);
% end

% this is the vesicle standardized based on another vesicle
Xrot = zeros(size(Xinit));
x = Xinit(1:end/2); y = Xinit(end/2+1:end);
Xrot(1:end/2) = (x-mean(x))*cos(rotate) - (y-mean(y))*sin(rotate) + mean(x) + trans(1);
Xrot(end/2+1:end) = (x-mean(x))*sin(rotate) + (y-mean(y))*cos(rotate) + mean(y) + trans(2);
% scale down 
Xrot = scaling * Xrot;

% Sort the points
Xstand = [Xrot(sortIdx);Xrot(sortIdx+N)];

end





end