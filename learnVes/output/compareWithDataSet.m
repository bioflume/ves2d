% clear; clc;
% addpath ../../src/
% 
% fileIDs = [ 58,  52,  89,  55,  18,  59,   2,  36,   3,  81,  73,  99,  76,        16,  73,   9,  78,  82,  41,  42,  11,  99,   3,  49,  46,  30,        62,  58,   5,  86,  85,  57, 100,  64,  51,  57,  85,  73,  90,        82,  72,  88,  32,   2,  69,  63,  99,  10,  17,  77,  28,  66,        29,  91,  67,  19,  50,  25,  34,  96,  64,  84,  44,  29,  76,       100,  57,  24,  37,   3,  66,  15,  36,  21,  30,  39,  40,  42,        61,   3, 100,  94,  41,  51,  19,  47,  60,  27,  55,  21,  35,         8,  29,  17,  98,  11,  10,   4,  57,  98];
% vesIDs = [14, 12, 31,  9,  1,  1,  9, 27,  1,  9,  1,  2,  2,  3,  2, 29,  3,       13,  2, 31,  2,  1,  2, 14,  9, 31, 13,  4, 13,  1, 27,  3, 13,  9,        3, 31, 30, 29,  4, 27, 13, 27, 14, 15, 13,  4, 29, 13, 31, 13, 27,        1,  3,  3, 17, 13,  7, 17, 29,  3,  7, 13, 31,  4, 13,  1,  0,  4,        3, 29, 27, 17,  3, 27, 13, 23, 27,  0, 13,  7,  9, 14, 30, 27, 31,       14, 31,  1, 13, 30,  0, 13, 19,  1, 14, 30, 23,  2, 29, 27] + 1;
% load ../dataPrepCodes/advectionNetInputX.mat
%%
for k = 2 : 10 %numel(fileIDs)
  k
  load(['~/Desktop/randFromGT/comparisonX' num2str(vesIDs(k)) '.mat'])
  [Xstd,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(X(:,vesIDs(k)),128);
  % for nv = 1 : size(XstandStore,2)
  %   diff(nv,1) = mean(sqrt((Xstd(1:end/2)-XstandStore(1:end/2,nv)).^2 + (Xstd(end/2+1:end)-XstandStore(end/2+1:end,nv)).^2))./...
  %       mean(sqrt(Xstd(1:end/2).^2 + Xstd(end/2+1:end).^2));
  % end
  % diffStore{k} = diff;
  % [maxDiff(k,1),maxID(k,1)] = max(diff);
  % [minDiff(k,1),minID(k,1)] = min(diff);
  % % maxID(k,1) = 
  [maxDiff(k,1), maxID(k,1), minDiff(k,1), minID(k,1)] = hausdorfDistance(Xstd,XstandStore);

end
%%
oc = curve;
for k = 1 : 10

figure(1);clf;
[Xstd,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(X(:,vesIDs(k)),128);
subplot(1,2,1)
plot(Xstd(1:end/2),Xstd(end/2+1:end),'k-o','linewidth',2)
hold on
plot(XstandStore(1:end/2,minID(k)),XstandStore(end/2+1:end,minID(k)),'r-o','linewidth',2)
centShape = oc.getPhysicalCenterShan(Xstd)
centData = oc.getPhysicalCenterShan(XstandStore(:,minID(k)))
axis equal
legend('Shape','Closest in data')
title(['Diff = ' num2str(minDiff(k))])


subplot(1,2,2)
histogram(diffStore{k})
axis square
title('Difference distribution')
% axis equal
% 
% subplot(1,3,3)
% plot(XstandStore(1:end/2,maxID(k)),XstandStore(end/2+1:end,maxID(k)),'k-o','linewidth',2)
% axis equal
% title(['Most different in data set, L2 = ' num2str(maxDiff(k))])
% axis equal

pause


end



%%
function [maxErr, maxID, minErr, minID] = hausdorfDistance(X1,X2)
N = numel(X1)/2;
nv = size(X2,2);

% Closest points on X2 to the points on X1
for j = 1 : nv
d1to2 = zeros(N,1);
for in = 1 : N
    d1to2(in) = min(((X1(in)-X2(1:end/2,j)).^2 + (X1(in+N)-X2(end/2+1:end,j)).^2).^0.5)/...
      sqrt(X1(in).^2+X1(in+N).^2);
end
hausErr(j,1) = max(d1to2);
end
[maxErr, maxID] = max(hausErr);
[minErr, minID] = min(hausErr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(Xin,Nnet)
oc = curve;
N = numel(Xin)/2;
if Nnet ~= N
  Xin = [interpft(Xin(1:end/2),Nnet);interpft(Xin(end/2+1:end),Nnet)];    
end

% Equally distribute points in arc-length
for iter = 1 : 10
  [Xin,~,~] = oc.redistributeArcLength(Xin);
end


X = Xin;
[trans,rotate,rotCent,scaling,sortIdx] = referenceValues(X);

% Fix misalignment in center and angle due to reparametrization
% X = oc.alignCenterAngle(Xin,X);

% standardize angle, center, scaling and point order

X = standardize(X,trans,rotate,rotCent,scaling,sortIdx);
end % standardizationStep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XrotSort = standardize(X,translation,rotation,rotCent,scaling,sortIdx)
N = numel(sortIdx);

% translate, rotate and scale configuration

Xrotated = rotationOperator(X,rotation,rotCent);   
Xrotated = translateOp(Xrotated,translation);

% now order the points
XrotSort = [Xrotated(sortIdx);Xrotated(sortIdx+N)];

XrotSort = scaling*XrotSort;

end % standardize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [translation,rotation,rotCent,scaling,sortIdx] = referenceValues(Xref)
oc = curve;
N = numel(Xref)/2;

% find translation, rotation and scaling
center = oc.getPhysicalCenterShan(Xref);
V = oc.getPrincAxesGivenCentroid(Xref,center);
% % find rotation angle
w = [0;1]; % y-axis
rotation = atan2(w(2)*V(1)-w(1)*V(2), w(1)*V(1)+w(2)*V(2));


% translation = [-mean(Xref(1:end/2));-mean(Xref(end/2+1:end))];
% rotation = pi/2-oc.getIncAngle2(Xref);
       
% find the ordering of the points
rotCent = center;
Xref = rotationOperator(Xref, rotation, center);
center = oc.getPhysicalCenterShan(Xref);
translation = -center;

Xref = translateOp(Xref, translation);

firstQuad = find(Xref(1:end/2)>=0 & Xref(end/2+1:end)>=0);
theta = atan2(Xref(end/2+1:end),Xref(1:end/2));
[~,idx]= min(theta(firstQuad));
sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];

% amount of scaling
[~,~,length] = oc.geomProp(Xref);
scaling = 1/length;
end % referenceValues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xrot = rotationOperator(X,theta, rotCent)
% Get x-y coordinates
Xrot = zeros(size(X));
x = X(1:end/2)-rotCent(1); y = X(end/2+1:end)-rotCent(2);

% Rotated shape
xrot = (x)*cos(theta) - (y)*sin(theta);
yrot = (x)*sin(theta) + (y)*cos(theta);

Xrot(1:end/2) = xrot+rotCent(1);
Xrot(end/2+1:end) = yrot+rotCent(2);
end % rotationOperator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateOp(X,transXY)
Xnew = zeros(size(X));
Xnew(1:end/2) = X(1:end/2)+transXY(1);
Xnew(end/2+1:end) = X(end/2+1:end)+transXY(2);
end  % translateOp  