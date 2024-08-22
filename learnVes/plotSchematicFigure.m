clear; clc;

addpath ../src/
oc = curve;

load shearInstants.mat
Nnet = 128;
X = [xs(:,:,4); ys(:,:,4)];
x1 = xs(:,1,4);
y1 = ys(:,1,4);
x2 = xs(:,2,4);
y2 = ys(:,2,4);


maxLayerDist = sqrt(1/128);
nlayers = 3;
dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist;

% standardize input and form layers in stand. setup
Xstand = zeros(size(X));
nv = numel(X(1,:));
scaling = zeros(nv,1);
rotate = zeros(nv,1);
rotCent = zeros(2,nv);
trans = zeros(2,nv);
sortIdx = zeros(Nnet,nv);

tracersX = zeros(2*Nnet,3,nv);
for k = 1 : nv
  [Xstand(:,k),scaling(k),rotate(k),rotCent(:,k),trans(:,k),sortIdx(:,k)] = standardizationStep(X(:,k),128);
  [~,tang] = oc.diffProp(Xstand(:,k));
  nx = tang(Nnet+1:2*Nnet);
  ny = -tang(1:Nnet);

  tracersX(:,1,k) = Xstand(:,k);
  for il = 2 : nlayers 
    tracersX(:,il,k) = [Xstand(1:end/2,k)+nx*dlayer(il); Xstand(end/2+1:end,k)+ny*dlayer(il)];
  end
end


% destand layers
xlayers = zeros(128,3,nv);
ylayers = zeros(128,3,nv);
for k = 1 : nv
  for il = 1 : 3
     Xl = destandardize(tracersX(:,il,k),trans(:,k),rotate(k),rotCent(:,k),scaling(k),sortIdx(:,k));
     xlayers(:,il,k) = Xl(1:end/2);
     ylayers(:,il,k) = Xl(end/2+1:end);
  end
end


% do ray casting
Xlarge = zeros(256,2);
Xlarge(:,1) = [xlayers(:,3,1); ylayers(:,3,1)];
Xlarge(:,2) = [xlayers(:,3,2); ylayers(:,3,2)];

iCallNear = zeros(nv,1);
for j = 1 : nv
  K = [(1:j-1) (j+1:nv)];
  
  S = zeros(256,1);
  
  S(1:2:end) = Xlarge(1:end/2,j);
  S(2:2:end) = Xlarge(end/2+1:end,j);
  for k = K
    queryX{k} = []; % k's points in j's near-field
    idsInStore{k} = [];

    % also store neighbor vesicles
    nearVesIds{k} = [];

    cnt = 1; 
    for p = 1 : 128
      flag = rayCasting([X(p,k);X(p+128,k)],S);  
      if flag
        idsInStore{k}(cnt,1) = p;
        % points where we need interpolation  
        queryX{k}(1,cnt) = X(p,k);
        queryX{k}(2,cnt) = X(p+128,k);
        nearVesIds{k}(cnt,1) = j; 
        cnt = cnt + 1;
        iCallNear(k) = 1;    
      end
    end
  end
end

figure(1); clf;
plot([x1;x1(1)],[y1;y1(1)],'Color',[202,0,32]/255,'linewidth',1.5)
hold on
hFill = fill([x1;x1(1)],[y1;y1(1)], [202,0,32]/255);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', [202,0,32]/255);

plot([x2;x2(1)],[y2;y2(1)], 'Color',[5,113,176]/255,'linewidth',1.5)
hFill = fill([x2;x2(1)],[y2;y2(1)], [5,113,176]/255);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', [5,113,176]/255);

xls = xlayers(:,:,1);
yls = ylayers(:,:,1);
p1 = scatter(xls(:),yls(:),15,'o');
xls = xlayers(:,:,2);
yls = ylayers(:,:,2);
p2 = scatter(xls(:),yls(:),15,'^');
p1.MarkerFaceColor = [202,0,32]/255;
p1.MarkerEdgeColor = [202,0,32]/255;
p1.MarkerFaceAlpha = 0;
p1.MarkerEdgeAlpha = 1;


p2.MarkerFaceColor = [5,113,176]/255;
p2.MarkerEdgeColor = [5,113,176]/255;
p2.MarkerFaceAlpha = 0;
p2.MarkerEdgeAlpha = 1;



p3 = scatter(queryX{1}(1,:),queryX{1}(2,:),40,'o');
p3.MarkerFaceColor = [202,0,32]/255;
p3.MarkerEdgeColor = [202,0,32]/255;
p3.MarkerFaceAlpha = 1;
p3.MarkerEdgeAlpha = 1;



p4 = scatter(queryX{2}(1,:),queryX{2}(2,:),40,'o');
p4.MarkerFaceColor = [5,113,176]/255;
p4.MarkerEdgeColor = [5,113,176]/255;
p4.MarkerFaceAlpha = 1;
p4.MarkerEdgeAlpha = 1;

axis equal

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);
        
set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')

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
function X = destandardize(XrotSort,translation,rotation,rotCent,scaling,sortIdx)

N = numel(sortIdx);    
    
% scaling back
XrotSort = XrotSort/scaling;

% change ordering back 
X = zeros(size(XrotSort));
X([sortIdx;sortIdx+N]) = XrotSort;

% take translation back
X = translateOp(X,-translation);

% take rotation back
X = rotationOperator(X,-rotation,rotCent);


end % destandardize

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



