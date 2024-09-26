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



% standardize input and form layers in stand. setup
Xstand = zeros(size(X));
nv = numel(X(1,:));
scaling = zeros(nv,1);
rotate = zeros(nv,1);
rotCent = zeros(2,nv);
trans = zeros(2,nv);
sortIdx = zeros(Nnet,nv);

vesicle = capsules(X,[],[],1,1,0);
h = vesicle.length/vesicle.N;
h = h(1);

for k = 1 : nv
  [Xstand(:,k),scaling(k),rotate(k),rotCent(:,k),trans(:,k),sortIdx(:,k)] = standardizationStep(X(:,k),128);
end


% destand ves (so equilidistributed points)
Xback = zeros(size(X));
outerLayer = zeros(size(X));
interLayer = zeros(size(X));
for k = 1 : nv
   [Xstand(:,k),scaling(k),rotate(k),rotCent(:,k),trans(:,k),sortIdx(:,k)] = standardizationStep(X(:,k),128);
   Xback(:,k) = destandardize(Xstand(:,k),trans(:,k),rotate(k),rotCent(:,k),scaling(k),sortIdx(:,k));
   [~,tang] = oc.diffProp(Xback(:,k));
   nx = tang(end/2+1:end);
   ny = -tang(1:end/2);
   outerLayer(:,k) = [Xback(1:end/2,k)+h*nx; Xback(end/2+1:end,k) + h*ny];
   interLayer(:,k) = [Xback(1:end/2,k)+h/2*nx; Xback(end/2+1:end,k) + h/2*ny];
end

vesicle = capsules(Xback,[],[],1,1,0);
NearV2V = vesicle.getZone([],1);

zone = NearV2V.zone;
dist = NearV2V.dist;
nearest = NearV2V.nearest;
icp = NearV2V.icp;
argnear = NearV2V.argnear;
beta = 1.1;
hves = vesicle.length/vesicle.N;

J = find(zone{1}(:,2) == 1);
i = 4;
XLag = zeros(2,6);
nx = (Xback(J(i),2) - nearest{1}(J(i),2))/dist{1}(J(i),2);
ny = (Xback(J(i)+128,2) - nearest{1}(J(i)+128,2))/dist{1}(J(i),2);
XLag(1,:) = nearest{1}(J(i),2) + beta*hves*nx*(1:6);
XLag(2,:) = nearest{1}(J(i)+128,2) + beta*hves*ny*(1:6);


%%

figure(1); clf;
plot([outerLayer(1:end/2,1);outerLayer(1,1)],[outerLayer(end/2+1:end,1);outerLayer(end/2+1,1)],'Color',[99,99,99]/255,'linewidth',2)
hold on
hFill = fill([outerLayer(1:end/2,1);outerLayer(1,1)],[outerLayer(end/2+1:end,1);outerLayer(end/2+1,1)], [99,99,99]/255);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', [99,99,99]/255);

plot([Xback(1:end/2,1);Xback(1,1)],[Xback(end/2+1:end,1);Xback(end/2+1,1)], 'Color',[165,15,21]/255,'linewidth',1.5)
hold on
hFill = fill([Xback(1:end/2,1);Xback(1,1)],[Xback(end/2+1:end,1);Xback(end/2+1,1)], [165,15,21]/255);
hFill.FaceAlpha = 1;
set(hFill,'EdgeColor', [165,15,21]/255);

% plot([Xback(1:end/2,2);Xback(1,2)],[Xback(end/2+1:end,2);Xback(end/2+1,2)], 'Color',[5,113,176]/255,'linewidth',1.5)
hFill = fill([Xback(1:end/2,2);Xback(1,2)],[Xback(end/2+1:end,2);Xback(end/2+1,2)], [5,113,176]/255);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', [5,113,176]/255);

% p1 = scatter(Xback(1:end/2,1),Xback(end/2+1:end,1),40,'^');
% p1.MarkerFaceColor = 'w';
% p1.MarkerEdgeColor = [202,0,32]/255;
% p1.MarkerFaceAlpha = 1;
% p1.MarkerEdgeAlpha = 1;

% p2 = scatter(Xback(1:end/2,2),Xback(end/2+1:end,2),40,'^');
% p2.MarkerFaceColor = [5,113,176]/255;
% p2.MarkerEdgeColor = [5,113,176]/255;
% p2.MarkerFaceAlpha = 0;
% p2.MarkerEdgeAlpha = 1;


p4 = scatter(Xback(J,2),Xback(J+128,2),30,'o');
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


figure(2); clf;
plot([outerLayer(1:end/2,1);outerLayer(1,1)],[outerLayer(end/2+1:end,1);outerLayer(end/2+1,1)],'Color',[99,99,99]/255,'linewidth',4)
hold on
plot([Xback(1:end/2,1);Xback(1,1)],[Xback(end/2+1:end,1);Xback(end/2+1,1)], 'Color',[165,15,21]/255,'linewidth',4)
plot([Xback(1:end/2,2);Xback(1,2)],[Xback(end/2+1:end,2);Xback(end/2+1,2)], 'Color',[5,113,176]/255,'linewidth',4)
p4 = scatter(Xback(1:end/2,2),Xback(end/2+1:end,2),130,'o');
p4.MarkerFaceColor = 'w';
p4.MarkerEdgeColor = [5,113,176]/255;
p4.MarkerFaceAlpha = 1;
p4.MarkerEdgeAlpha = 1;

p4 = scatter(Xback(J,2),Xback(J+128,2),130,'o');
p4.MarkerFaceColor = [5,113,176]/255;
p4.MarkerEdgeColor = [5,113,176]/255;
p4.MarkerFaceAlpha = 0.5;
p4.MarkerEdgeAlpha = 1;

p4 = scatter(Xback(29,2),Xback(29+128,2),140,'o');
p4.MarkerFaceColor = [5,113,176]/255;
p4.LineWidth = 1;
p4.MarkerEdgeColor = [5,113,176]/255;
p4.MarkerFaceAlpha = 1;
p4.MarkerEdgeAlpha = 1;


p1 = scatter(Xback(1:end/2,1),Xback(end/2+1:end,1),130,'^');
p1.MarkerFaceColor = 'w';
p1.MarkerEdgeColor = [165,15,21]/255;
p1.MarkerFaceAlpha = 1;
p1.MarkerEdgeAlpha = 1;

p2 = scatter(Xback(1,1),Xback(end/2+1,1),150,'^');
p2.MarkerFaceColor = [165,15,21]/255;
p2.MarkerEdgeColor = [165,15,21]/255;
p2.MarkerFaceAlpha = 1;
p2.MarkerEdgeAlpha = 1;


p3 = scatter(nearest{1}(J(4),2), nearest{1}(J(4)+128,2), 150, 's');
p3.MarkerFaceColor = 'k';
p3.MarkerEdgeColor = 'k';
p3.MarkerFaceAlpha = 1;
p3.MarkerEdgeAlpha = 1;


p3 = scatter(XLag(1,:), XLag(2,:), 150, 's');
p3.MarkerFaceColor = 'k';
p3.MarkerEdgeColor = 'k';
p3.MarkerFaceAlpha = 1;
p3.MarkerEdgeAlpha = 1;


axis equal

xlim([-0.1 0])
ylim([0 0.09])



set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')

% 

figure(3);clf;
plot([outerLayer(1:end/2,1);outerLayer(1,1)],[outerLayer(end/2+1:end,1);outerLayer(end/2+1,1)],'Color',[99,99,99]/255,'linewidth',4)
hold on
plot([interLayer(1:end/2,1);interLayer(1,1)],[interLayer(end/2+1:end,1);interLayer(end/2+1,1)],'Color',[204,204,204]/255,'linewidth',4)
plot([Xback(1:end/2,1);Xback(1,1)],[Xback(end/2+1:end,1);Xback(end/2+1,1)], 'Color',[165,15,21]/255,'linewidth',4)

p1 = scatter(Xback(1:end/2,1),Xback(end/2+1:end,1),100,'^');
p1.MarkerFaceColor = [165,15,21]/255;
p1.MarkerEdgeColor = [165,15,21]/255;
p1.MarkerFaceAlpha = 1;
p1.MarkerEdgeAlpha = 1;

p1 = scatter(outerLayer(1:end/2,1),outerLayer(end/2+1:end,1),100,'^');
p1.MarkerFaceColor = [99,99,99]/255;
p1.MarkerEdgeColor = [99,99,99]/255;
p1.MarkerFaceAlpha = 1;
p1.MarkerEdgeAlpha = 1;

p1 = scatter(interLayer(1:end/2,1),interLayer(end/2+1:end,1),100,'^');
p1.MarkerFaceColor = [204,204,204]/255;
p1.MarkerEdgeColor = [204,204,204]/255;
p1.MarkerFaceAlpha = 1;
p1.MarkerEdgeAlpha = 1;


plot([Xback(1:end/2,2);Xback(1,2)],[Xback(end/2+1:end,2);Xback(end/2+1,2)], 'Color',[5,113,176]/255,'linewidth',4)
p4 = scatter(Xback(1:end/2,2),Xback(end/2+1:end,2),100,'o');
p4.MarkerFaceColor = 'w';
p4.MarkerEdgeColor = [5,113,176]/255;
p4.MarkerFaceAlpha = 1;
p4.MarkerEdgeAlpha = 1;

p4 = scatter(Xback(J,2),Xback(J+128,2),100,'o');
p4.MarkerFaceColor = [5,113,176]/255;
p4.MarkerEdgeColor = [5,113,176]/255;
p4.MarkerFaceAlpha = 0.5;
p4.MarkerEdgeAlpha = 1;

p4 = scatter(Xback(29,2),Xback(29+128,2),100,'o');
p4.MarkerFaceColor = [5,113,176]/255;
p4.LineWidth = 1;
p4.MarkerEdgeColor = [5,113,176]/255;
p4.MarkerFaceAlpha = 1;
p4.MarkerEdgeAlpha = 1;

axis equal

xlim([-0.1 0])
ylim([0 0.09])

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')


figure(1); 
ax = gca;
exportgraphics(ax,'~/Desktop/nearSchematic1.png','Resolution',300)


figure(2); 
ax = gca;
exportgraphics(ax,'~/Desktop/nearSchematic2.png','Resolution',300)


figure(3); 
ax = gca;
exportgraphics(ax,'~/Desktop/nearSchematic3.png','Resolution',300)


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



