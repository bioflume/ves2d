clear; 

% load someSymmVes.mat
load initShapes.mat

addpath ../../src/
oc = curve;
w = [0;1];

theta = pi/3;

X1 = Xics(:,1);
X2 = Xics(:,1);

% First center their physical centers
c1 = [-0.1;0.1];
c2 = [0.1; -0.1];

cP = oc.getPhysicalCenter(X1);
X1 = [X1(1:end/2)-cP(1)+c1(1);X1(end/2+1:end)-cP(2)+c1(2)];

cP = oc.getPhysicalCenter(X2);
X2 = [X2(1:end/2)-cP(1)+c2(1);X2(end/2+1:end)-cP(2)+c2(2)];

cP = oc.getPhysicalCenter(X1);
disp(['New physical CoM of X1:' num2str(cP')])

cP = oc.getPhysicalCenter(X2);
disp(['New physical CoM of X2:' num2str(cP')])
X2 = rotateVes(X2,cP,pi/6);
cP = oc.getPhysicalCenter(X2);
disp(['New physical CoM of X2 after rotation:' num2str(cP')])
X2 = [X2(1:end/2)-cP(1)+c2(1);X2(end/2+1:end)-cP(2)+c2(2)];
cP = oc.getPhysicalCenter(X2);
disp(['Physical center is fixed:' num2str(cP')])



% NOW LET'S ADJUST X1 -- THEN ADJUST X2 BASED ON THAT
cP1 = oc.getPhysicalCenter(X1);
cP2 = oc.getPhysicalCenter(X2);

d21 = cP2 - cP1;
disp(['Distance of cP2 to cP1:' num2str(d21')])
V1 = oc.getPrincAxesGivenCentroid(X1,cP1);
V2 = oc.getPrincAxesGivenCentroid(X2,cP2);

theta1 = atan2(w(2)*V1(1)-w(1)*V1(2), w(1)*V1(1)+w(2)*V1(2));
theta2 = atan2(w(2)*V2(1)-w(1)*V2(2), w(1)*V2(1)+w(2)*V2(2));

X1rot = rotateVes(X1,[mean(X1(1:end/2));mean(X1(end/2+1:end))],pi/3);
cP = oc.getPhysicalCenter(X1rot);
X1rot = [X1rot(1:end/2)-cP(1); X1rot(end/2+1:end)-cP(2)];

X2rot = rotateVes(X2,[mean(X2(1:end/2));mean(X2(end/2+1:end))],pi/3);
X2rot = [X2rot(1:end/2)-cP(1); X2rot(end/2+1:end)-cP(2)];

cP1 = oc.getPhysicalCenter(X1rot);
cP2 = oc.getPhysicalCenter(X2rot);
d21rot = cP2 - cP1;

figure(1); clf;
plot(X1(1:end/2),X1(end/2+1:end),'r','linewidth',2)
hold on
plot(X2(1:end/2),X2(end/2+1:end),'b','linewidth',2)
axis equal
legend('X(t)','X(t+dt)')
legend boxoff
ax = gca;
exportgraphics(ax,'~/Desktop/orig.png','Resolution',300)

figure(2); clf;
plot(X1rot(1:end/2),X1rot(end/2+1:end),'r','linewidth',2)
hold on
plot(X2rot(1:end/2),X2rot(end/2+1:end),'b','linewidth',2)
axis equal
legend('X(t)','X(t+dt)')
legend boxoff
ax = gca;
exportgraphics(ax,'~/Desktop/alt1.png','Resolution',300)

%% 
DX = X2 - X1;

DXrot = rotateVes(DX,[0;0],pi/3);
X2rot2 = [X1rot(1:end/2)+DXrot(1:end/2);X1rot(end/2+1:end)+DXrot(end/2+1:end)];

figure(3); clf;
plot(X1rot(1:end/2),X1rot(end/2+1:end),'r','linewidth',2)
hold on
plot(X2rot2(1:end/2),X2rot2(end/2+1:end),'b','linewidth',2)
axis equal
legend('X(t)','X(t+dt)')
legend boxoff
ax = gca;
exportgraphics(ax,'~/Desktop/alt3.png','Resolution',300)


function Xrot = rotateVes(X,center,theta)

 Xrot = zeros(size(X));
 x = X(1:end/2); y = X(end/2+1:end);
 Xrot(1:end/2) = (x-center(1))*cos(theta) - (y-center(2))*sin(theta) + center(1);
 Xrot(end/2+1:end) = (x-center(1))*sin(theta) + (y-center(2))*cos(theta) + center(2);

end
  