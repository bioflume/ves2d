clear;

% load someVesConfigs.mat
load someSymmVes.mat

addpath ../../src/
oc = curve;
w = [0;1]; % y-axis

for k = 1 : numel(Xics(1,:))
  X1 = Xics(:,k);
  X2 = Xics2(:,k);

  center1 = oc.getPhysicalCenterShan(X1);
  V = oc.getPrincAxesGivenCentroid(X1,center1);
  disp(['Center1 from data set: ' num2str(center1')])
  
  theta = atan2(w(2)*V(1)-w(1)*V(2), w(1)*V(1)+w(2)*V(2));
  disp(['Angle in the data set: ' num2str(theta)])

  % Rotate the vesicle (with centering/decentering)
  Xrot1 = rotateVes(X1,center1,theta);
  center2 = oc.getPhysicalCenterShan(X2);
  Xrot2 = rotateVes(X2,center2,theta);
  
  disp(['Center2 from data set: ' num2str(center2')])

  % Center after rotating
  center1rot = oc.getPhysicalCenterShan(Xrot1);
  V = oc.getPrincAxesGivenCentroid(Xrot1,center1rot);
  disp(['Center1 after rotation: ' num2str(center1rot')])
  
  thetaRot = atan2(w(2)*V(1)-w(1)*V(2), w(1)*V(1)+w(2)*V(2));
  disp(['Angle1 after rot.: ' num2str(thetaRot)])

  center2rot = oc.getPhysicalCenterShan(Xrot2);
  disp(['Center2 after rotation: ' num2str(center2rot')])

  disp('*************************')
  disp(' ')
  figure(1);clf;
  plot(X1(1:end/2),X1(end/2+1:end),'linewidth',2)
  hold on
  plot(Xrot1(1:end/2),Xrot1(end/2+1:end),'linewidth',2)
  legend('Before','After')
  axis equal

  figure(2);clf;
  plot(X2(1:end/2),X2(end/2+1:end),'linewidth',2)
  hold on
  plot(Xrot2(1:end/2),Xrot2(end/2+1:end),'linewidth',2)
  legend('Before','After')
  axis equal

  pause

end


function Xrot = rotateVes(X,center,theta)

 Xrot = zeros(size(X));
 x = X(1:end/2); y = X(end/2+1:end);
 Xrot(1:end/2) = (x-center(1))*cos(theta) - (y-center(2))*sin(theta) + center(1);
 Xrot(end/2+1:end) = (x-center(1))*sin(theta) + (y-center(2))*cos(theta) + center(2);

end
