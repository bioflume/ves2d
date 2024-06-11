clear; 

load initShapes.mat

addpath ../../../src/
oc = curve;
w = [0;1];

theta = pi/2;
for k = 1 : numel(Xics(1,:))
  X = Xics(:,k);
  X2 = [X(1:end/2)+0.2;X(end/2+1:end)+0.2];
  cPhys = oc.getPhysicalCenterShan(X2);
  
  disp('***************************')
  disp(['Physical Center:' num2str(cPhys')])
  
  

  Vvec = [1;0]; % at (1,1)

  Xrot = rotateVes(X,cPhys,theta);
  X2rot = rotateVes(X2,cPhys,theta);

  VvecRot0 = rotateVes(Vvec,[0;0],theta);
  VvecRotC = rotateVes(Vvec,cPhys,theta);


  figure(1); clf;
  plot(X(1:end/2),X(end/2+1:end),'k','linewidth',2)
  hold on
  plot(X2(1:end/2),X2(end/2+1:end),'b','linewidth',2)
  quiver(X2(1),X2(end/2+1),Vvec(1),Vvec(2),'r')
  axis equal
  grid

  figure(2); clf;
  plot(Xrot(1:end/2),Xrot(end/2+1:end),'k','linewidth',2)
  hold on
  plot(X2rot(1:end/2),X2rot(end/2+1:end),'b','linewidth',2)
  quiver(X2rot(1),X2rot(end/2+1),VvecRot0(1),VvecRot0(2),'r')
  quiver(X2rot(1),X2rot(end/2+1),VvecRotC(1),VvecRotC(2),'b--')
  axis equal
  grid

  pause

end

function Xrot = rotateVes(X,center,theta)

 Xrot = zeros(size(X));
 x = X(1:end/2); y = X(end/2+1:end);
 Xrot(1:end/2) = (x-center(1))*cos(theta) - (y-center(2))*sin(theta) + center(1);
 Xrot(end/2+1:end) = (x-center(1))*sin(theta) + (y-center(2))*cos(theta) + center(2);

end
  