clear; 

% load someSymmVes.mat
load initShapes.mat

addpath ../../src/
oc = curve;
w = [0;1];

theta = pi/3;
for k = 1 : numel(Xics(1,:))
  X = Xics(:,k);
  cPhys = oc.getPhysicalCenterShan(X);
  cGeom = [mean(X(1:end/2));mean(X(end/2+1:end))];
  disp('***************************')
  disp(['Physical Center:' num2str(cPhys')])
  disp(['Geometrical Center:' num2str(cGeom')])

  Vphs = oc.getPrincAxesGivenCentroid(X,cPhys);
  Vgeo = oc.getPrincAxesGivenCentroid(X,cGeom);
   
  disp(['PA w.r.t. Physical Center:' num2str(Vphs')])
  disp(['PA w.r.t. Geometrical Center:' num2str(Vgeo')])

  thetaPhs = atan2(w(2)*Vphs(1)-w(1)*Vphs(2), w(1)*Vphs(1)+w(2)*Vphs(2));
  thetaGeo = atan2(w(2)*Vgeo(1)-w(1)*Vgeo(2), w(1)*Vgeo(1)+w(2)*Vgeo(2));
  
  Xphs = rotateVes(X,cPhys,theta);
  Xgeo = rotateVes(X,cGeom,theta);

  % for Xphs
  cPhys_Xphs = oc.getPhysicalCenterShan(Xphs);
  cGeom_Xphs = [mean(Xphs(1:end/2));mean(Xphs(end/2+1:end))];
  
  disp(['Physical Center of Xphs:' num2str(cPhys_Xphs')])
  disp(['Geometrical Center of Xphs:' num2str(cGeom_Xphs')])

  Vphs_Xphs = oc.getPrincAxesGivenCentroid(Xphs,cPhys_Xphs);
  disp(['PA of Xphs w.r.t. Physical Center:' num2str(Vphs_Xphs')])
  
  thetaPhs_Xphs = atan2(w(2)*Vphs_Xphs(1)-w(1)*Vphs_Xphs(2), w(1)*Vphs_Xphs(1)+w(2)*Vphs_Xphs(2));
  Xphs = [Xphs(1:end/2)-cPhys_Xphs(1); Xphs(end/2+1:end)-cPhys_Xphs(2)];
  
  % for Xgeo
  cPhys_Xgeo = oc.getPhysicalCenterShan(Xgeo);
  cGeom_Xgeo = [mean(Xgeo(1:end/2));mean(Xgeo(end/2+1:end))];
  
  disp(['Physical Center of Xgeo:' num2str(cPhys_Xgeo')])
  disp(['Geometrical Center of Xgeo:' num2str(cGeom_Xgeo')])

  Vphs_Xgeo = oc.getPrincAxesGivenCentroid(Xgeo,cPhys_Xphs);
  disp(['PA of Xgeo w.r.t. Physical Center:' num2str(Vphs_Xgeo')])

  thetaPhs_Xgeo = atan2(w(2)*Vphs_Xgeo(1)-w(1)*Vphs_Xgeo(2), w(1)*Vphs_Xgeo(1)+w(2)*Vphs_Xgeo(2));

  Xgeo = [Xgeo(1:end/2)-cPhys_Xgeo(1); Xgeo(end/2+1:end)-cPhys_Xgeo(2)];

  figure(1); clf;
  plot(X(1:end/2),X(end/2+1:end),'k','linewidth',2)
  hold on
  plot(Xphs(1:end/2),Xphs(end/2+1:end),'r','linewidth',2)
  plot(Xgeo(1:end/2),Xgeo(end/2+1:end),'b','linewidth',2)
  axis equal
  legend('X0','w.r.t physical','w.r.t. geometrical')

  
  pause

end

function Xrot = rotateVes(X,center,theta)

 Xrot = zeros(size(X));
 x = X(1:end/2); y = X(end/2+1:end);
 Xrot(1:end/2) = (x-center(1))*cos(theta) - (y-center(2))*sin(theta) + center(1);
 Xrot(end/2+1:end) = (x-center(1))*sin(theta) + (y-center(2))*cos(theta) + center(2);

end
  