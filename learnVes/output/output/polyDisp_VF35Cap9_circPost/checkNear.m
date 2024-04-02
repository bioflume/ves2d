% clear;
% addpath ../../src/
% oc = curve;
% N = 32;
% nv = 1;
% Nbd = 512;
% nvbd = 2;
% kappa = 2e-1;
% 
% % load alternTstep_MLARMcouette_nves1_Nves32_Kappa0.2_Dt0.0001_Speed100.mat
% load nv1N32DNNlikeTrueNearDNNlikeExactTension_Kb2E1Dt1E4_1ves_bgFlowcouette_speed100
% XDNN = Xhist(:,:,1:it);
% sigDNN = sigStore(:,:,1:it);
% etaDNN = etaStore(:,:,1:it);
% RSDNN = RSstore(:,:,1:it);
% it1 = it;
% 
% % load nv1N32DNNTrueNear_Kb2E1Dt1E4_1ves_bgFlowcouette_speed100
% load nv1N32DNNlikeTrueNearDNNlikeExact_Kb2E1Dt1E4_1ves_bgFlowcouette_speed100
% XtrueNear = Xhist(:,:,1:it);
% sigTrueNear = sigStore(:,:,1:it);
% etaTrueNear = etaStore(:,:,1:it);
% RStrueNear = RSstore(:,:,1:it);
% it2 = it;
% 
% load TrueN32_1vesKb2e1Dt1e4_bgFlowcouette_speed100
% Xtrue = XhistTrue(:,:,:);
% sigTrue = sigStore(:,:,:);
% etaTrue = etaStore(:,:,:);
% RStrue = RSstore(:,:,:);
% 
% it = min([it1;it2]);
% it = it2;

% load singleVesVelocity.mat

load singleVesVelocityDNNlikes.mat
if 1
  centRadDNN = zeros(it,1); centRadTrueNear = centRadDNN;
  centRadTrue = centRadDNN;
  for k = 1 : it
    
  
  cx = mean(XDNN(1:end/2,k)); cy = mean(XDNN(end/2+1:end,k));
  centRadDNN(k) = sqrt(cx^2+cy^2);
  
  cx = mean(XtrueNear(1:end/2,k)); cy = mean(XtrueNear(end/2+1:end,k));
  centRadTrueNear(k) = sqrt(cx^2+cy^2);
  
  cx = mean(Xtrue(1:end/2,k)); cy = mean(Xtrue(end/2+1:end,k));
  centRadTrue(k) = sqrt(cx^2+cy^2);
  end
end
    
if 0
centRadDNN = zeros(it,1); centRadTrueNear = centRadDNN;
centRadTrue = centRadDNN;
% BUILD WALLS
walls = capsules(Xwalls,[],[],zeros(nvbd,1),zeros(nvbd,1),1);
walls.setUpRate();
opWall = poten(Nbd);
NearW2W = walls.getZone([],1);
wallDLP = opWall.stokesDLmatrix(walls);
wallDLPnoCorr = opWall.stokesDLmatrixNoCorr(walls);
wallN0 = opWall.stokesN0matrix(walls);

Fwall2VesDNN = zeros(2*N,it);
Fwall2VesTN = Fwall2VesDNN; Fwall2VesTR = Fwall2VesDNN;

for k = 2 : it
    
  if 0
  cx = mean(XDNN(1:end/2,k)); cy = mean(XDNN(end/2+1:end,k));
  centRadDNN(k) = sqrt(cx^2+cy^2);
  
  cx = mean(XtrueNear(1:end/2,k)); cy = mean(XtrueNear(end/2+1:end,k));
  centRadTrueNear(k) = sqrt(cx^2+cy^2);
  
  cx = mean(Xtrue(1:end/2,k)); cy = mean(Xtrue(end/2+1:end,k));
  centRadTrue(k) = sqrt(cx^2+cy^2);
  end
  
  % BUILD VESICLES
  vesicleDNN = capsules(XDNN(:,:,k),[],[],kappa,1,1);
  vesicleDNN.setUpRate();
  
  vesicleTN = capsules(XtrueNear(:,:,k),[],[],kappa,1,1);
  vesicleTN.setUpRate();
  
  vesicleTR = capsules(Xtrue(:,:,k),[],[],kappa,1,1);
  vesicleTR.setUpRate();
  
  % Get near structure
  [NearV2Vdnn,NearV2Wdnn] = vesicleDNN.getZone(walls,3);
  [~,NearW2Vdnn] = walls.getZone(vesicleDNN,2);
  
  [NearV2Vtn,NearV2Wtn] = vesicleTN.getZone(walls,3);
  [~,NearW2Vtn] = walls.getZone(vesicleTN,2);
  
  [NearV2Vtr,NearV2Wtr] = vesicleTR.getZone(walls,3);
  [~,NearW2Vtr] = walls.getZone(vesicleTR,2);
  
  % Compute wall-2-vesicle interactions
  kernel = @opWall.exactStokesDL;
  kernelDirect = @opWall.exactStokesDL;
  
  DLP = @(X) -1/2*X + opWall.exactStokesDLdiag(walls,wallDLP,X);
  Fwall2VesDNN(:,k) = opWall.nearSingInt(walls,etaTrue(:,:,k),DLP,[],...
      NearW2Vdnn,kernel,kernelDirect,vesicleDNN,false,false);
  
  Fwall2VesTN(:,k) = opWall.nearSingInt(walls,etaTrue(:,:,k),DLP,[],...
      NearW2Vtn,kernel,kernelDirect,vesicleTN,false,false);
  
  Fwall2VesTR(:,k) = opWall.nearSingInt(walls,etaTrue(:,:,k),DLP,[],...
      NearW2Vtr,kernel,kernelDirect,vesicleTR,false,false);
  
  stokeslet = RStrue(1:2,2,k);
  rotlet = RStrue(3,2,k);
  % the center of the rotlet/stokeslet terms
  [cx,cy] = oc.getXY(walls.center(:,2));
  
  % set of points where we are evaluating the velocity
  [x,y] = oc.getXY(XDNN(:,k));
  
  % distance squared
  rho2 = (x-cx).^2 + (y-cy).^2;

  % x component of velocity due to the stokeslet and rotlet
  LogTerm = -0.5*log(rho2)*stokeslet(1);
  rorTerm = 1./rho2.*((x-cx).*(x-cx)*stokeslet(1) + ...
   (x-cx).*(y-cy)*stokeslet(2));
  RotTerm = (y-cy)./rho2*rotlet;
  velx = 1/4/pi*(LogTerm + rorTerm) + RotTerm;

  % y component of velocity due to the stokeslet and rotlet
  LogTerm = -0.5*log(rho2)*stokeslet(2);
  rorTerm = 1./rho2.*((y-cy).*(x-cx)*stokeslet(1) + ...
   (y-cy).*(y-cy)*stokeslet(2));
  RotTerm = -(x-cx)./rho2*rotlet;
  vely = 1/4/pi*(LogTerm + rorTerm) + RotTerm;

   % velocity
  wallLets_DNN(:,k) = [velx;vely];
   
   
   % set of points where we are evaluating the velocity
  [x,y] = oc.getXY(XtrueNear(:,k));
  
  % distance squared
  rho2 = (x-cx).^2 + (y-cy).^2;

  % x component of velocity due to the stokeslet and rotlet
  LogTerm = -0.5*log(rho2)*stokeslet(1);
  rorTerm = 1./rho2.*((x-cx).*(x-cx)*stokeslet(1) + ...
   (x-cx).*(y-cy)*stokeslet(2));
  RotTerm = (y-cy)./rho2*rotlet;
  velx = 1/4/pi*(LogTerm + rorTerm) + RotTerm;

  % y component of velocity due to the stokeslet and rotlet
  LogTerm = -0.5*log(rho2)*stokeslet(2);
  rorTerm = 1./rho2.*((y-cy).*(x-cx)*stokeslet(1) + ...
   (y-cy).*(y-cy)*stokeslet(2));
  RotTerm = -(x-cx)./rho2*rotlet;
  vely = 1/4/pi*(LogTerm + rorTerm) + RotTerm;

  % velocity
  wallLets_TN(:,k) = [velx;vely];
      
  % set of points where we are evaluating the velocity
  [x,y] = oc.getXY(Xtrue(:,k));
  
  % distance squared
  rho2 = (x-cx).^2 + (y-cy).^2;

  % x component of velocity due to the stokeslet and rotlet
  LogTerm = -0.5*log(rho2)*stokeslet(1);
  rorTerm = 1./rho2.*((x-cx).*(x-cx)*stokeslet(1) + ...
   (x-cx).*(y-cy)*stokeslet(2));
  RotTerm = (y-cy)./rho2*rotlet;
  velx = 1/4/pi*(LogTerm + rorTerm) + RotTerm;

  % y component of velocity due to the stokeslet and rotlet
  LogTerm = -0.5*log(rho2)*stokeslet(2);
  rorTerm = 1./rho2.*((y-cy).*(x-cx)*stokeslet(1) + ...
   (y-cy).*(y-cy)*stokeslet(2));
  RotTerm = -(x-cx)./rho2*rotlet;
  vely = 1/4/pi*(LogTerm + rorTerm) + RotTerm;

   % velocity
   wallLets_TR(:,k) = [velx;vely];   
   
   k
end

Fwall2VesDNN = wallLets_DNN + Fwall2VesDNN;
Fwall2VesTN = wallLets_TN + Fwall2VesTN;
Fwall2VesTR = wallLets_TR + Fwall2VesTR;
end

if 1
velR_DNN = zeros(it,1); velR_TN = velR_DNN; velR_TR = velR_DNN;
velT_DNN = zeros(it,1); velT_TN = velT_DNN; velT_TR = velT_DNN;
velR_DNN_x = zeros(it,1); velR_DNN_y = zeros(it,1);
velR_TN_x = zeros(it,1); velR_TN_y = zeros(it,1);
velR_TR_x = zeros(it,1); velR_TR_y = zeros(it,1);
for k = 2 : it
  velx = mean(Fwall2VesDNN(1:end/2,k)); vely = mean(Fwall2VesDNN(end/2+1:end,k));
  alpha = atan2(mean(XDNN(end/2+1:end,k)),mean(XDNN(1:end/2,k)));
  velR_DNN(k) = velx*cos(alpha)+vely*sin(alpha);
  velT_DNN(k) = vely*cos(alpha)-velx*sin(alpha);
  velR_DNN_x(k) = velR_DNN(k)*cos(alpha);
  velR_DNN_y(k) = velR_DNN(k)*sin(alpha);
  
  velx = mean(Fwall2VesTN(1:end/2,k)); vely = mean(Fwall2VesTN(end/2+1:end,k));
  alpha = atan2(mean(XtrueNear(end/2+1:end,k)),mean(XtrueNear(1:end/2,k)));
  velR_TN(k) = velx*cos(alpha)+vely*sin(alpha);
  velT_TN(k) = vely*cos(alpha)-velx*sin(alpha);
  velR_TN_x(k) = velR_TN(k)*cos(alpha);
  velR_TN_y(k) = velR_TN(k)*sin(alpha);
  
  velx = mean(Fwall2VesTR(1:end/2,k)); vely = mean(Fwall2VesTR(end/2+1:end,k));
  alpha = atan2(mean(Xtrue(end/2+1:end,k)),mean(Xtrue(1:end/2,k)));
  velR_TR(k) = velx*cos(alpha)+vely*sin(alpha);
  velT_TR(k) = vely*cos(alpha)-velx*sin(alpha);
  velR_TR_x(k) = velR_TR(k)*cos(alpha);
  velR_TR_y(k) = velR_TR(k)*sin(alpha);
    
end
end

% save singleVesVelocityDNNlikes Fwall2VesDNN Fwall2VesTN Fwall2VesTR XDNN Xtrue XtrueNear Xwalls it time

if 1
titY = 2.4; titX = -0.6;
xlimits = [min(min(Xwalls(1:end/2,:))) max(max(Xwalls(1:end/2,:)))];
ylimits = [min(min(Xwalls(1+end/2:end,:))) max(max(Xwalls(1+end/2:end,:)))];
xwalls = [Xwalls(1:end/2,:);Xwalls(1,:)]; 
ywalls = [Xwalls(end/2+1:end,:);Xwalls(end/2+1,:)];
count = 1;

for k = 1 : 1 : it
  figure(1);clf
  
%   subplot(1,2,1);
  plot(xwalls,ywalls,'Color',[.5 .5 .5],'linewidth',2)
  hold on;
  plot([interpft(XDNN(1:end/2,k),128);XDNN(1,k)],[interpft(XDNN(end/2+1:end,k),128);XDNN(end/2+1,k)],'r','linewidth',1.5)
%   quiver(XDNN(1:end/2,k),XDNN(end/2+1:end,k),Fwall2VesDNN(1:end/2,k),Fwall2VesDNN(end/2+1:end,k),'r')
%   mvelx = mean(Fwall2VesDNN(1:end/2,k)); mvely = mean(Fwall2VesDNN(end/2+1:end,k));
  quiver(mean(XDNN(1:end/2,k)),mean(XDNN(end/2+1:end,k)),velR_DNN_x(k),velR_DNN_y(k),4e-1,'r','linewidth',1.5)
%   quiver(mean(XDNN(1:end/2,k)),mean(XDNN(end/2+1:end,k)),mean(XDNN(1:end/2,k)),mean(XDNN(end/2+1:end,k)),3e-1,'Color',[0 .5 0])

  plot([interpft(XtrueNear(1:end/2,k),128);XtrueNear(1,k)],[interpft(XtrueNear(end/2+1:end,k),128);XtrueNear(end/2+1,k)],'b','linewidth',1.5)
%   quiver(XtrueNear(1:end/2,k),XtrueNear(end/2+1:end,k),Fwall2VesTN(1:end/2,k),Fwall2VesTN(end/2+1:end,k),'b')
%   mvelx = mean(Fwall2VesTN(1:end/2,k)); mvely = mean(Fwall2VesTN(end/2+1:end,k));
  quiver(mean(XtrueNear(1:end/2,k)),mean(XtrueNear(end/2+1:end,k)),velR_TN_x(k),velR_TN_y(k),4e-1,'b','linewidth',1.5)
%   quiver(mean(XtrueNear(1:end/2,k)),mean(XtrueNear(end/2+1:end,k)),mean(XtrueNear(1:end/2,k)),mean(XtrueNear(end/2+1:end,k)),3e-1,'Color',[0 .5 0])
  
  plot([interpft(Xtrue(1:end/2,k),128);Xtrue(1,k)],[interpft(Xtrue(end/2+1:end,k),128);Xtrue(end/2+1,k)],'k','linewidth',1.5)
%   quiver(Xtrue(1:end/2,k),Xtrue(end/2+1:end,k),Fwall2VesTR(1:end/2,k),Fwall2VesTR(end/2+1:end,k),'k')
%   mvelx = mean(Fwall2VesTR(1:end/2,k)); mvely = mean(Fwall2VesTR(end/2+1:end,k));
  quiver(mean(Xtrue(1:end/2,k)),mean(Xtrue(end/2+1:end,k)),velR_TR_x(k),velR_TR_y(k),4e-1,'k','linewidth',1.5)
%   quiver(mean(Xtrue(1:end/2,k)),mean(Xtrue(end/2+1:end,k)),mean(Xtrue(1:end/2,k)),mean(Xtrue(end/2+1:end,k)),3e-1,'Color',[0 .5 0])

  axis equal
  
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'ztick',[]);

  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  set(gca,'zcolor','w');
  box off
  set(gca,'visible','off')
  titleStr = ['t = ' num2str(time(k),'%.2f')];
  text(titX,titY,titleStr,'FontSize',28,'FontName','Palatino') 
  
  filename = ['./frames/image', sprintf('%04d',count),'.png'];
  count = count+1;
  figure(1);
  pause;
%   print(gcf,'-dpng','-r300',filename);
  
  
  
%   subplot(1,2,2); hold on;
%   plot(time(1:it),velR_TR,'k','linewidth',2)
%   plot(time(1:it),velR_TN,'b','linewidth',2)
%   plot(time(1:it),velR_DNN,'r','linewidth',2)
%   plot(time(k),velR_TR(k),'ko','markersize',8,'markerfacecolor','k')
%   plot(time(k),velR_TN(k),'bo','markersize',8,'markerfacecolor','b')
%   plot(time(k),velR_DNN(k),'ro','markersize',8,'markerfacecolor','r')
%   axis square
%   box on
%   xlabel('Time')
%   ylabel('Migration velocity')
%   legend('Truth','MLARM + true near','MLARM')
%   legend boxoff
%   xlim([time(1) time(it)])
%   ylim([-4 4.5])
%   grid on
%   pause(0.1)
end

end



if 0
  alphaDNN = zeros(it,1); alphaTN = zeros(it,1); alphaTR = zeros(it,1);
  IADNN = zeros(it,1); IATN = zeros(it,1); IATR = zeros(it,1);
  for k = 1 : it
      
    cx = mean(XDNN(end/2+1:end,k)); cy = mean(XDNN(1:end/2,k));
    alphaDNN(k) = atan2(abs(cy),abs(cx));        
    if cx < 0 && cy > 0
      alphaDNN(k) = pi - alphaDNN(k);
    elseif cx < 0 && cy < 0
      alphaDNN(k) = pi + alphaDNN(k);  
    elseif cx > 0 && cy < 0
      alphaDNN(k) = 2*pi - alphaDNN(k);  
    end
      
    IADNN(k) = oc.getIncAngle(XDNN(:,k));
    if k >= 2
    if abs(IADNN(k)-IADNN(k-1)) > 0.5*pi 
      IADNN(k) = IADNN(k) + pi;
    end
    end
    
    cx = mean(XtrueNear(end/2+1:end,k)); cy = mean(XtrueNear(1:end/2,k));
    alphaTN(k) = atan2(abs(cy),abs(cx));              
    if cx < 0 && cy > 0
      alphaTN(k) = pi - alphaTN(k);
    elseif cx < 0 && cy < 0
      alphaTN(k) = pi + alphaTN(k);  
    elseif cx > 0 && cy < 0
      alphaTN(k) = 2*pi - alphaTN(k);  
    end
    
    IATN(k) = oc.getIncAngle(XtrueNear(:,k));
    if k >= 2
    if abs(IATN(k)-IATN(k-1)) > 0.5*pi 
      IATN(k) = IATN(k) + pi;
    end    
    end
    
    cx = mean(Xtrue(end/2+1:end,k)); cy = mean(Xtrue(1:end/2,k));
    alphaTR(k) = atan2(abs(cy),abs(cx));             
    if cx < 0 && cy > 0
      alphaTR(k) = pi - alphaTR(k);
    elseif cx < 0 && cy < 0
      alphaTR(k) = pi + alphaTR(k);  
    elseif cx > 0 && cy < 0
      alphaTR(k) = 2*pi - alphaTR(k);  
    end
    
    IATR(k) = oc.getIncAngle(Xtrue(:,k));
    if k >= 2
    if abs(IATR(k)-IATR(k-1)) > 0.5*pi 
      IATR(k) = IATR(k) + pi;
    end
    end
    
  end

  alphaTanDNN = rem(alphaDNN + pi/2,2*pi);
  IADNN_Tan = alphaTanDNN-IADNN;  
  
  alphaTanTN = rem(alphaTN + pi/2,2*pi);
  IATN_Tan = alphaTanTN-IATN;  
  
  alphaTanTR = rem(alphaTR + pi/2,2*pi);  
  IATR_Tan = alphaTanTR-IATR;  
    
end