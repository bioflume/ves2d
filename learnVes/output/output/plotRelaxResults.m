clear;clc;

imovie  = ~false;
isaveImages = false;
skip = 1;
imeasureError = ~false;
if isaveImages
  count = 1;
  mkdir frames
end

% The WORST CASES
% 9309 -- 8%
% 7001 -- 6%
% 5816 -- 6%
% 47409 -- 5%

% AVERAGE CASES (~1%)
% 12722; 63348; 9772; 91724 

vesID = [];
% curly, star, ellipse, openStar
initShape = 'openStar';
% load(['coarseTrueSimVesShapecurly' '_bgFlowrelax_speed0.mat'])
% XhistCoarse = XhistTrue; timeCoarse = timeTrue;

if ~isempty(vesID)
load(['N256xLDNNsimVesID' num2str(vesID) '_bgFlowrelax_speed0.mat'])
else
load(['N256xLDNNsimVesShape' initShape '_bgFlowrelax_speed0.mat'])    
end

% XhistTrue = XhistTrue(:,1:101);
% Xhist = Xhist(:,1:101); time = time(1:101);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if imovie
for k = 1 : skip : numel(time)
    
figure(1); clf;
plot([XhistTrue(1:end/2,k);XhistTrue(1,k)],[XhistTrue(end/2+1:end,k);...
    XhistTrue(end/2+1,k)],'Color',[.5 .5 .5],'linewidth',2)
hold on
plot([Xhist(1:end/2,k);Xhist(1,k)],[Xhist(end/2+1:end,k);...
    Xhist(end/2+1,k)],'r','linewidth',2)
axis equal
xlim([-0.3 0.3])
ylim([-0.3 0.3])
title(time(k));
legend('True','DNN')
legend boxoff
box on

if isaveImages
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'ztick',[]);

  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  set(gca,'zcolor','w');
  box off
  set(gca,'visible','off')
  titleStr = ['t = ' num2str(time(k),'%.3f')];
  text(-0.1,0.3,titleStr,'FontSize',28,'FontName','Palatino') 
end

if isaveImages
  filename = ['./frames/image', sprintf('%04d',count),'.png'];
  count = count+1;
  figure(1);
  print(gcf,'-dpng','-r300',filename); 
else
  pause(0.1);   
end       

end % for k = 1 : skip : numel(time)  
end % if imovie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if imeasureError
addpath ~/padas/Ves2Dn/src/
oc = curve;

trueIAstore = zeros(numel(time),1);
trueCentStore = zeros(numel(time),2);
trueBendEgyStore = zeros(numel(time),1);
predIAstore = zeros(numel(time),1);
predCentStore = zeros(numel(time),2);
predBendEgyStore = zeros(numel(time),1);

errInShapeStore = zeros(numel(time),1);
Ntrue = numel(XhistTrue(1:end/2,1)); Npred = numel(Xhist(1:end/2,1));
for k = 1 : numel(time)
  Xtrue = XhistTrue(:,k); Xpred = Xhist(:,k);
  
  trueIAstore(k) = oc.getIncAngle2(Xtrue);
  trueCentStore(k,:) = mean([Xtrue(1:end/2) Xtrue(end/2+1:end)]);
  [jac,tan,curv] = oc.diffProp(Xtrue);
  trueBendEgyStore(k) = sum(0.5*jac.*curv.^2)*2*pi/Ntrue;
  
  predIAstore(k) = oc.getIncAngle2(Xpred);
  predCentStore(k,:) = mean([Xpred(1:end/2) Xpred(end/2+1:end)]);
  [jac,tan,curv] = oc.diffProp(Xpred);
  predBendEgyStore(k) = sum(0.5*jac.*curv.^2)*2*pi/Npred;
  
  if Ntrue ~= Npred
    Xpred = [interpft(Xpred(1:end/2),Ntrue);interpft(Xpred(end/2+1:end),Ntrue)];    
  end
  
  errInShapeStore(k) = sqrt(1/Ntrue*sum((Xtrue(1:end/2)-Xpred(1:end/2)).^2+...
      (Xtrue(end/2+1:end)-Xpred(end/2+1:end)).^2))./sqrt(1/Ntrue*sum(Xtrue.^2));
    
end    
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% Error in IA
errIA = abs(trueIAstore(end)-predIAstore(end));
if abs(trueIAstore(end))<1e-10
  errIA = errIA/1e-10;
else
  errIA = errIA/abs(trueIAstore(end));
end
disp(['Error in IA in equil. : ' num2str(errIA)])
% Error in Center
errCent = sqrt((trueCentStore(end,1)-predCentStore(end,1))^2+...
    (trueCentStore(end,2)-predCentStore(end,2))^2)/sqrt(sum(trueCentStore(end,:).^2));
disp(['Error in center in equil. : ' num2str(errCent)])

% Error in shape
figure(2);clf;
plot(time,errInShapeStore,'b','linewidth',2)
axis square
xlabel('Time')
ylabel('Error in shape')
title(['VesID: ' num2str(vesID)])

% Error in Bending Energy
figure(3); clf;
plot(time,trueBendEgyStore,'Color',[.5 .5 .5],'linewidth',2)
hold on
plot(time,predBendEgyStore,'r','linewidth',2)
axis square
xlabel('Time')
ylabel('Bending energy')
title(['VesID: ' num2str(vesID)])

errInBendEgy = abs(trueBendEgyStore(end)-predBendEgyStore(end))/abs(trueBendEgyStore(end));
disp(['Error in bending energy in equil. : ' num2str(errInBendEgy)])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
% [jac,tan,curv] = oc.diffProp(XhistCoarse(:,end));
% coarseBendEgy = sum(0.5*jac.*curv.^2)*2*pi/96;    
% errInCoarseBendEgy = abs(trueBendEgyStore(end)-coarseBendEgy)/abs(coarseBendEgy);
% disp(['Error in coarse bending energy in equil. : ' num2str(errInCoarseBendEgy)])
% coarseIA = oc.getIncAngle(XhistCoarse(:,end));
% errCoarseIA = abs(trueIAstore(end)-coarseIA);
% if abs(trueIAstore(end))<1e-10
%   errCoarseIA = errCoarseIA/1e-10;
% else
%   errCoarseIA = errCoarseIA/abs(trueIAstore(end));
% end
% disp(['Error in Coarse IA in equil. : ' num2str(errCoarseIA)])
end % if imeasureError
