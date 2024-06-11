clear; clc;
load diverseParabolicDataSet.mat

N = numel(Xstore(:,1))/2;
nInstances = numel(Xstore(1,:));

XstoreMirrY = zeros(2*N,nInstances);
XstoreMirrYX = zeros(2*N,nInstances);
XstoreMirrX = zeros(2*N,nInstances);

for k = 1 : nInstances
  
  XstoreMirrY(1:end/2,k) = -Xstore(1:end/2,k);
  XstoreMirrY(end/2+1:end,k) = Xstore(end/2+1:end,k);

  XstoreMirrYX(1:end/2,k) = XstoreMirrY(1:end/2,k);
  XstoreMirrYX(end/2+1:end,k) = -XstoreMirrY(end/2+1:end,k);

  XstoreMirrX(1:end/2,k) = -XstoreMirrYX(1:end/2,k);
  XstoreMirrX(end/2+1:end,k) = XstoreMirrYX(end/2+1:end,k);

  % figure(1); clf;
  % plot(Xstore(1:end/2,k),Xstore(end/2+1:end,k),'k')
  % hold on
  % plot(XstoreMirrY(1:end/2,k),XstoreMirrY(end/2+1:end,k),'r')
  % plot(XstoreMirrYX(1:end/2,k),XstoreMirrYX(end/2+1:end,k),'b')
  % plot(XstoreMirrX(1:end/2,k),XstoreMirrX(end/2+1:end,k),'g')
  % axis equal
  % title(k)
  % pause
  
end

Xstore = [Xstore XstoreMirrY XstoreMirrYX XstoreMirrX];

save MirrorredParabolicDataSet Xstore


% save n128Dt1e-06RelaxMirrdDataSet.mat XstandStore XnewStandStore nInstances