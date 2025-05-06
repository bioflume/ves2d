clear; clc;

load selfTensionDataSet.mat

for iter = 1 : 21

  load(['~/Desktop/networkBIEMsolveCompareDataForSet/comparisonX' num2str(iter) '.mat'])

  figure(1); clf;
  subplot(2,2,1)
  plot(linspace(0,1,32)',trueSelfBendSolve,'linewidth',2)
  hold on
  plot(linspace(0,1,32)',interpft(trueSelfBendSolve128,32),'linewidth',2)
  plot(linspace(0,1,32)',interpft(Toutput(:,iter),32),'--','linewidth',2)
  plot(linspace(0,1,32)',selfBendSolve,'linewidth',2)
  axis square
  legend('True solve (N = 32)','True solve (N = 128) - Down to N = 32','Dataset Down to N = 32','Network solve (N = 32)','location','north')
  grid
  box on
  xlabel('s')
  ylabel('Tension')

  z = Xinput(1:end/2,iter) + 1i*Xinput(end/2+1:end,iter);
  zh = fft(z);

  subplot(2,2,2)
  plot(Xinput(1:end/2,iter),Xinput(end/2+1:end,iter),'linewidth',2)
  axis equal

  subplot(2,2,3)
  plot(real(zh))
  axis square
  title('Real components of fft(X)')
  ylim([-20 20])

  subplot(2,2,4)
  plot(imag(zh))
  axis square
  title('Imaginary components of fft(X)')
  ylim([-2 2])

  % figure(1);
  % ax = gca;
  % exportgraphics(ax,['~/Desktop/exampleShape' num2str(iter) '.png'],'Resolution',300)
  pause




end

%%
% 
% load selfTensionDataSet.mat
% zImag = zeros(32,32,21);
% zReal = zeros(32,32,21);
% for imode = 1 : 32
%   load(['~/Desktop/tensionAdvectNet_modes32/TensionAdvectNetFFTBasisN32_mode' num2str(imode) 'Data.mat'])
%   for iter = 1 : 21
%     zImag(:,imode,iter) = zImagStore(:,iter);
%     zReal(:,imode,iter) = zRealStore(:,iter);
%   end
% end
% 
% vinf = @(X) [X(end/2+1:end); zeros(size(X(1:end/2)))];
% 
% for iter = 1 : 21
%   v = vinf([interpft(Xinput(1:end/2,iter),32);interpft(Xinput(end/2+1:end,iter),32)]);
%   z = v(1:end/2)+1i*v(end/2+1:end);
%   zh = fft(z);
%   V1 = real(zh); V2 = imag(zh);
%   tension = zReal(:,:,iter)*V1 + zImag(:,:,iter)*V2;
% 
%   figure(1); clf;
%   subplot(1,2,1)
%   plot(linspace(0,1,32)',tension,'linewidth',2)
%   axis square
%   grid
%   box on
%   xlabel('s')
%   ylabel('Tension')
%   subplot(1,2,2)
%   plot(Xinput(1:end/2,iter),Xinput(end/2+1:end,iter),'linewidth',2)
%   axis equal
% 
%   pause
% 
% end