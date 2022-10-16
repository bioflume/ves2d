clear; clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

fileGT = 'nv70N64VF30GT_bgFlowcouette_speed100';
load(fileGT)
XTrue = XhistTrue;
tTrue = timeTrue;

fileDNN = './repulsion/repStr2000_repScale1.5VF30.bin';
[vesx, vesy, ten, etax, etay, time, N, nv, xinit, yinit, xwalls, ywalls, ncountNN, ncountExact] = loadManyVesFile(fileDNN);
X = [vesx;vesy];
tDNN = time;

ntime = numel(tDNN); nv = numel(X(1,:,1));
rdist = zeros(nv,ntime);
rdistTrue = rdist;
for it = 1 : 1: ntime
  idx = find(abs(tTrue-70/100*tDNN(it)) <= 1e-3);
  if ~isempty(idx)
    itT = idx(1);
  else
    disp('Could not find corresponding time')
    pause
  end
  for k = 1 : nv
    cx = mean(X(1:end/2,k,it)); cy = mean(X(end/2+1:end,k,it));    
    rdist(k,it) = sqrt(cx.^2 + cy.^2);
    
    cx = mean(XTrue(1:end/2,k,itT)); cy = mean(XTrue(end/2+1:end,k,itT));    
    rdistTrue(k,it) = sqrt(cx.^2 + cy.^2);
  end
end
%%
pts = linspace(1,2.2,100)';
count = 1;
for it = 1 : 50: ntime
  [f, xi] = ksdensity(rdist(:,it), pts);
  
  f = f./trapz(xi,f);
  
  [fT, xiT] = ksdensity(rdistTrue(:,it), pts);
  
  fT = fT./trapz(xiT,fT);
  
  figure(1); clf; hold on;
  plot(xi, f, 'r','linewidth',2)
  hold on
  plot(xiT, fT, 'k','linewidth',2)
  title(['Time = ' num2str(tDNN(it))])
  legend('DNN','True')
  legend boxoff
  axis square
  grid
  xlabel('Radial position')
  box on
  ylim([0 1.5])
  xlim([1 2.2])
  filename = ['./frames/image', sprintf('%04d',count),'.png'];
  figure(1);
  print(gcf,'-dpng','-r300',filename);
  count = count + 1;
end
  
