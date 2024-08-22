set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

fileNameTR = 'relax_N128_EXACT_dt1e-05_speed2000.bin';
fileNameNN = 'more_smooth_relax_N128_Network_dt1e-05_speed2000.bin';


[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameNN);

tsteps = [1; 2; 5; 10; 100; 200];
cx = [0; 0.25; 0.5; 0.75; 1.0; 1.25];
if 1
figure(1); clf;
for ik = 1 : numel(tsteps)
  k = tsteps(ik);  
  xvecT = [vesxT(:,1,k);vesxT(1,1,k)] + cx(ik);
  yvecT = [vesyT(:,1,k);vesyT(1,1,k)];
  
  xvecN = [vesxN(:,1,k);vesxN(1,1,k)]+ cx(ik);
  yvecN = [vesyN(:,1,k);vesyN(1,1,k)];

  h = plot(xvecT, yvecT, 'Color',[0 0 0 1],'linewidth',2);
  set(h,'Color',[h.Color, 1],'linewidth',2)
  hold on

  h2 = plot(xvecN, yvecN, 'Color',[1 0 0 0.65],'linewidth',2);
  set(h2,'Color',[h2.Color, 0.65],'linewidth',2)
  

  

end
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);
    
set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')

axis equal

figName = ['All.png'];
ax = gca;
exportgraphics(ax,['~/Desktop/relax' figName],'Resolution',300)

end

if 1
oc = curve;
bendingT = zeros(numel(timeT),1);
bendingN = zeros(numel(timeN),1);
for k = 1 : numel(timeT)
[jac, tan, curv] = oc.diffProp([vesxT(:,k);vesyT(:,k)]);
bendingT(k) = 1/2 * 1/128 * sum(curv.^2 .* jac);

[jac, tan, curv] = oc.diffProp([vesxN(:,k);vesyN(:,k)]);
bendingN(k) = 1/2 * 1/128 * sum(curv.^2 .* jac);
end

figure(1);clf;
loglog(1:numel(timeT),bendingT,'k','linewidth',2)
hold on
loglog(1:numel(timeN),bendingN,'r','linewidth',2)
axis square
box on
ylim([1 200])

 ax = gca;
exportgraphics(ax,['~/Desktop/relaxBendingEnergy.png'],'Resolution',300)
end

