clear all;clc

n = 64; 
X = boundary(n);

prams.T = 20;
prams.m = 200;
prams.kappa = 1e0;  
prams.order = 3;    
prams.Incompressibility = 1;
prams.vInf = @(X) farFieldVel(X,'shear',1);
                            
options.usePlot = 0;                   
options.findSteadyState = 1;                                       

vcRange = [.001 .01 1 2 4 6 8 10 100 1000];

figure;hAll = gca; figure;
lineColor = {'red','green','blue','cyan','magenta','yellow'}

for ii = 1:length(vcRange) 
  vc = vcRange(ii)
  clear functions global;clf
  prams.viscCont = vc;
  [Xfinal status] = Ves2D(X,prams,options);
  title(['Inclination angle vs. time (\lambda =' num2str(vc) ')']);
  xlabel('Non-dimensional time');
  ylabel('\theta (rad)');
  fileName = ['../results/incAngleTime_viscCont' num2str(vc) '.fig'];
  set(findobj(gca,'Type','line'),'Color',lineColor{mod(ii,6)+1}, 'LineWidth',2, ...
                    'DisplayName',['\lambda =' num2str(vc)]);
  legend show;
  saveas(gcf,fileName);
  hc = get(gca,'children');
  copyobj(hc,hAll);
  pause(1);
end

close;
title('Inclination angle vs. time');
xlabel('Non-dimensional time');
ylabel('\theta (rad)');
fileName = ['../results/incAngleTime_all.fig'];
legend show;
saveas(gcf,fileName);
