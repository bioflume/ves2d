clear; clc;

R = 0.1291;
chanWidths = R./[0.2; 0.4; 0.6; 0.75];
speeds = [6000 7500;
          3000 3750;
          2000 2500;
          1600 2000];


count = 1;
for iw = 1 : numel(chanWidths)
    for is = 1 : 2
      w = chanWidths(iw);
      vmax = speeds(iw,is);
      Cks(count,1) = 2*vmax*R^3/w;
      Cns(count,1) = R/w;
      CnVal(iw, is) = Cns(count,1);
      CkVal(iw, is) = Cks(count,1);
      count = count + 1;
    end
end


%% Load final vesicle configurations


% True ones
runName = 'resume_poisTrueRuns_dt5e-06_speed6000_width0.6455.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfTrue{1,1} = [vesxN(:,end);vesyN(:,end)];

runName = 'resume_poisTrueRuns_dt5e-06_speed7500_width0.6455.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfTrue{1,2} = [vesxN(:,end);vesyN(:,end)];

runName = 'resume_poisTrueRuns_dt5e-06_speed3000_width0.32275.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfTrue{2,1} = [vesxN(:,end);vesyN(:,end)];

runName = 'resume_poisTrueRuns_dt5e-06_speed3750_width0.32275.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfTrue{2,2} = [vesxN(:,end);vesyN(:,end)];

runName = 'poisTrueRuns_dt2.5e-06_speed2000_width0.21517.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfTrue{3,1} = [vesxN(:,end);vesyN(:,end)];

runName = 'poisTrueRuns_dt2.5e-06_speed2500_width0.21517.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfTrue{3,2} = [vesxN(:,end);vesyN(:,end)];

runName = 'poisTrueRuns_dt2.5e-06_speed1600_width0.17213.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfTrue{4,1} = [vesxN(:,end);vesyN(:,end)];

runName = 'poisTrueRuns_dt2.5e-06_speed2000_width0.17213.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfTrue{4,2} = [vesxN(:,end);vesyN(:,end)];

% DNN ones
runName = 'test_exAdv_diff625kNetJune8_dt5e-06poisRuns_speed6000_width0.6455.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfDNN{1,1} = [vesxN(:,end-1200);vesyN(:,end-1200)];

runName = 'resume_test_exAdv_diff625kNetJune8_dt5e-06poisRuns_speed7500_width0.6455.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfDNN{1,2} = [vesxN(:,end-200);vesyN(:,end-200)];

runName = 'resume_test_exAdv_diff625kNetJune8_dt5e-06poisRuns_speed3000_width0.32275.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfDNN{2,1} = [vesxN(:,end-400);vesyN(:,end-400)];

runName = 'resume_test_exAdv_diff625kNetJune8_dt5e-06poisRuns_speed3750_width0.32275.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfDNN{2,2} = [vesxN(:,end-10);vesyN(:,end-10)];

runName = 'test_exAdv_diff625kNetJune8_dt2.5e-06poisRuns_speed2000_width0.21517.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfDNN{3,1} = [vesxN(:,end);vesyN(:,end)];

runName = 'test_exAdv_diff625kNetJune8_dt2.5e-06poisRuns_speed2500_width0.21517.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfDNN{3,2} = [vesxN(:,end);vesyN(:,end)];

runName = 'test_exAdv_diff625kNetJune8_dt2.5e-06poisRuns_speed1600_width0.17213.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfDNN{4,1} = [vesxN(:,end-60);vesyN(:,end-60)];

runName = 'test_exAdv_diff625kNetJune8_dt2.5e-06poisRuns_speed2000_width0.17213.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runName);
XfDNN{4,2} = [vesxN(:,end);vesyN(:,end)];


%% 
figure(1); clf; hold on; % True
figure(2); clf; hold on; % DNN

for iw = 1 : 4
  for is = 1 : 2
    cxNew = (iw-1)*0.5;
    cyNew = (is-1)*0.5;
    
    
    figure(1);
    xvec = XfTrue{iw,is}(1:end/2);
    yvec = XfTrue{iw,is}(end/2+1:end);

    cxTrue(iw,is) = mean(xvec);
    cyTrue(iw,is) = mean(yvec);

    plot([xvec;xvec(1)]-mean(xvec)+cxNew,[yvec;yvec(1)]-mean(yvec)+cyNew,'r','linewidth',2)
    xT = [xvec;xvec(1)]-mean(xvec)+cxNew;
    yT = [yvec;yvec(1)]-mean(yvec)+cyNew;
    hFill = fill(xT, yT,'r');
    set(hFill,'EdgeColor','r')


    figure(2);
    xvec = XfDNN{iw,is}(1:end/2);
    yvec = XfDNN{iw,is}(end/2+1:end);
    
    cxDNN(iw,is) = mean(xvec);
    cyDNN(iw,is) = mean(yvec);

    plot([xvec;xvec(1)]-mean(xvec)+cxNew,[yvec;yvec(1)]-mean(yvec)+cyNew,'r','linewidth',2)
    xT = [xvec;xvec(1)]-mean(xvec)+cxNew;
    yT = [yvec;yvec(1)]-mean(yvec)+cyNew;
    hFill = fill(xT, yT,'r');
    set(hFill,'EdgeColor','r')


  end
end

figure(1);
axis equal
xticks(([0:numel(chanWidths)-1]*0.5))
xticklabels({'0.2','0.4','0.6','0.75'})
yticks(([0;0.5]))
yticklabels({'20','25'})
box on

figure(2);
axis equal
xticks(([0:numel(chanWidths)-1]*0.5))
xticklabels({'0.2','0.4','0.6','0.75'})
yticks(([0;0.5]))
yticklabels({'20','25'})
box on

%% 
figure(3); clf; 
plot(CnVal(:,1),cyTrue(:,1)./chanWidths,'k-s','linewidth',2,'markersize',10,'markerfacecolor','k')
hold on
plot(CnVal(:,1),cyDNN(:,1)./chanWidths,'r-o','linewidth',2,'markersize',10,'markerfacecolor','r')
ylim([-0.05 0.15])
xlim([0.1 0.9])
axis square
grid 
lgd = legend('Numerical','Network');
lgd.Location = 'northoutside';
lgd.Orientation = 'horizontal';
lgd.Box = 'off';

figure(4); clf; 
plot(CnVal(:,2),cyTrue(:,2)./chanWidths,'k-s','linewidth',2,'markersize',10,'markerfacecolor','k')
hold on
plot(CnVal(:,2),cyDNN(:,2)./chanWidths,'r-o','linewidth',2,'markersize',10,'markerfacecolor','r')
ylim([-0.05 0.15])
xlim([0.1 0.9])
axis square
grid 
lgd = legend('Numerical','Network');
lgd.Location = 'northoutside';
lgd.Orientation = 'horizontal';
lgd.Box = 'off';


