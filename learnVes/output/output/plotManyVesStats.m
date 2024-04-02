clear; clc;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

% if file does not exist, leave it []
% fileDNN = 'nv47N32DNNwNewJiggNewLP_VF20NoLoop_bgFlowcouette_speed100';
% fileDoF = 'nv47N32VF20DOF_bgFlowcouette_speed100';
% fileCost = 'nv47N16VF20Cost_bgFlowcouette_speed100';
% fileGT = 'nv47N64VF20TrueLoadIC_bgFlowcouette_speed100';

% fileDNN = 'nv70N32DNNwNewJiggNewLP_VF30NoLoop_bgFlowcouette_speed100';
% fileDoF = 'nv70N32VF30TrueLoadIC_bgFlowcouette_speed100';
% fileCost = [];
% fileGT = 'nv70N64VF30GT_bgFlowcouette_speed100.mat';

% fileDNN = 'nv81N32DNNwNewJiggNewLP_VF35_NoLoop_bgFlowcouette_speed100';
fileDNN = 'nv81N32DNNwFMMmoreNwallsMoreTol_VF35_bgFlowcouette_speed100';
fileDoF = 'nv81N32VF35DoF_bgFlowcouette_speed100';
fileCost = [];
fileGT = 'nv81N48VF35TrueLoadIC_bgFlowcouette_speed100';

% load files
if ~isempty(fileDNN)
  load(fileDNN);
  ntime = numel(time); nv = numel(Xhist(1,:,1));
  cxDNN = zeros(nv,ntime); cyDNN = cxDNN; crDNN = cxDNN;
  for k = 1 : ntime
    x = interpft(Xhist(1:end/2,:,k),256);
    y = interpft(Xhist(end/2+1:end,:,k),256);
    cxDNN(:,k) = mean(x)'; cyDNN(:,k) = mean(y)';
    crDNN(:,k) = sqrt(cxDNN(:,k).^2+cyDNN(:,k).^2);
  end
  crDNNrange = linspace(1,2.2,500);
  pdfCrDNN = ksdensity(crDNN(:),crDNNrange);
else
  crDNNrange = [];
  pdfCrDNN = [];
end


if ~isempty(fileDoF)
  load(fileDoF);
  ntime = numel(timeTrue); nv = numel(XhistTrue(1,:,1));
  cxDoF = zeros(nv,ntime); cyDoF = cxDoF; crDoF = cxDoF;
  for k = 1 : ntime
    x = interpft(XhistTrue(1:end/2,:,k),256);
    y = interpft(XhistTrue(end/2+1:end,:,k),256);
    cxDoF(:,k) = mean(x)'; cyDoF(:,k) = mean(y)';
    crDoF(:,k) = sqrt(cxDoF(:,k).^2+cyDoF(:,k).^2);
  end
  crDoFrange = linspace(1,2.2,500);
  pdfCrDoF = ksdensity(crDoF(:),crDoFrange);
else
  crDoFrange = [];
  pdfCrDoF = [];
end


if ~isempty(fileCost)
  load(fileCost);
  ntime = numel(timeTrue); nv = numel(XhistTrue(1,:,1));
  cxCost = zeros(nv,ntime); cyCost = cxCost; crCost = cxCost;
  for k = 1 : ntime
    x = interpft(XhistTrue(1:end/2,:,k),256);
    y = interpft(XhistTrue(end/2+1:end,:,k),256);
    cxCost(:,k) = mean(x)'; cyCost(:,k) = mean(y)';
    crCost(:,k) = sqrt(cxCost(:,k).^2+cyCost(:,k).^2);
  end
  crCostrange = linspace(1,2.2,500);
  pdfCrCost = ksdensity(crCost(:),crCostrange);
else
  crCostrange = [];
  pdfCrCost = [];
end

if ~isempty(fileGT)
  load(fileGT);
  ntime = numel(timeTrue); nv = numel(XhistTrue(1,:,1));
  cxGT = zeros(nv,ntime); cyGT = cxGT; crGT = cxGT;
  for k = 1 : ntime
    x = interpft(XhistTrue(1:end/2,:,k),256);
    y = interpft(XhistTrue(end/2+1:end,:,k),256);
    cxGT(:,k) = mean(x)'; cyGT(:,k) = mean(y)';
    crGT(:,k) = sqrt(cxGT(:,k).^2+cyGT(:,k).^2);
  end
  crGTrange = linspace(1,2.2,500);
  pdfCrGT = ksdensity(crGT(:),crGTrange);
else
  crGTrange = [];
  pdfCrGT = [];
end

% plot
figure(1); clf; hold on;
if ~isempty(pdfCrGT)
  plot(crGTrange,pdfCrGT,'k','linewidth',2)
end
if ~isempty(pdfCrDoF)
  plot(crDoFrange,pdfCrDoF,'Color',[0 0.45 0.74],'linewidth',2)    
end
if ~isempty(pdfCrCost)
  plot(crCostrange,pdfCrCost,'Color',[0 0.5 0],'linewidth',2)
end
if ~isempty(pdfCrDNN)
  plot(crDNNrange,pdfCrDNN,'r','linewidth',2)
end

axis square
grid on
box on
xlabel('$\| \mathbf{c} \|$')
ylabel('$p(\| \mathbf{c} \|)$')
if ~isempty(pdfCrGT) && ~isempty(pdfCrDoF) && ~isempty(pdfCrCost) && ~isempty(pdfCrDNN)
legend('True','Same-DOF','Same-cost','MLARM','location','northwest')
legend boxoff
end
xlim([1 2.2])
    