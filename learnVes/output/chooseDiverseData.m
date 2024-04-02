clear; clc;

addpath ../../src/

oc = curve;

fileNames{1} = './output/nv70N32DNN5thOrder_Kb2e1Dt2E5_VF30_bgFlowcouette_speed100.mat';
fileNames{2} = './output/nv81N32DNN5thOrderNEWNETS_Kb2e1Dt2E5_VF35_bgFlowcouette_speed100.mat';
fileNames{3} = './output/initialManyVesRuns/nv70N32DNN32PCAmodes_wNewDistJiggVF30_bgFlowcouette_speed100.mat';
fileNames{4} = './output/initialManyVesRuns/nv47N32DNNprobICwJiggVF20_bgFlowcouette_speed100.mat';

fileNames{5} = './output/initialManyVesRuns/nv47N32DNN32PCAmodes_wNewDistJiggVF20_bgFlowcouette_speed200.mat';
fileNames{6} = './output/initialManyVesRuns/nv47N32DNN32PCAmodes_wNewDistJiggVF20_bgFlowcouette_speed200.mat';
fileNames{7} = './output/initialManyVesRuns/nv24N64VF10TrueLoadIC_bgFlowcouette_speed100.mat';
fileNames{8} = './output/initialManyVesRuns/nv24N64VF10TrueLoadIC_bgFlowcouette_speed100.mat';

fileNames{9} = './output/newNetRuns/nv70N32DNN5thOrder_Kb2e1Dt2E5_VF30_bgFlowcouette_speed100.mat';
fileNames{10} = './output/newNetRuns/nv81N32DNN5thOrderNEWNETS_Kb2e1Dt2E5_VF35_bgFlowcouette_speed100.mat';

fileNames{11} = './output/paperDataFiles/nv70N64DNNwNewJiggNewLP_VF30NoLoop_bgFlowcouette_speed100.mat';
fileNames{12} = './output/paperDataFiles/nv70N64VF30GT_bgFlowcouette_speed100.mat';
fileNames{13} = './output/paperDataFiles/nv81N32DNNVF35UseRegulNear_hHalfTol_bgFlowcouette_speed100.mat';
fileNames{14} = './output/paperDataFiles/nv81N32DNNwNewJiggNewLP_VF35_NoLoop_bgFlowcouette_speed100.mat';
fileNames{15} = './output/paperDataFiles/nv81N32DNNwNewJiggNewLP_VF35UseRegulNear_NoLoop_bgFlowcouette_speed100.mat';
fileNames{16} = './output/paperDataFiles/nv81N48VF35TrueLoadIC_bgFlowcouette_speed100.mat';
fileNames{17} = './output/paperDataFiles/nv70N32DNNwNewJiggNewLP_VF30NoLoop_bgFlowcouette_speed100.mat';

nsteps = 20;
shapeTol = 1e-2;


vesCount = 0;
Xstore = [];

%appTenStore = {};
%trueTenStore = {};
%trueVselfStore = {};
%trueVnearStore = {};
%trueVfarStore = {};
%appVselfStore = {};
%appVnearStore = {};
%appVfarStore = {};
%vBackStore = {};

for irun = 2 : numel(fileNames)
  file1 = fileNames{irun};
  load(file1)
  nsteps = numel(time);
  nves = size(Xhist,2);

  for itime = 1:1000:nsteps
    for ives = 1 : nves
      X = Xhist(:,ives,itime);   
      if norm(X)>0
      [Xstand,~,~,~,~] = standardizationStep(X,oc);
        
      iadd = false;
      if vesCount == 0
        iadd = true;
      else
        if vesCount < 11
          err = zeros(vesCount,1);
          for is = 1 : vesCount
            err(is) = hausdorfDistance(Xstand,Xstore(:,is));
          end
        else
          err = zeros(10,1); idx = 1;
          for is = vesCount-9 : vesCount
            err(idx) = hausdorfDistance(Xstand,Xstore(:,is));
            idx = idx + 1;
          end
        end
        if min(err) >= shapeTol
          iadd = true;
        end
      end % vesCount

      if iadd
        vesCount = vesCount+1;
        Xstore(:,vesCount) = Xstand;
        figure(1); clf;
        plot(Xstand(1:end/2),Xstand(end/2+1:end))
        axis equal
        pause(0.1)
      end % iadd

      if rem(vesCount,1000) == 0
        disp('Saving data set...')
        save('diverseNetworkSet','Xstore','-v7.3') 
      end
      end
    end %ives
    disp([num2str(itime) 'th time step in ' num2str(irun) 'th run is done'])
    disp(['Total number of different vesicle configs is ' num2str(vesCount)])
  end %itime
end % irun

disp('Saving data set...')
save('diverseNetworkSet','Xstore','-v7.3') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hausErr = hausdorfDistance(X1,X2)
N = numel(X1)/2;

% Closest points on X2 to the points on X1
d1to2 = zeros(N,1);
for in = 1 : N
  d1to2(in) = min(((X1(in)-X2(1:end/2)).^2 + (X1(in+N)-X2(end/2+1:end)).^2).^0.5)/...
      sqrt(X1(in).^2+X1(in+N).^2);
end

hausErr = max(d1to2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,scaling,rotate,trans,sortIdx] = standardizationStep(Xin,oc)
Xin = [interpft(Xin(1:end/2),128);interpft(Xin(end/2+1:end),128)];
N = numel(Xin)/2;
X = Xin;
% Equally distribute points in arc-length
for iter = 1 : 3
  [X,~,~] = oc.redistributeArcLength(X);
end
% Fix misalignment in center and angle due to reparametrization
X = oc.alignCenterAngle(Xin,X);

% standardize angle, center, scaling and point order
[trans,rotate,scaling,sortIdx] = referenceValues(X,oc);
X = standardize(X,trans,rotate,scaling,sortIdx,N);
end % standardizationStep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XrotSort = standardize(X,translation,rotation,scaling,sortIdx,N)

% translate, rotate and scale configuration
Xrotated = scaling*rotationOperator(translateOp(X,translation),rotation);   

% now order the points
XrotSort = [Xrotated(sortIdx);Xrotated(sortIdx+N)];

end % standardize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [translation,rotation,scaling,sortIdx] = referenceValues(Xref,oc)

N = numel(Xref)/2;

% find translation, rotation and scaling
translation = [-mean(Xref(1:end/2));-mean(Xref(end/2+1:end))];
rotation = pi/2-oc.getIncAngle2(Xref);
    
% amount of scaling
[~,~,length] = oc.geomProp(Xref);
scaling = 1/length;
    
% find the ordering of the points
Xref = scaling*rotationOperator(translateOp(Xref,translation),rotation);

firstQuad = find(Xref(1:end/2)>=0 & Xref(end/2+1:end)>=0);
theta = atan2(Xref(end/2+1:end),Xref(1:end/2));
[~,idx]= min(theta(firstQuad));
sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];

end % referenceValues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xrot = rotationOperator(X,theta)
% Get x-y coordinates
Xrot = zeros(size(X));
x = X(1:end/2); y = X(end/2+1:end);

% Rotated shape
xrot = (x)*cos(theta) - (y)*sin(theta);
yrot = (x)*sin(theta) + (y)*cos(theta);

Xrot(1:end/2) = xrot;
Xrot(end/2+1:end) = yrot;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateOp(X,transXY)
Xnew = zeros(size(X));
Xnew(1:end/2) = X(1:end/2)+transXY(1);
Xnew(end/2+1:end) = X(end/2+1:end)+transXY(2);

end
    
   
