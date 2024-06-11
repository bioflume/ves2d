clear; clc;

addpath ../../src/

oc = curve;

shapeTol = 1e-2;

vesCount = 0;
%Xstore = [];

load diverseParabolic1DataSet.mat
% load Xstore
vesCount = numel(Xstore(1,:));

% fileIDs = (1:150)';
fileICs = (1:5)';
fileIVs = (1:9)';

fileNameDir = './simpleFlowData/parabolicRuns/';

% for irun = 1 : numel(fileIDs)
for ic = 1 : numel(fileICs)
  for iv = 1 : numel(fileIVs)

  file1 = [fileNameDir 'RunIC' num2str(ic) '_IV' num2str(iv) '.bin'];
  [vesxN, vesyN, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(file1);
  
  nsteps = numel(time);

  for itime = 1:100:nsteps
    
    X = [vesxN(:,itime);vesyN(:,itime)];   
    if norm(X)>0
      [Xstand,~,~,~,~] = standardizationStep(X,oc);
        
      iadd = false;
      if vesCount == 0
        iadd = true;
      else
        if vesCount < 1000
          err = zeros(vesCount,1);
          for is = 1 : vesCount
            err(is) = hausdorfDistance(Xstand,Xstore(:,is));
          end
        else
          err = zeros(1000,1); idx = 1;
          for is = vesCount-999 : vesCount
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
        save('diverseParabolicDataSet','Xstore','-v7.3') 
      end
    end
   
    disp([num2str(itime) 'th time step in ' num2str(iv*ic) 'th run is done'])
    disp(['Total number of different vesicle configs is ' num2str(vesCount)])
  end %itime
  end % irun
end

disp('Saving data set...')
save('diverseParabolicDataSet','Xstore','-v7.3') 

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
for iter = 1 : 10
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
    
