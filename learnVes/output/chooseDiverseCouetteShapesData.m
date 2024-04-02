clear; clc;

addpath ~/codes/Ves2Dn/src/
addpath /workspace/gokberk/couette150VesData

oc = curve;

fileNames{1} = 'couette150VelData_Set1Part';
fileNames{2} = 'couette150VelData_Set2Part';
fileNames{3} = 'couette150VelData_Set3Part';
fileNames{4} = 'couette150VelData_Set4Part';
fileNames{5} = 'couette150VelData_Set5Part';


shapeTol = 1e-2;

% load current set Xstore, XorigStore, vesCount
load diverseMoreCouetteRuns1Data
vesCount = numel(Xstore(1,:));

%vesCount = 0;
%Xstore = [];
%XorigStore = [];

for irun = 2 : numel(fileNames)
  file1 = [fileNames{irun} '1.mat'];
  file2 = [fileNames{irun} '2.mat'];
  
  Xshapes = zeros(64,150,350);

  load(file1)
  Xshapes(:,:,1:175) = Xparam;
  
  
  load(file2)
  Xshapes(:,:,176:350) = Xparam;
  
  
  for itime = 1:350
    for ives = 1 : 150
      X = Xshapes(:,ives,itime);    
      [translation,rotation,scaling,sortIdx] = referenceValues(X,oc);
      Xstand = standardize(X,translation,rotation,scaling,sortIdx,32);
      % Equi-distribute in arc-length
      for iter = 1 : 5
        Xstand = oc.redistributeArcLength(Xstand);
      end
        
      iadd = false;
      if vesCount == 0
        iadd = true;
      else
        err = zeros(vesCount,1);
        for is = 1 : vesCount
          err(is) = hausdorfDistance(Xstand,Xstore(:,is));
        end
        
        if min(err) >= shapeTol
          iadd = true;
        end
      end % vesCount

      if iadd
        vesCount = vesCount+1;
        Xstore = [Xstore Xstand];
        XorigStore = [XorigStore X];
      end % iadd
    end %ives
    disp([num2str(itime) 'th time step in ' num2str(irun) 'th run is done'])
    disp(['Total number of different vesicle configs is ' num2str(vesCount)])
    if rem(vesCount,5000)==0
      save('diverseMoreCouetteData','Xstore','XorigStore','vesCount')    
    end
  end %itime
end % irun

disp('Saving data set...')
save('diverseMoreCouetteData','Xstore','XorigStore','-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hausErr = hausdorfDistance(X1,X2)
N = numel(X1)/2;

x1 = X1(1:end/2); y1 = X1(end/2+1:end);
x2 = X2(1:end/2); y2 = X2(end/2+1:end);

xx1 = x1(:,ones(N,1)); yy1 = y1(:,ones(N,1));
xx2 = x2(:,ones(N,1))'; yy2 = y2(:,ones(N,1))';

d1to2 = min(sqrt((xx1-xx2).^2+(yy1-yy2).^2),[],2)./sqrt(x1.^2+y1.^2);

hausErr = max(d1to2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XrotSort = standardize(X,translation,rotation,scaling,sortIdx,N)

% translate, rotate and scale configuration
Xrotated = scaling*rotationOperator(translateOp(X,translation),rotation);   

% now order the points
sorting = [(sortIdx(1):N)';(1:sortIdx(1)-1)'];
XrotSort = [Xrotated(sorting);Xrotated(sorting+N)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [translation,rotation,scaling,sortIdx] = referenceValues(Xref,oc)
% find translation, rotation and scaling
translation = [-mean(Xref(1:end/2));-mean(Xref(end/2+1:end))];
rotation = pi/2-oc.getIncAngle(Xref);
    
% amount of scaling
[~,~,length] = oc.geomProp(Xref);
scaling = 1/length;
    
% find the ordering of the points
XrefRot = scaling*rotationOperator(translateOp(Xref,translation),rotation);
theta = atan2(XrefRot(end/2+1:end),XrefRot(1:end/2));
[~,sortIdx] = sort(theta);
end
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
    
   
