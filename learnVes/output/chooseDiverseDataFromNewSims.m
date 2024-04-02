clear; clc;

addpath ../../src/

oc = curve;

shapeTol = 1e-2;

vesCount = 0;
Xstore = [];

fileIDs = (1:55)';

fileNameDir = './newDataGenRuns/dataGenSpeed_Id';

for irun = 1 : numel(fileIDs)
  file1 = [fileNameDir num2str(irun) '.bin'];
  [Xhist, time, ~, N, nves, ~, ~, ~, ~, ~] = loadNewData(file1);
  nsteps = numel(time);

  for itime = 1:100:nsteps
    for ives = 1 : nves
      X = Xhist(:,ives,itime);   
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
        save('diverseNetworkSet','Xstore','-v7.3') 
      end
      end
    end %ives
    disp([num2str(itime) 'th time step in ' num2str(irun) 'th run is done'])
    disp(['Total number of different vesicle configs is ' num2str(vesCount)])
  end %itime
end % irun

disp('Saving data set...')
save('diverseNewDataSet','Xstore','-v7.3') 

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function [X, time, sigma, N, nv, XwallsExt, XwallsInt, etaExt, etaInt, RS] = loadNewData(fileName)

fid = fopen(fileName,'r');
val = fread(fid,'double');
fclose(fid);

% Some important numbers
dt = val(1);
nv = val(2);
N = val(3);
nvbdExt = val(4);
NbdExt = val(5);
nvbdInt = val(6);
NbdInt = val(7);

% Remove them and move on with rigid boundaries
val = val(8:end);
XwallsExt = val(1:2*NbdExt);
val = val(2*NbdExt+1:end);
XwallsInt = val(1:2*nvbdInt*NbdInt);
XwallsInt = reshape(XwallsInt,[2*NbdInt nvbdInt]);
val = val(2*nvbdInt*NbdInt+1:end);

% Now read in time steps
X = []; sigma = []; time = [];
etaExt = []; etaInt = []; RS = [];
nEntry2Read = 1 + 3 * N * nv + 2 * NbdExt + 2 * nvbdInt * NbdInt + 3*nvbdInt;
ist = 1;
while numel(val)>=nEntry2Read

  time(ist) = val(1);
  Xst = reshape(val(2:2*N*nv+1),[2*N nv]);
  X(:,:,ist) = Xst;

  val = val(2*N*nv+2:end);
  sigma(:,:,ist) = reshape(val(1:N*nv), [N nv]);

  val = val(N*nv+1:end);
  etaExt(:,ist) = val(1:2*NbdExt);
  val = val(2*NbdExt+1:end);
  etaInt(:,:,ist) = reshape(val(1:2*nvbdInt*NbdInt), [2*NbdInt nvbdInt]);
  val = val(2*NbdInt*nvbdInt+1:end);
  RS(:,:,ist) = reshape(val(1:3*nvbdInt+3), [3 nvbdInt+1]);
  val = val(3*nvbdInt+4:end);

  ist = ist + 1;
end


end