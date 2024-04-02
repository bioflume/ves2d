clear; clc;
addpath ../src/

oc = curve;

load ./necessaryMatFiles/X100KinitShapes.mat
Xstore1 = Xstore;

load ./output/diverseNewDataSet.mat
Xstore1 = [interpft(Xstore1(1:end/2,:),128);interpft(Xstore1(end/2+1:end,:),128)];
Xstore = [Xstore1 Xstore];


% number of vesicles
nves = size(Xstore,2);
% num. points per vesicle
N = 128;

XstandStore = [];
cnt = 0;
dnn = dnnTools;
% Increase data set by the mirror images and pi degree rotations
for ives = 1 : nves
  disp(ives)
  % Standardize if not done yet
  Xinit = Xstore(:,ives);
  [Xinit,scaling,rotate,trans,sortIdx] = dnn.standardizationStep(Xinit,oc);
   
  Xinit2other = zeros(2*N,2);

  Xinit2other(:,1) = [-Xinit(1:end/2);Xinit(end/2+1:end)];

  Xinit2other(1:end/2,2) = Xinit(1:end/2)*cos(pi) - Xinit(end/2+1:end)*sin(pi);
  Xinit2other(end/2+1:end,2) = Xinit(end/2+1:end)*sin(pi) + Xinit(end/2+1:end)*cos(pi);

  % Find ordering for XinitMirr and XinitRot
  
  for k = 1 : 2
    cnt = cnt + 1;
    Xref = Xinit2other(:,k);
    firstQuad = find(Xref(1:end/2)>=0 & Xref(end/2+1:end)>=0);
    theta = atan2(Xref(end/2+1:end),Xref(1:end/2));
    [~,idx]= min(theta(firstQuad));
    sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];
    XstandStore(:,cnt) = [Xref(sortIdx);Xref(sortIdx+N)];
  end
  disp(['cnt = ' num2str(cnt)])
end
