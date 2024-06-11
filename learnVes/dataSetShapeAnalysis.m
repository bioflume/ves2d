clear; clc;

load ./workingOnDataSet/newStand625K_n128Dt1e-05RelaxDataSet_IT3

Xstore = XstandStore(:,1:nInstances/4);
nInstances = nInstances/4;

addpath ../src/
op = poten(N);
oc = curve;

bendEnergy = zeros(nInstances,1);
condL = zeros(nInstances,1);
meanCurv = zeros(nInstances, 1);

for k = 1 : nInstances
  disp(k)
  X = Xstore(:,k);
  [X,~] = oc.reparametrize(X,[],6,20);
  
  vesicle = capsules(X,[],[],1,1,0);
  G = op.stokesSLmatrix(vesicle);
  [Ben, Ten, Div] = vesicle.computeDerivs;
  L = Div*G*Ten;
  condL(k) = cond(L);

  [jac, tan, curv] = oc.diffProp(X);
  bendEnergy(k) = sum(curv.^2 .* jac) * pi/N;
  meanCurv(k) = sum(curv.*jac) * 2 * pi /N;
end
