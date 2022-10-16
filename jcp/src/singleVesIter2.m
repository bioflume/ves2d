clear all; clc

prams.kappa = 1e0;  
prams.order = 1; 
prams.Case = 'schurComp';

options.usePlot = 0;                    
options.progressBar = 0;
options.verbose = 0;
options.GmresTol = 10^-4;    
load('../svo1.mat')
svo1 = svo1(:,1:4);
svo1 = svo1([1:20 1:40],:);
svo1(1:20,4) = 9;

load('../../jcp/results/iterCount.mat')
TMavgIter(:,[4 5]) = 0;
LHSavgIter(:,[4 5]) = 0;
TsaveFlag = 1; 
LsaveFlag = 1;
save('../../jcp/results/iterCount.mat','TMavgIter','LHSavgIter','TsaveFlag','LsaveFlag')

for ii=1:60
  shear = svo1(ii,1);
  vc = svo1(ii,2);
  n = svo1(ii,3);
  m = svo1(ii,4); 
  
  prams.T = 6/m/shear;
  prams.vInf = @(X) farFieldVel(X,'shear',shear);
  prams.viscCont = vc;
  prams.n = n;  
  
  load('../../jcp/results/iterCount.mat')
  TsaveFlag = 1;
  LsaveFlag = 1;
  save('../../jcp/results/iterCount.mat','TMavgIter','LHSavgIter','TsaveFlag','LsaveFlag')
  
  X = boundary(prams.n);
  clear functions global;
  [Xfinal status] = Ves2D(X,prams,options); 
  disp([ii shear vc n m status.flag]);
end