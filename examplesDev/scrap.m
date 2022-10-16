clear all; clc;
%dbstop in updateFirstM at 92
n = 64;
nv = 1;
% Setting modeling parameters
prams.n = n;
prams.nv = nv;
prams.T = 180;
prams.m = 1800;
prams.kappa = 1e-1;
prams.viscCont = ones(nv,1);
prams.order = 2;   
prams.ts = prams.T/prams.m;
prams.Incompressibility = 1;
prams.flowType = 'confined';
prams.M = [128 64];         
prams.bd = @(ind,m) sampleBd(ind,m,1,'couette');
prams.bc = @(x) forcing(x,'couette');
prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,prams);


options.usePlot = 0;
options.progressBar = 0;
options.AxesHandle = gca;
options.saveData = 1;
options.dataStride = 20;

[trash prams options] = monitor(prams,options);
%% 
%% scaleRange = 2:2:20;
%% for ii = scaleRange
%%   disp(ii);
%%   clear functions global;
%%   scale = ii/10;
%%   options.fileName = ['../docs/testRuns/vesSize' num2str(ii) '.bin'];
%%   X = boundary(n,'couette','nv',nv,'scale',scale);
%%   status(ii).flag = []; status(ii).fileName = [];
%%   [XF status(ii)] = Ves2D(X,prams,options,@monitor);
%% end
%% 
%% scaleRange = 2:2:30;
%% for ii = scaleRange
%%   disp(ii);
%%   clear functions global;
%%   scale = ii/10;
%%   options.fileName = ['../docs/testRuns/vesSizeTilted' num2str(ii) '.bin'];
%%   X = boundary(n,'couette','nv',nv,'scale',scale,'angle',pi/2);
%%   [XF statusTilted(ii)] = Ves2D(X,prams,options,@monitor);
%% end

% Reading the saved data from the file.
Result = []; jj = 1;
for ii = 6:2:20
  fileName = ['../docs/testRuns/vesSize' num2str(ii) '.bin'];
  fileId = fopen(fileName,'r');
  R = fread(fileId,'double');
  fclose(fileId);

  R = reshape(R,321,[]); 
  Result(1:2*n,:,jj) = R(1:2*n,1:end);
  jj = jj+1;
end

Result = permute(Result,[1 3 2]);

for jj = 1:2:size(Result,3)
  pause(.1);cla;
  viewer(Result(1:2*n,:,jj),[],prams,options);       
  figName = ['../docs/testRuns/vesSize' num2str(jj) '.jpg'];
  saveas(gca,figName);
end


Result = []; jj = 1;
for ii = 6:2:30
  fileName = ['../docs/testRuns/vesSizeTilted' num2str(ii) '.bin'];
  fileId = fopen(fileName,'r');
  R = fread(fileId,'double');
  fclose(fileId);

  R = reshape(R,321,[]); 
  Result(1:2*n,:,jj) = R(1:2*n,1:end);
  jj = jj+1;
end

Result = permute(Result,[1 3 2]);

for jj = 1:2:size(Result,3)
  pause(.1);cla;
  viewer(Result(1:2*n,:,jj),[],prams,options);       
  figName = ['../docs/testRuns/vesSizeTilted' num2str(jj) '.jpg'];
  saveas(gca,figName);
end
