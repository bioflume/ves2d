addpath ../src/

postProccessing = true;
n = 128; nv = 10; reducedArea = .65; flowCurv = 20;

prams.T = 100;                         
prams.ts = 1e-3;
prams.m = prams.T/prams.ts;

prams.kappa = 1e-1;
prams.viscCont = ones(1,nv);
prams.order = 1;   
prams.Incompressibility = 1;
prams.vInf = @(X) farFieldVel(X,'invParabolic','k',flowCurv);

options.usePlot = true;                        
options.progressBar = 0;                      
options.saveData = 1;                       
options.dataStride = ceil(prams.m/500);     
options.useGPU = false;                         
options.showError = true;                   
options.GmresTol = 1e-8;                    
options.verbose = true;

fileName = ['./fileOfVes_' num2str(nv) '_' num2str(100*reducedArea) '_' ...
            num2str(flowCurv) '_' num2str(prams.T)  '_' num2str(1000*prams.ts) ...
            '_' num2str(n)];
options.fileName = [fileName '.bin'];      

%-If the file exists continue from file
if ( exist(options.fileName) )
  fileId = fopen(options.fileName,'r');
  Result = fread(fileId,'double');
  fclose(fileId);
 
  Result = reshape(Result,5*nv*n+1,[]); 
  X = Result(1:2*nv*n,:);
  X = reshape(X(:,end),[],nv);

  clear Result;
else
  X = boundary(n,'nv',nv,'scale',1,'reducedArea', reducedArea);  
  X = makeGridUniform(X);
end

if (~postProccessing)
  [trash prams options] = monitor(prams,options);
  viewer(X,[],prams,options);
 
  save([fileName '.mat']);
  [XF status] = Ves2D(X,prams,options,@monitor);
else

  load(fileName);
  
  fileId = fopen(options.fileName,'r');
  Result = fread(fileId,'double');
  fclose(fileId);

  Result = reshape(Result,5*nv*n+1,[]); 
  
  Xv   = Result(1       :2*nv*n,:);
  Time = Result(5*nv*n+1       ,:);
      
  [trash params options] = monitor(prams,options);
  options.axisOn=1;
  for ii=1:size(Xv,2)
    viewer(reshape(Xv(:,ii),[],nv),[],params,options);
    title(num2str(ii*options.dataStride*params.ts));
    pause(.2);
  end
end