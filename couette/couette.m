addpath ../src/

if(~any(strcmp(who,'RandSeed')) || isempty(RandSeed))
  RandSeed = 10;
end

rand('state',RandSeed);

vesDist = 'volFrac'; 
vesSize = 1;         
useGPU = false;

prams.T = 200;                               
prams.ts = ts;                              
prams.kappa = 1e-1;                            
prams.n = n;                     
prams.m = prams.T/prams.ts;                   
prams.order = 1;                    

vf = [2 3 4  5  6  7  8  9  10]/100;
nv = [6 9 12 16 18 21 24 27 30];
nv = nv(vf == volFrac);
if( isempty(nv) ), nv = round(3 * volFrac * 100);end

%%-- Setting up the boundary
Ri = 10; Ro = 20;             
omegaIn = 1.0; omegaOut = 0;        
prams.flowType = 'confined';             
prams.M = 2*[128 64];                     
prams.bd = @(ind,m) sampleBd(ind,m,1,'couette','Ri',Ri,'Ro',Ro);
prams.bc = @(x) forcing(x,'couette','Ri',Ri,'Ro',Ro,'omegaIn',omegaIn,'omegaOut',omegaOut);;
prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,prams,'direct',useGPU);

if( viscCont == 1 )
  viscString = '';
else
  viscString = ['_' num2str(viscCont)];
end

fileName = ['./couetteFlow_' num2str(100*volFrac) '_' num2str(100*reducedArea) ...
            viscString '_' num2str(prams.T) '_' num2str(1000*prams.ts) '_' ...
            num2str(prams.n)];

options.usePlot = 0;                          
options.progressBar = 0;                      
options.saveData = 1;                         
options.dataStride = ceil(prams.m/1000);      
options.fileName = [fileName '.bin'];         
options.useGPU = useGPU;                             
options.showError = true;                     
options.GmresTol = 1e-8;                     
options.verbose = true;

%-If the file exists continue from file
if ( exist(options.fileName) )
  fileId = fopen(options.fileName,'r');
  Result = fread(fileId,'double');
  fclose(fileId);
 
  Result = reshape(Result,5*nv*prams.n+1+2*sum(prams.M)+3,[]); 
  X = Result(1:2*nv*prams.n,:);
  X = reshape(X(:,end),[],nv);

  clear Result;
else
  %%-- Generating the boundary
  domain = fixedBound(prams.M,prams.bd,1);
  if (nv < 22)
    X = boundary(prams.n,'couette',domain,vesDist,volFrac,'scale',vesSize, ...
                 'reducedArea',reducedArea);
  else
    ds = 4.5;
    allCents = [];
    angle = [];
    
    rBase = (Ri + Ro)/2;
    nEach = floor(2*pi*rBase/ds);
    nTier = ceil(nv/nEach);
    rBase = Ri + (1:nTier) * (Ro - Ri)/(nTier+1);

    nEach = zeros(1,nTier);
    nRemain = nv;
    for kk=nTier:-1:1
      nEach(kk) = floor(2*pi*rBase(kk)/ds);
      if(nEach(kk) > nRemain), nEach(kk) = nRemain; end
      nRemain = nRemain - nEach(kk);
    end
    
    for kk=1:nTier
      tt = linspace(0,2*pi,nEach(kk)+1);
      tt = tt(1:end-1) + (rand(1,nEach(kk))- .5)*2*pi/nEach(kk)/5;
      rr = rBase(kk) + (rand(1,nEach(kk))- .5)*(Ro-Ri)/(nTier+1)/10;
      cen = zeros(2,nEach(kk));
      [cen(1,:) cen(2,:)] = pol2cart(tt,rr);
      allCents = [allCents cen];
      angle = [angle pi/2+tt];
     end
    X = boundary(prams.n,'nv',nv,'scale',vesSize, 'reducedArea', ...
                 reducedArea,'center',allCents,'angle',angle);
    clear tt rr cen angle allCents nTier nEach;

%     tt = linspace(0,2*pi,nv+1);tt=tt(1:end-1);
%     rr = (Ri + Ro)/2 + (rand(1,nv)-.5)*(Ro-Ri)*.6;
%     cen = zeros(2,nv);
%     [cen(1,:) cen(2,:)] = pol2cart(tt,rr);
    
%     X = boundary(prams.n,'nv',nv,'scale',vesSize, 'reducedArea', ...
%                  reducedArea,'center',cen,'angle',pi/3);
%     clear tt rr cen;
  end
  
  X = makeGridUniform(X);
end

[trash prams options] = monitor(prams,options);
viewer(X,[],prams,options);
prams.viscCont = viscCont * ones(1,nv);                           
 
save([fileName '.mat']);
[XF status] = Ves2D(X,prams,options,@monitor);