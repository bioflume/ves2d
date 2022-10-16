addpath ../src/ 

ReducedArea = [.6 .7 .8 .9];
ViscCont = [.25 .5 1 2 4 8 16];
%ShearRate = [.01 .1 1 10];

np          = 32;                     
prams.order = 1;                              
prams.kappa = 1;                            

options.usePlot = 1;                          
options.progressBar = 0;                      
options.saveData = 1;                         
options.useGPU = false;                             
options.showError = true;                     
options.GmresTol = 1e-8;                     
options.verbose = true;


for sr = ShearRate
  
  prams.T  = 50/sr;                               
  prams.m = 50000;
  prams.ts = prams.T/prams.m;                              
  options.dataStride = ceil(prams.m/500);      
  
  for ra = ReducedArea
    for vc = ViscCont

      clear functions;

      fileName = ['InclinationAngle_' num2str(sr) '_' num2str(ra) '_' num2str(vc)];
      options.fileName = [fileName '.bin'];         
      
      %- Background flow
      prams.vInf = @(X) farFieldVel(X,'shear',sr);
      
      %- Generating the boundary
      X = boundary(np,'scale',1,'reducedArea',ra);
      X = makeGridUniform(X);
      [Ra A L] = reducedVolume(X);

      prams.viscCont = vc;              
      NonDimShear = 1 * (L / 2 / pi)^3 / prams.kappa;

      save([fileName '.mat']);
      [XF status] = Ves2D(X,prams,options,@monitor);    
      save([fileName '.mat']);
    end
  end
end

% %- post process
% for sr = ShearRate
%   for ra = ReducedArea
%     for vc = ViscCont
%       fileName = ['InclinationAngle_' num2str(sr) '_' num2str(ra) '_' num2str(vc)];
      
%       %load(fileName);
%       fileId = fopen([fileName '.bin'],'r');
%       Result = fread(fileId,'double');
%       fclose(fileId);

%       nv = 1;
%       Result = reshape(Result,5*nv*np+1,[]); 
%       Xv   = Result(1       :2*nv*np,:);
%       Time = Result(5*nv*np+1       ,:);

%       for ii=1:size(Xv,2)
%         II = secondMoment(Xv(:,ii));
%         e = min(eig(II));
%         inclinationAngle(ii) = atan(II(1,2)/(II(2,2)-e));
%       end
%       %[Ra A L] = reducedVolume(X);
%       %NonDimShear = 1 * (L / 2 / pi)^3 / prams.kappa;

%     end
%   end
% end