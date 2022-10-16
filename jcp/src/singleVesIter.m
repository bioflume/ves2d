addpath ../../src;

clear all;% clc

prams.kappa = 1;  
prams.order = 1; 
prams.Case = 'schurComp';

options.usePlot = 0;                    
options.progressBar = 0;
options.verbose = 0;
options.GmresTol = 10^-8;    

format compact;

TT = 2;
vc = [.01 .1 1 10 100];
shear = [0 1 10];

TITER = zeros(4,length(shear),length(vc),length(TT));
LITER = TITER;

for ii = 1:length(TT)
  disp('====================================================');
  disp(['                 TT = ' num2str(TT(ii))]);
  for jj = 1:length(vc)
    for kk = 1:length(shear)
      disp('====================================================');
      for a = 5:8
        prams.T = TT(ii)/(shear(kk) + (shear(kk)==0));
        prams.vInf = @(x) farFieldVel(x,'shear',shear(kk));
        prams.viscCont = vc(jj);
        prams.n = 2^a;  
        prams.m = 100;
        prams.ts = prams.T/prams.m;
        
        X = boundary(prams.n);
        clear functions global;
        
        [Xfinal status] = Ves2D(X,prams,options); 
        [xx ss ff tIter] = timeStepper();
        [xx ss ff lIter] = TimeMatVec();
        
        TITER(a-4,kk,jj,ii) = ceil(mean(tIter));
        LITER(a-4,kk,jj,ii) = ceil(mean(lIter));

        disp([vc(jj) shear(kk) prams.n prams.ts status.flag ceil(mean(tIter)) ceil(mean(lIter))]);
      end
    end
  end
end
save('../iterationsNoPrecond.mat','TITER','LITER');

for ii=1:length(TT)
  for a=5:8
    fprintf('%3.0f &',2^a)
    for jj=1:length(vc)
      for kk=1:length(shear)
        fprintf('%3.0f &',TITER(a-4,kk,jj,ii));
      end
    end
    fprintf('\n');
  end
end
disp('=====================================================')
for ii=1:length(TT)
  for a=5:8
    fprintf('%3.0f &',2^a)
    for jj=1:length(vc)
      for kk=1:length(shear)
        fprintf('%3.f &',LITER(a-4,kk,jj,ii));
      end
    end
    fprintf('\n');
  end
end

% clear all; clc
% 
% prams.m = 20;
% prams.kappa = 1;  
% prams.order = 1; 
% prams.Case = 'schurComp';
% 
% options.usePlot = 0;                    
% options.progressBar = 0;
% options.verbose = 0;
% options.GmresTol = 10^-8;    
% load('../results/stableStepSize_Ra75_T6_shear_vc_n_m')
% 
% profile on
% for ii=1:3
%   for jj= 1:20
%     shear = sss(jj,1,ii);
%     vc = sss(jj,2,ii);
%     N = sss(jj,3,ii);
%     M = sss(jj,4,ii); 
%     
%     prams.T = 6/M;
%     prams.vInf = @(X) farFieldVel(X,'shear',shear);
%     prams.viscCont = vc;
%     prams.n = N;  
%     
%     X = boundary(prams.n);
%     clear functions global;
% 
%     %dbstop in PrecondTM at 32
%     profile clear
%     [Xfinal status] = Ves2D(X,prams,options); 
%     
%     stats.(genvarname(['ij' num2str(ii) '_' num2str(jj)])) = profile('info');
%     disp([ii jj status.flag]);
%   end
% end
% profile off;
% fileName = '/home/abtin/works/Ves2D/jcp/results/singleVesIterProfile_8_12';   
% save(fileName,'stats');
% 
% % load(fileName); 
% s = fieldnames(stats);
% 
% for  ii = 1:size(s)
%   f = s{ii};
%   p = stats.(f);
%   if isfield(p,'FunctionTable')
%     
%     for jj = 1:size(p.FunctionTable)
%       if strcmp(p.FunctionTable(jj).FunctionName,'TimeMatVec')
%         Tcalls(ii) = p.FunctionTable(jj).NumCalls;   
%         %Ttime(ii) = p.FunctionTable(jj).TotalTime - ...
%         % sum([p.FunctionTable(jj).Children.TotalTime]);
%       end
%       if strcmp(p.FunctionTable(jj).FunctionName,'LhsMatVec')
%         Lcalls(ii) = p.FunctionTable(jj).NumCalls;
%         %Ltime(ii) = p.FunctionTable(jj).TotalTime - ...
%         % sum([p.FunctionTable(jj).Children.TotalTime]);
%       end
%     end
%   end
% end
% Lcalls = Lcalls./(Tcalls+1); %one Lcall is outside Tcall
% 
% Tcalls = floor(reshape(Tcalls,20,3)/20-4); 
% Lcalls = floor(reshape(Lcalls,20,3));
% 
% % %Ttime = reshape(Ttime,20,3)/20;
% % %Ltime = reshape(Ltime,20,3)/20;
% % 
% fileName = '/home/abtin/works/Ves2D/jcp/results/singleVesIter.txt';   
% fileId = fopen(fileName,'a');
% if(fileId== -1);return;end
% 
% midLine = [repmat('-',1,75) '\n'];
% fprintf(fileId,midLine);  
% printString = [' %-3.0f & %-3.0f & %-3.0f&'];
% printString = [repmat(printString,1,5) '\n'];
% indSet = 1:4:20;
% 
% for ii = 1:2
%   for jj = 0:3
%     switch ii
%      case 1
%       MAT = Tcalls;
%      case 2
%       MAT = Lcalls;
%     end
%     printVec = MAT(indSet+jj,:);
%     printVec = printVec'; printVec=printVec(:)';
%     fprintf(fileId,printString,printVec);  
%   end
%   fprintf(fileId,midLine);  
% end
% 
% fclose(fileId);




