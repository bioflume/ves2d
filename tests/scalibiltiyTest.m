% P = path; i = find(pwd==filesep); i = i(end);
% subPath = pwd; subPath = [subPath(1:i) 'src' filesep 'base'];
% if isempty(strfind(P, subPath)),addpath(subPath);end
% 
% clear functions
% prams.T = 1e0;
% prams.m = 10;
% prams.viscCont = 2;
% prams.Case = 6;
% 
% prams.kappa = 1;
% prams.mu0 = 1;
% prams.rhoIn = 1;
% prams.rhoOut = 1;
% prams.order = 1;
% prams.ts = prams.T/prams.m;
% prams.Incompressibility = 1;
% prams.vInf = @(X) farFieldVel(X,'shear',1);
% prams.g = [0;-1];
% 
% options.GmresTol = 10^-14;
% options.Lp =  2*pi;
% options.usePlot = 1;
% options.AxesHandle = [];
% options.track = 1;
% options.ProgressBar = 1;
% options.verbose =0;
% 
% n = 64; X = boundary(n);
% Ves2D(X,prams,options);
% 

% profile on
% for lamb = 1
%     for n = 64
%         X = boundary(n);
%         for c = 2:6
%             clear functions
%             profile on
%             prams.viscCont = lamb;
%             prams.Case = c;
%             Ves2D(X,prams,options);
%             profile off;
%             stats.(genvarname(['LNC' num2str(lamb)...
%                 '_' num2str(n) '_' num2str(c)])) = profile('info') ;
%         end
%     end
% end
%
% save('profStats','stats');
% 
clear all
clc
load profStats0
s = fieldnames(stats);
for  ii = 1:size(s)
    f = s{ii};
    p = stats.(f);
    if isfield(p,'FunctionTable')
        
        for jj = 1:size(p.FunctionTable)
            if strcmp(p.FunctionTable(jj).FunctionName,'TimeMatVec')
                Tcalls(ii) = p.FunctionTable(jj).NumCalls;
                
                Ttime(ii) = p.FunctionTable(jj).TotalTime - ...
                    sum([p.FunctionTable(jj).Children.TotalTime]);
                
            end
            if strcmp(p.FunctionTable(jj).FunctionName,'LhsMatVec')
                Lcalls(ii) = p.FunctionTable(jj).NumCalls;
                Ltime(ii) = p.FunctionTable(jj).TotalTime - ...
                    sum([p.FunctionTable(jj).Children.TotalTime]);
            end
        end
    end
end

Tcalls = ceil(reshape(Tcalls,3,3,8)/20); 
Lcalls = ceil(reshape(Lcalls,3,3,8)./Tcalls/20);
Ttime = reshape(Ttime,3,3,8)/20;
Ltime = reshape(Ltime,3,3,8)/20;

clc
fprintf('     \t  Outer loop \t\t | \t Inner loop \n')
fprintf('------------------------------------------------------------\n');
fprintf('    \t %4.0f \t %4.0f \t %4.0f \t |\t %4.0f \t %4.0f \t %4.0f \n',...
    64,128,256,64,128,256);
fprintf('-------------------------------------------------------------');
for ii=1:8
    for jj =1:3
        fprintf('\n %2.0f \t %4.0f \t %4.0f \t %4.0f \t | \t %2.0f \t %4.0f \t %4.0f \t %4.0f',2*jj,...
            Tcalls(jj,:,ii),Lcalls(jj,:,ii));
    end
    fprintf('\n-------------------------------------------------------------');
end
fprintf('\n');










