clear all;clc

n = 64;
X0 = boundary(n);

prams.T = 6;
prams.kappa = 1;
prams.order = 1;   
prams.Incompressibility = 1;

options.usePlot = 0;        
options.track = 0;          
options.axis = [-2 2 -2 2];

options.findSteadyState = 1;
options.saveFig = 0;                   
options.velProfile = 0 ;               
options.progressBar = 0;
options.saveData = 0;
options.verbose = 0;

shear = [1 2 3 4 5 6 7 8 9 10];
nu = [.01 .1 .2:.2:5];
M = 2048;

for jj = 8:length(shear) 
  prams.vInf = @(X) farFieldVel(X,'shear',shear(jj));
 
  for ii = 1:length(nu)
    prams.m = M;
    prams.viscCont = nu(ii);
    clear functions global;
    global state;
    [X status] = Ves2D(X0,prams,options);
    disp([jj ii status.flag]); viewer(X);
    T(:,ii,jj) = state.timeSteps(2:8:end)';
    B(:,ii,jj) = state.beta(1:8:end)';
  end 
end

% nu = [.01 .1 .2:.2:5];
% shear = [1 2 3 4 5 6 7 8 9 10];
% load '../results/incAngle.mat';
% M = size(B,1);

% for ii=1:size(B,3)
%   DB(2:M,:,ii) = (B(2:M,:,ii)-B(1:M-1,:,ii))./(T(2:M,:,ii)-T(1:M-1,:,ii));
%   DB(1,:,ii) = DB(2,:,ii);
% end

% cond = abs(DB./B) <5e-2;
% for ii=1:size(B,3)
%   for jj =1:size(B,2)
%     betaSteady(ii,jj) = mean(B(cond(:,jj,ii),jj,ii));
%   end
% end

% lineSpec = {'-+k','-ok','-sk','-dk','-hk'};
% subplot(1,2,2);hold on;
% for ii=1:size(betaSteady,1)/2
%   plot(nu,betaSteady(2*ii,:),lineSpec{ii},'LineWidth',1.4); 
% end
% hold off;
% sCell = cellstr(num2str(shear(2:2:end)'));
% for ii=1:size(sCell)
%   sCell{ii} = ['\chi = ' sCell{ii}];
% end
% legend(sCell{:}); 
% axis([min(nu) max(nu(1:3:end)) 0 pi/5]);
% set(gca,'xtick',nu(1:3:end));
% set(gca,'ytick',(0:2)*pi/16);
% set(gca,'yminortick','on');
% set(gca,'YTickLabel',{'0';'\pi/16';'\pi/8'});
% xlabel('\nu');
% ylabel('Tank treading inclination angle');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ind = find(nu<4,1,'last');
% nuCell = cellstr(num2str(nu(1:end)'));
% for ii=1:size(nuCell)
%   nuCell{ii} = ['\nu = ' nuCell{ii}];
% end
  
% subplot(1,2,1);hold on;
% for ii=1:3:ind
%   plot(shear,betaSteady(:,ii),'-k','LineWidth',2); 
%   text(shear(3),betaSteady(3,ii),nuCell{ii},'VerticalAlignment','bottom');
% end
% hold off;

% axis([min(shear) max(shear) 0 pi/5]);
% set(gca,'ytick',(0:2)*pi/16);
% set(gca,'yminortick','on');
% set(gca,'YTickLabel',{'0';'\pi/16';'\pi/8'});
% xlabel('\chi');
% ylabel('Tank treading inclination angle');
