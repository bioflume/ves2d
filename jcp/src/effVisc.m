clear all; clc

n = 64;
%X = boundary(n,'reducedArea',.9);

prams.T = 10;                            
prams.ts = 1e-1;                  
prams.kappa = 1;
prams.order = 1;                        
prams.vInf = @(X) farFieldVel(X,'shear',2);
prams.volRatio = 1;

options.usePlot = 1;                    
options.axis = [-5 5 -5 5];
options.axisOn = 1;
options.track = 0;                      
options.progressBar = 0;
options.verbose = 0;
options.AxesHandle = [];%subplot(1,2,1);
options.showError = 1;

% vc = 10.^(.6:.1:1);
% for ii = 1:length(vc)
%   disp(vc(ii));
%   clear functions global;
%   prams.viscCont = vc(ii);                     
% 
%   [Xfinal status avgS(:,:,ii)] = Ves2D(X,prams,options);
% end

%% Post processing
load ../results/effVisc6.mat

muA = @(nu) 2.5-(23*nu+32)*.6/16/pi;

sigma21 = squeeze(AVGS(3,:,:));
VC = VC(1:size(AVGS,3));

efV = mean(sigma21);
xx = [10.^(-1:.05:1)];
yy = spline(VC,efV,xx);

[bifPnt ind] = min(yy);ind=ind-5;

subplot(3,1,1);
semilogx(xx,10*smooth(yy)','-k','linewidth',1.5)
hold; semilogx(xx(1:ind),muA(xx(1: ind)),'-.k','linewidth',1.5);
axis([.1 10 1 2.5])
legend('[\mu]','[\mu]_a')
ylabel('[\mu]')
text(.2,2.2,'\Delta = .6')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
load ../results/effVisc.mat

muA = @(nu) 2.5-(23*nu+32)*.75/16/pi;

efV = mean(sigma21);
xx = [10.^(-1:.05:1)];
yy = spline(VC,efV,xx);

[bifPnt ind] = min(yy);ind=ind-5;

subplot(3,1,2);
semilogx(xx,10*smooth(yy)','-k','linewidth',1.5)
hold; semilogx(xx(1:ind),muA(xx(1: ind)),'-.k','linewidth',1.5);
axis([.1 10 1 2.5])
ylabel('[\mu]')
text(.2,2.2,'\Delta = .75')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load ../results/effVisc9.mat

muA = @(nu) 2.5-(23*nu+32)*.9/16/pi;

sigma21 = squeeze(AVGS(3,:,:));
VC = VC(1:size(AVGS,3));

efV = mean(sigma21);
xx = [10.^(-1:.05:1)];
yy = spline(VC,efV,xx);

[bifPnt ind] = min(yy);ind=ind-5;

subplot(3,1,3);
semilogx(xx,10*smooth(yy)','-k','linewidth',1.5)
hold; semilogx(xx(1:ind),muA(xx(1: ind)),'-.k','linewidth',1.5);
axis([.1 10 1 2.5])
ylabel('[\mu]')
xlabel('\nu')
text(.2,2.2,'\Delta = .9')

close all; drawnow
load ../results/effVisc.mat
Delta = .75;
muA = @(nu) 2.5-(23*nu+32)*Delta/16/pi;

efV = 10*mean(squeeze(AVGS(3,:,:)));
xx = [10.^(-1:.05:1)];
yy = spline(VC,efV,xx);

[bifPnt ind] = min(yy);

subplot('position',[.1 .55 .8 .3]);
semilogx(xx,smooth(yy)','-k','linewidth',1.5);
hold; semilogx(xx(1:ind),muA(xx(1: ind)),'-.k','linewidth',1.5);
axis([.1 10 1 2.5]);
%set(gca,'xtick',[]);
ylabel('[\mu]')
text(.25,2.2,['\Delta = ' num2str(Delta)],'VerticalAlignment','bottom','HorizontalAlignment', 'center');
legend('[\mu]','[\mu]_a')

load ../results/effVisc9.mat
Delta = .9;
muA = @(nu) 2.5-(23*nu+32)*Delta/16/pi;

efV = 10*mean(squeeze(AVGS(3,:,:)));

yy = spline(VC,efV,xx);

[bifPnt ind] = min(yy);

subplot('position',[.1 .2 .8 .3]);
semilogx(xx,smooth(yy)','-k','linewidth',1.5);
hold; semilogx(xx(1:ind),muA(xx(1: ind)),'-.k','linewidth',1.5);
axis([.1 10 1 2.5]);
text(.25,2.2,['\Delta = ' num2str(Delta)],'VerticalAlignment','bottom','HorizontalAlignment', 'center');

xlabel('\nu')
ylabel('[\mu]')

