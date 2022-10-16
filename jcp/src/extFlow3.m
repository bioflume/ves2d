% clear all;clc
% 
% n = 128;nv =2;
% options.fileName = ['../results/extFlow_3.mat'];
% fileId = fopen(options.fileName,'r');
% Result = fread(fileId,'double');
% fclose(fileId);
% Result = reshape(Result,5*nv*n+1,[]); 
% Xf = reshape(Result(1:2*n*nv,end),[],2);
%  
% XX = reshape(Xf(2*n+1:end),[],2);
% 
% [r ind] = min(XX(1:n/2,1));
%  
% nn = 32; gg = (1:nn)'*2*pi/nn; 
% capCCW = [cos(gg) sin(gg)];
% capCW  = [cos(gg) -sin(gg)];
% 
% bigCut = ind-2;
% smallCut = ind+2;
% 
% bCap = XX(bigCut,1)*capCW(1:nn/4,:); 
% bCap(:,2) = bCap(:,2)+XX(bigCut,2);
% 
% sCap = XX(smallCut,1)*capCCW(1:nn/4,:);
% sCap(:,2) = sCap(:,2)+XX(smallCut,2);
% 
% XB = [XX([end 1:bigCut],:);bCap];
% XS = [XX(n/2:-1:smallCut,:);sCap];
% 
% XB = [XB;[-XB(end-1:-1:1,1) XB(end-1:-1:1,2)]];
% XB = [XB;[XB(end-1:-1:2,1) -XB(end-1:-1:2,2)]];
% 
% XS = [XS;[-XS(end-1:-1:1,1) XS(end-1:-1:1,2)]];
% XS = [XS;[XS(end-1:-1:2,1) -XS(end-1:-1:2,2)]];
% 
% XB = interpft(XB,512);XB = XB(1:4:end,:);
% XS = interpft(XS,512);XS = XS(1:4:end,:);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prams.nv = 1;
% prams.kappa = .1;
% prams.order = 1;
% prams.vInf = @(X) farFieldVel(X,'extensional',0);
% prams.Case = 'fullSystem';
% 
% options.verbose = 0;
% options.usePlot = 1;   
% options.AxesHandle = gca;
% options.progressBar = 1;
% options.axisOn = 1;
% 
% clear functions global;
% prams.T = 2;
% prams.ts = .02;
% prams.m = prams.T/prams.ts;
% options.axis = [-2 2 -2 2];      
% prams.viscCont = 10;
% [XBf status] = Ves2D(XB(:),prams,options);
% 
% clear functions global;
% prams.T = .1;
% prams.ts = .005;
% prams.m = prams.T/prams.ts;
% options.axis = [-.5 .5 -.7 .7];      
% prams.viscCont = .1;
% [XSf status] = Ves2D(XS(:),prams,options);

XBf = reshape(XBf,[],2);
XSf = reshape(XSf,[],2);

subplot(3,1,1)
hold on;
plot(XX(:,1),XX(:,2),'-k');
plot(-XX(:,1),XX(:,2),'-k');
hold off; axis([-2 2 -1 1]);
title(['bending energy = ' num2str(2*Energy(XX(:)))]);

subplot(3,1,2)
hold on;
plot(XB(:,1),XB(:,2),'-k');
plot(-XS(:,1),XS(:,2),'-k');
hold off; axis([-2 2 -1 1]);
title(['bending energy = ' num2str(Energy(XB(:))) ',  ' num2str(Energy(XS(:)))]);

subplot(3,1,3)
hold on;
plot(XBf(:,1),XBf(:,2),'-k');
plot(XSf(:,1),XSf(:,2),'-k');
hold off; axis([-2 2 -1 1]);
title(['bending energy = ' num2str(Energy(XBf(:))) ',  ' num2str(Energy(XSf(:)))]);