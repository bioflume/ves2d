 clear all;clc
% 
% nv = 2;   
% n = 128;   
% angle = 0;
% cc = 1.2; ar = 1/2;
% gg = (1:n)'*2*pi/n; X = [[-cc+cos(gg);ar*sin(gg)] [cc+cos(gg);ar*sin(gg)]];
% %[RA A L] = reducedVolume(X);
% %disp(RA);
% 
% prams.nv = nv;
% prams.T = 12; 
% prams.ts = .0001;
% prams.m = prams.T/prams.ts;
% prams.kappa = .1;
% prams.order = 1;
% prams.vInf = @(X) farFieldVel(X,'extensional',2);
% prams.Case = 'fullSystem';
% 
% options.verbose = 0;
% options.usePlot = 1;   
% options.AxesHandle = gca;
% options.axis = [-2 2 -2 2];      
% options.progressBar = 1;
% options.axisOn = 1;
% 
% options.saveData = 1;
% options.dataStride = 10;%prams.m/5;               
% 
% vc = [.1 1 10];
% for jj=3
%   clear functions global;
%   options.fileName = ['../results/extFlow_h_' num2str(jj) '.mat'];
%   prams.viscCont = vc(jj)*ones(1,nv);
%   [Xfinal status] = Ves2D(X,prams,options);
% end

% Potting the results
n = 128;nv =2;
ww = .43; hh = .35;
W = .45;
H = .45;
COL = [1 1 0 1];
ROW = [1 1 0 0];

figure('Units','inches','position',[3 3 5 4]);
NN = [4 4 6];
for jj=1:3
  options.fileName = ['../results/extFlow_' num2str(jj) '.mat'];
  fileId = fopen(options.fileName,'r');
  Result = fread(fileId,'double');
  fclose(fileId);
  
  Result = reshape(Result,5*nv*n+1,[]); 
  
  Time = Result(end,:);
  X = Result(1:2*nv*n,:);
  sig = Result(2*nv*n+1:3*nv*n,:);
  U = Result(3*nv*n+1:5*nv*n,:);
  
  col = COL(jj);
  row = ROW(jj);
  
  position = [.08+W*col .1+H*row ww hh];
  h = subplot('position',position); 
  hold on;
  for ii =1:NN(jj)
    XX = reshape(X(:,ii),[],2);    
    if(jj<3), Plots = jj;end
    if(jj==3), Plots = [1 2];end
    for kk = Plots
      plot(XX([1:n 1],kk),XX([n+1:2*n n+1],kk),'k','LineWidth',1.5);
    end
  end
  if(jj==1)
    textPos = mean(X(1:n,1)); 
    text(textPos,1,'\nu = 0.1','VerticalAlignment', ...
         'bottom','HorizontalAlignment','center');
  end
  set(gca,'xtick',[]); set(gca,'ytick',[]);
  if(jj==2)
    text(-textPos,1,'\nu = 1','VerticalAlignment', 'bottom','HorizontalAlignment','center');
    set(gca,'xtick',0); grid on
  end
  hold off; box on;
  axis([-2.3 2.3 -1 1]);axis equal
  
  
  if(jj==3) 
    col = COL(jj+1);
    row = ROW(jj+1);
  
    position = [.08+W*col .1+H*row ww hh];
    h = subplot('position',position); 
    hold on;
    for ii =1:NN(jj)
      XX = reshape(X(:,ii),[],2);
      for kk = 1:nv
        plot(XX([1:n 1],kk),XX([n+1:2*n n+1],kk),'k','LineWidth',1.5);
      end
    end
    hold off; 
    axis([-1 1 -.5 .5]);axis equal
    set(gca,'xtick',[]); set(gca,'ytick',[]); box on
  end
end

%% Finding min distance and plotting
col = 0;
row = 1;
position = [.08+W*col .1+H*row ww hh];
h = subplot('position',position); 
hold on;

lineSpec = {'--k','-k','-.k'};
nu = [.1 1 10];
for jj=1:3
  options.fileName = ['../results/extFlow_h_' num2str(jj) '.mat'];
  fileId = fopen(options.fileName,'r');
  Result = fread(fileId,'double');
  fclose(fileId);
  
  Result = reshape(Result,2,[]); 
  r = Result(1,:);
  Time = Result(2,:);
  Ind = Time<10;
  plot(Time(Ind),smooth(r(Ind)),lineSpec{jj},'LineWidth',1.4);
  lambdaCell{jj} = ['\nu = ' num2str(nu(jj))];
end
xlabel('time');
set(gca,'ytick',linspace(0,.2,5));
ylabel('Vesicle''s distance ');
legend(lambdaCell{:}); 
hold off; %axis equal
box on