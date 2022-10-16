toclear = [];
toclear = who;
ind = ~strcmp(toclear,'testCase');
clear(toclear{ind});

addpath ../src/

fileSpec.prefix  = '';
fileSpec.perturb = [];
fileSpec.volFrac = [];
fileSpec.reducedArea = [];
fileSpec.viscCont = [];
fileSpec.T       = 200;
fileSpec.ts      = 1e-3;
fileSpec.n       = 64;

if(~any(strcmp(who,'testCase')) || isempty(testCase)), testCase = 5;end

switch testCase
 case 1
  fileSpec.prefix = './results/couetteFlow';
  fileSpec.volFrac = 3;
  [Xv Time C effVisc] = postProcess(fileSpec,'TimePlot','Trajectory', 'ErrorPlot', 'EffectiveViscosity');

 case 2
  fileSpec.prefix = './results/uniformCouetteFlowCPU';
  fileSpec.perturb = .01;
  fileSpec.volFrac = 19;
  fileSpec.T = 100;
  [Xv Time C effVisc] = postProcess(fileSpec,'TimePlot','Trajectory', 'ErrorPlot', 'EffectiveViscosity');
  
 case 3
  fileSpec.prefix = './results/uniformCouetteFlowCPU_tiered';
  fileSpec.perturb = .04;
  fileSpec.volFrac = 31;
  fileSpec.T = 100;
  [Xv Time C effVisc] = postProcess(fileSpec,'TimePlot','Trajectory', 'ErrorPlot', 'EffectiveViscosity');
  
 case 4
  fileSpec.prefix = './results/uniformCouetteFlowCPU_radial';
  fileSpec.perturb = .15;
  fileSpec.volFrac = 19;
  fileSpec.T = 100;
  [Xv Time C effVisc] = postProcess(fileSpec,'TimePlot','Trajectory', 'ErrorPlot', 'EffectiveViscosity');

 case 5 
  fileSpec.prefix  = './couetteFlow';
  fileSpec.perturb = [];
  fileSpec.viscCont = [];
  fileSpec.T       = 200;
  fileSpec.ts      = 1e-3;
  fileSpec.n       = 64;

  volFrac = [2 3 4 5 6 7 8 9 10 20 30 40];
  reducedArea = [.65 .8];
  f1 = figure;
  f2 = figure;
  
  for ii = 1:length(volFrac)
    for jj = 1:length(reducedArea)
      fileSpec.volFrac = volFrac(ii);
      fileSpec.reducedArea = reducedArea(jj);
  
      figure(f1);
      h1 = subplot(length(reducedArea),length(volFrac),length(volFrac)*(jj-1)+ii);
      
      figure(f2);
      h2 = subplot(length(reducedArea),length(volFrac),length(volFrac)*(jj-1)+ii);
  
      postProcess(fileSpec,'Trajectory', h1,'EffectiveViscosity',h2);
      title(h1,num2str(volFrac(ii)));
      title(h2,num2str(volFrac(ii)));
      %axis(h1,[0 500 10 20]);
      %axis(h2,[0 500 1 1.3]);
    end
  end

 case 6 
  fileSpec.prefix  = './couetteFlow';
  fileSpec.perturb = [];
  fileSpec.viscCont = [];
  fileSpec.T       = 200;
  fileSpec.ts      = 1e-3;
  fileSpec.n       = 64;

  volFrac = [6 7 8 9 10];
  reducedArea = .8;
  f1 = figure; h1 = gca;
  f2 = figure;

  for ii = 1:length(volFrac)
    for jj = 1:length(reducedArea)
      fileSpec.volFrac = volFrac(ii);
      fileSpec.reducedArea = reducedArea(jj);
  
      figure(f2);
      h2 = subplot(length(reducedArea),length(volFrac),length(volFrac)*(jj-1)+ii);
      postProcess(fileSpec,'TimePlot',h1,'Trajectory', h2);
      title(h2,num2str(volFrac(ii)));
    end
  end

 case 7 %Two colored plot of 6% volume fraction
  fileSpec.prefix  = './couetteFlow';
  fileSpec.perturb = [];
  fileSpec.viscCont = [];
  fileSpec.T       = 200;
  fileSpec.ts      = 1e-3;
  fileSpec.n       = 64;
  fileSpec.volFrac = 6;
  fileSpec.reducedArea = .65;
  
  nv = 18;
  idx1 = [1 4 7 8 10 12 14 18];
  idx2 = [2 3 5 6 9 11 13 15 16 17];

  [Xv Time C effVisc options prams] = postProcess(fileSpec);
  domain = fixedBound(prams.M,prams.bd,1);
  
  nameformat = ceil(log10(length(Time)));
  nameformat = ['%0' num2str(nameformat) 'u'];

  for jj=1:length(Time)
    clf;
    hold on;
    XX = domain(1).X; XX = [XX;XX(1,:)];
    plot(XX(:,1),XX(:,2),'k','LineWidth',2);

    XX = domain(2).X; XX = [XX;XX(1,:)];
    plot(XX(:,1),XX(:,2),'k','LineWidth',2);
    
    
    X = reshape(Xv(:,jj),[],2*nv);
    X = [X; X(1,:)];
    X = reshape(X,[],nv);

    X1 = X(:,idx1);
    X2 = X(:,idx2);
    for ii=1:size(X1,2)
      XX = reshape(X1(:,ii),[],2);
      plot(XX(:,1),XX(:,2),'r','Line Width',1.2);
    end

    for ii=1:size(X2,2)
      XX = reshape(X2(:,ii),[],2);
      plot(XX(:,1),XX(:,2),'b','LineWidth',1.2);
    end
    
    axis equal
    axis off
    title(num2str(Time(jj),'%02.2f'));
    drawnow;

    %saveas(gcf,['Couette_6_' num2str(jj,nameformat)],'jpg');

    pause(.1);
  end

 case 8
  fileSpec.prefix  = './couetteFlow';
  fileSpec.perturb = [];
  fileSpec.T       = 200;
  fileSpec.ts      = 5e-5;
  fileSpec.n       = 128;
  fileSpec.volFrac = 2;

  viscCont = [2 4];
  reducedArea = [.65 .8];
  f1 = figure; h1 = gca;
  f2 = figure;
  for ii = 1:length(viscCont)
    for jj = 1:length(reducedArea)
      fileSpec.viscCont = viscCont(ii);      
      fileSpec.reducedArea = reducedArea(jj);
  
      figure(f2);
      h2 = subplot(length(reducedArea),length(viscCont),length(viscCont)*(jj-1)+ii);
      postProcess(fileSpec,'TimePlot',h1,'Trajectory', h2);
      title(h2,num2str(viscCont(ii)));
      %postProcess(fileSpec,'Trajectory', h2);
      ylim([10 20]);
    end
  end
 case 9
  fileSpec.prefix  = './couetteFlow';
  fileSpec.perturb = [];
  fileSpec.T       = 200;
  fileSpec.ts      = 1e-3;
  fileSpec.n       = 64;
  fileSpec.viscCont = [];
  fileSpec.reducedArea = .65;

  VF = [20];
  for vf = VF
    fileSpec.volFrac=vf;
    %postProcess(fileSpec,'TimePlot','Trajectory','EffectiveViscosity');
    postProcess(fileSpec,'Trajectory','EffectiveViscosity');
    pause;
  end
end
