clear; clc;
addpath ./necessaryMatFiles/

%load /workspace/gokberk/only1timeStepPCAData.mat
%load /workspace/gokberk/only1timeStepTrainPCAData.mat
%load /workspace/gokberk/only1timeStepMoreModesPCAData.mat 
load /workspace/gokberk/n256AllDtTimesKappasPCAData.mat
load ./necessaryMatFiles/pcaCoeffsBasis1step.mat
% PREDICT NEW SHAPE

nComp = 16; % use 16 principal components
compStart = 17;
compEnd = compStart-1+16;

netType = 'resNet'; % fcNet, resNet

fileName = ['relaxNetAllDtTimesKbs_' netType ...
    '_nModes' num2str(compStart) 'to' num2str(compEnd)];

nves = nInstances;

nTrain = floor(0.95*nves);
nTest = nves-nTrain;
Xtrain = zeros(nComp+1,1,1,nTrain);
Ytrain = zeros(nTrain,nComp);
Xtest = zeros(nComp+1,1,1,nTest);
Ytest = zeros(nTest,nComp);


if 0 % randomly choose idcs or load already randomly chosen one
  testIdcs = randperm(nves,nTest)';
  trainIdcs = zeros(nTrain,1);

  idx = 1;
  for k = 1 : nves
    if ~any(k==testIdcs)
      trainIdcs(idx) = k;
      idx = idx+1;
    end
  end
  save setOfTrainTestIdcsAllDtTimesKs testIdcs trainIdcs
else
  load setOfTrainTestIdcsAllDtTimesKs
end

% PREPARE INPUT AND OUTPUT

for i = 1 : nTrain
  Xtrain(1:nComp,1,1,i) = reshape(XoldCoeffs(trainIdcs(i),compStart:compEnd),...
      nComp,1,1);
  Xtrain(nComp+1,1,1,i) = DtTimesKappa(trainIdcs(i));
  Ytrain(i,:) = XnewCoeffs(trainIdcs(i),compStart:compEnd);    
end
for i = 1:nTest
  Xtest(1:nComp,1,1,i) = reshape(XoldCoeffs(testIdcs(i),compStart:compEnd),...
      nComp,1,1);
  Xtest(nComp+1,1,1,i) = DtTimesKappa(testIdcs(i));
  Ytest(i,:) = XnewCoeffs(testIdcs(i),compStart:compEnd);
end


% Normalize or scale inputs to a certain range
% Compute mean and stdev of channels
chan1vals = Xtrain(1:16,1,1,:); 
muChan1 = mean(chan1vals(:)); sdevChan1 = std(chan1vals(:));   
% Normalize
scale = 1; offset = 0;

Xtrain(1:16,1,1,:) = scale*(Xtrain(1:16,1,1,:)-muChan1)/sdevChan1+offset;
Ytrain = scale*(Ytrain-muChan1)/sdevChan1+offset;

Xtest(1:16,1,1,:) = scale*(Xtest(1:16,1,1,:)-muChan1)/sdevChan1+offset;
Ytest = scale*(Ytest-muChan1)/sdevChan1+offset;

Xtrain(17,1,1,:) = 1000 * Xtrain(17,1,1,:); % scale it to match PCA coeffs 
Xtest(17,1,1,:) = 1000 * Xtest(17,1,1,:); % scale it to match PCA coeffs 
% Normalize
scale = 1; offset = 0;

Xtrain(1:16,1,1,:) = scale*(Xtrain(1:16,1,1,:)-muChan1)/sdevChan1+offset;
Xtest(1:16,1,1,:) = scale*(Xtest(1:16,1,1,:)-muChan1)/sdevChan1+offset;

% LOAD NETWORK
[lgraph,numParams] = chooseNetwork(netType,evects,colMeans',nComp);

disp(['There are ' num2str(numParams) ' # of parameters.'])
disp(['There are ' num2str(nTrain) ' # of training instances.'])

miniBatchSize  = 256;
validationFrequency = floor(numel(Ytrain)/miniBatchSize);
options = trainingOptions('sgdm', ...
    'Momentum',0.9,...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',25, ...
    'InitialLearnRate',5e-4, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.2, ...
    'LearnRateDropPeriod',10, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{Xtest,Ytest},...
    'ValidationFrequency',validationFrequency,...
    'Verbose',true,...
    'VerboseFrequency',5);

% TRAIN NETWORK
net = trainNetwork(Xtrain,Ytrain,lgraph,options);


% PREDICT TEST DATA
Ypred = predict(net,Xtest); % compare with Ytest
YtrainPred = predict(net,Xtrain);

YOPred = (Ypred-offset)*sdevChan1/scale+muChan1;
YOTrainPred = (YtrainPred-offset)*sdevChan1/scale+muChan1;
Ytest = (Ytest-offset)*sdevChan1/scale+muChan1;
Ytrain = (Ytrain-offset)*sdevChan1/scale+muChan1;

err = sqrt(1/nComp*sum((YOPred-Ytest).^2,2))./...
    sqrt(1/nComp*sum(Ytest.^2,2));
disp(['Max error is ' num2str(max(err))])
disp(['Min error is ' num2str(min(err))])
disp(['Mean error is ' num2str(mean(err))])

netFileName = [fileName '.mat'];
save(netFileName,'net','muChan1','sdevChan1',...
'scale','offset','nComp')

%save(['/workspace/gokberk/' fileName '_w1step.mat'],'YOPred','Ytest',...
%    'Xtrain','Xtest','YOTrainPred','Ytrain')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lgraph,numParams] = chooseNetwork(netType,M,mu,nComp)

if strcmp(netType,'resNet') % 335,296
    
filterSizes = [];
filterSizes(1,:) = [3 1 1];
filterSizes(2,:) = [3 1 2];
filterSizes(3,:) = [3 1 2];
filterSizes(4,:) = [3 1 4];
filterSizes(5,:) = [3 1 4];
filterSizes(6,:) = [3 1 8];
filterSizes(7,:) = [3 1 8];
filterSizes(8,:) = [3 1 8];
filterSizes(9,:) = [3 1 16];
filterSizes(10,:) = [3 1 16];
filterSizes(11,:) = [3 1 16];
filterSizes(12,:) = [3 1 32];
filterSizes(13,:) = [3 1 32];
filterSizes(14,:) = [3 1 32];
filterSizes(15,:) = [3 1 64];
filterSizes(16,:) = [3 1 64];
filterSizes(17,:) = [3 1 64];
filterSizes(18,:) = [3 1 64];
filterSizes(19,:) = [3 1 96];
filterSizes(20,:) = [3 1 96];
filterSizes(21,:) = [3 1 96];
filterSizes(22,:) = [3 1 96];
filterSizes(23,:) = [3 1 128];
filterSizes(24,:) = [3 1 128];
filterSizes(25,:) = [3 1 128];

layers = [imageInputLayer([17 1 1],'Name','input','Normalization','none')

    convolution2dLayer(filterSizes(1,1:2),filterSizes(1,3),'Name','conv1','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm1')
    leakyReluLayer(0.1,'Name','relu1')
    
    additionLayer(2,'Name','add1')
    
    convolution2dLayer(filterSizes(2,1:2),filterSizes(2,3),'Name','conv2','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm2')
    leakyReluLayer(0.1,'Name','relu2')
    
    convolution2dLayer(filterSizes(3,1:2),filterSizes(3,3),'Name','conv3','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm3')
    leakyReluLayer(0.1,'Name','relu3')
    
    additionLayer(2,'Name','add3')
    
    convolution2dLayer(filterSizes(4,1:2),filterSizes(4,3),'Name','conv4','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm4')
    leakyReluLayer(0.1,'Name','relu4')
    
    convolution2dLayer(filterSizes(5,1:2),filterSizes(5,3),'Name','conv5','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm5')
    leakyReluLayer(0.1,'Name','relu5')
    
    additionLayer(2,'Name','add5')
    
    convolution2dLayer(filterSizes(6,1:2),filterSizes(6,3),'Name','conv6','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm6')
    leakyReluLayer(0.1,'Name','relu6')
    
    
    convolution2dLayer(filterSizes(7,1:2),filterSizes(7,3),'Name','conv7','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm7')
    leakyReluLayer(0.1,'Name','relu7')
    
    additionLayer(2,'Name','add7')
    
    convolution2dLayer(filterSizes(8,1:2),filterSizes(8,3),'Name','conv8','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm8')
    leakyReluLayer(0.1,'Name','relu8')
    
    additionLayer(2,'Name','add8')
    
    convolution2dLayer(filterSizes(9,1:2),filterSizes(9,3),'Name','conv9','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm9')
    leakyReluLayer(0.1,'Name','relu9')
    
    convolution2dLayer(filterSizes(10,1:2),filterSizes(10,3),'Name','conv10','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm10')
    leakyReluLayer(0.1,'Name','relu10')
    
    additionLayer(2,'Name','add10')
    
    convolution2dLayer(filterSizes(11,1:2),filterSizes(11,3),'Name','conv11','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm11')
    leakyReluLayer(0.1,'Name','relu11')
    
    additionLayer(2,'Name','add11')
    
    convolution2dLayer(filterSizes(12,1:2),filterSizes(12,3),'Name','conv12','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm12')
    leakyReluLayer(0.1,'Name','relu12')
    
    convolution2dLayer(filterSizes(13,1:2),filterSizes(13,3),'Name','conv13','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm13')
    leakyReluLayer(0.1,'Name','relu13')
    
    additionLayer(2,'Name','add13')
    
    convolution2dLayer(filterSizes(14,1:2),filterSizes(14,3),'Name','conv14','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm14')
    leakyReluLayer(0.1,'Name','relu14')
    
    additionLayer(2,'Name','add14')
    
    convolution2dLayer(filterSizes(15,1:2),filterSizes(15,3),'Name','conv15','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm15')
    leakyReluLayer(0.1,'Name','relu15')
    
    convolution2dLayer(filterSizes(16,1:2),filterSizes(16,3),'Name','conv16','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm16')
    leakyReluLayer(0.1,'Name','relu16')
    
    additionLayer(2,'Name','add16')
    
    convolution2dLayer(filterSizes(17,1:2),filterSizes(17,3),'Name','conv17','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm17')
    leakyReluLayer(0.1,'Name','relu17')
    
    additionLayer(2,'Name','add17')
    
    convolution2dLayer(filterSizes(18,1:2),filterSizes(18,3),'Name','conv18','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm18')
    leakyReluLayer(0.1,'Name','relu18')
    
    additionLayer(2,'Name','add18')
    
    convolution2dLayer(filterSizes(19,1:2),filterSizes(19,3),'Name','conv19')
    batchNormalizationLayer('Name','bnorm19')
    leakyReluLayer(0.1,'Name','relu19')
    
    convolution2dLayer(filterSizes(20,1:2),filterSizes(20,3),'Name','conv20')
    batchNormalizationLayer('Name','bnorm20')
    leakyReluLayer(0.1,'Name','relu20')
    
    convolution2dLayer(filterSizes(21,1:2),filterSizes(21,3),'Name','conv21')
    batchNormalizationLayer('Name','bnorm21')
    leakyReluLayer(0.1,'Name','relu21')
    
    convolution2dLayer(filterSizes(22,1:2),filterSizes(22,3),'Name','conv22')
    batchNormalizationLayer('Name','bnorm22')
    leakyReluLayer(0.1,'Name','relu22')
    
    convolution2dLayer(filterSizes(23,1:2),filterSizes(23,3),'Name','conv23')
    batchNormalizationLayer('Name','bnorm23')
    leakyReluLayer(0.1,'Name','relu23')
    
    convolution2dLayer(filterSizes(24,1:2),filterSizes(24,3),'Name','conv24')
    batchNormalizationLayer('Name','bnorm24')
    leakyReluLayer(0.1,'Name','relu24')
    
    convolution2dLayer(filterSizes(25,1:2),filterSizes(25,3),'Name','conv25')
    batchNormalizationLayer('Name','bnorm25')
    leakyReluLayer(0.1,'Name','relu25')
    
    fullyConnectedLayer(96,'Name','fc26')
    batchNormalizationLayer('Name','bnorm26')
    leakyReluLayer(0.1,'Name','relu26')
    
    fullyConnectedLayer(64,'Name','fc27')
    batchNormalizationLayer('Name','bnorm27')
    leakyReluLayer(0.1,'Name','relu27')
    
    fullyConnectedLayer(32,'Name','fc28')
    batchNormalizationLayer('Name','bnorm28')
    leakyReluLayer(0.1,'Name','relu28')
    
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(16,'Name','fcLast')
    regressionLayer('Name','reg')
];
% skipConv1 = convolution2dLayer([3 1],3,'Padding','Same','Name','skipConv1');

lgraph = layerGraph(layers); 
lgraph = connectLayers(lgraph,'input','add1/in2');
lgraph = connectLayers(lgraph,'relu2','add3/in2');
lgraph = connectLayers(lgraph,'relu4','add5/in2');
lgraph = connectLayers(lgraph,'relu6','add7/in2');
lgraph = connectLayers(lgraph,'add7','add8/in2');
lgraph = connectLayers(lgraph,'relu9','add10/in2');
lgraph = connectLayers(lgraph,'add10','add11/in2');
lgraph = connectLayers(lgraph,'relu12','add13/in2');
lgraph = connectLayers(lgraph,'add13','add14/in2');
lgraph = connectLayers(lgraph,'relu15','add16/in2');
lgraph = connectLayers(lgraph,'add16','add17/in2');
lgraph = connectLayers(lgraph,'add17','add18/in2');

nConvLayers = numel(filterSizes(:,1));
numParams = 0;
for k = 1 : nConvLayers
  if k == 1
    numParams = numParams + 2*prod(filterSizes(k,:));    
  else
    numParams = numParams + filterSizes(k-1,3)*prod(filterSizes(k,:));  
  end  
end

numParams = numParams + 128*3*96 + 96*64 + 64*32 + 16*32;


elseif strcmp(netType,'fcNet')
numParams = 17*32+32*64*2+64*128*2+128*256*2+256*512*2+16*32;

layers = [imageInputLayer([17 1 1],'Name','input','Normalization','none')

    fullyConnectedLayer(32,'Name','fc1')
    batchNormalizationLayer('Name','bnorm1')
    leakyReluLayer(0.1,'Name','relu1')
    
    fullyConnectedLayer(64,'Name','fc2')
    batchNormalizationLayer('Name','bnorm2')
    leakyReluLayer(0.1,'Name','relu2')

    fullyConnectedLayer(128,'Name','fc3')
    batchNormalizationLayer('Name','bnorm3')
    leakyReluLayer(0.1,'Name','relu3')
    
    fullyConnectedLayer(256,'Name','fc4')
    batchNormalizationLayer('Name','bnorm4')
    leakyReluLayer(0.1,'Name','relu4')
    
    fullyConnectedLayer(512,'Name','fc5')
    batchNormalizationLayer('Name','bnorm5')
    leakyReluLayer(0.1,'Name','relu5')
    
    fullyConnectedLayer(256,'Name','fc6')
    batchNormalizationLayer('Name','bnorm6')
    leakyReluLayer(0.1,'Name','relu6')
    
    fullyConnectedLayer(128,'Name','fc7')
    batchNormalizationLayer('Name','bnorm7')
    leakyReluLayer(0.1,'Name','relu7')
    
    fullyConnectedLayer(64,'Name','fc8')
    batchNormalizationLayer('Name','bnorm8')
    leakyReluLayer(0.1,'Name','relu8')
    
    fullyConnectedLayer(32,'Name','fc9')
    batchNormalizationLayer('Name','bnorm9')
    leakyReluLayer(0.1,'Name','relu9')

    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(16,'Name','fc10')
    regressionLayer('Name','reg')
];

lgraph = layerGraph(layers); 


end % if strcmp(netType,...)

end % function chooseNetwork
