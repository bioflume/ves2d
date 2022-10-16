clear; clc;
addpath ./necessaryMatFiles/

for idt = 1 : 13
mul = 2.^(0:12)';
dts = 1E-6*mul;
fileName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/workspace/completeData/n256Dt' num2str(dts(idt)) 'Kb1nModes64PCAData.mat'];
load(fileName)

load ./necessaryMatFiles/pcaBasisNewest

% PREDICT NEW SHAPE

nComp = 16; % use 16 principal components
compStart = 1;
compEnd = compStart-1+16;

netType = 'fcLarge2'; 


fileName = ['fcPCArelaxN256Dt' num2str(dts(idt)) 'Kb1nModes1to16_fcMedium'];

nves = nInstances;

nTrain = floor(0.95*nves);
nTest = nves-nTrain;
Xtrain = zeros(nComp,1,1,nTrain);
Ytrain = zeros(nTrain,nComp);
Xtest = zeros(nComp,1,1,nTest);
Ytest = zeros(nTest,nComp);

if nInstances > 0
if 1 % randomly choose idcs or load already randomly chosen one
  testIdcs = randperm(nves,nTest)';
  trainIdcs = zeros(nTrain,1);

  idx = 1;
  for k = 1 : nves
    if ~any(k==testIdcs)
      trainIdcs(idx) = k;
      idx = idx+1;
    end
  end
  save setOfTrainTestIdcsDt1E3 testIdcs trainIdcs
else
  load setOfTrainTestIdcsStand1step
  %load setOfTrainTestIdcsDt5E4
  %load setOfTrainTestIdcsDt1E3
end

% PREPARE INPUT AND OUTPUT

for i = 1 : nTrain
  Xtrain(:,:,:,i) = reshape(XoldCoeffs(trainIdcs(i),compStart:compEnd),...
      nComp,1,1);
  Ytrain(i,:) = XnewCoeffs(trainIdcs(i),compStart:compEnd);    
end
for i = 1:nTest
  Xtest(:,:,:,i) = reshape(XoldCoeffs(testIdcs(i),compStart:compEnd),...
      nComp,1,1);
  Ytest(i,:) = XnewCoeffs(testIdcs(i),compStart:compEnd);
end


% Normalize or scale inputs to a certain range
% Compute mean and stdev of channels
chan1vals = Xtrain(:,1,1,:); 
muChan1 = mean(chan1vals(:)); sdevChan1 = std(chan1vals(:));   
% Normalize
scale = 1; offset = 0;

Xtrain(:,1,1,:) = scale*(Xtrain(:,1,1,:)-muChan1)/sdevChan1+offset;
Ytrain = scale*(Ytrain-muChan1)/sdevChan1+offset;

Xtest(:,1,1,:) = scale*(Xtest(:,1,1,:)-muChan1)/sdevChan1+offset;
Ytest = scale*(Ytest-muChan1)/sdevChan1+offset;

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
    'InitialLearnRate',5e-5, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.2, ...
    'LearnRateDropPeriod',25, ...
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

%YOPred = Ypred;
%YOTrainPred = YtrainPred;

YOPred = (Ypred-offset)*sdevChan1/scale+muChan1;
YOTrainPred = (YtrainPred-offset)*sdevChan1/scale+muChan1;
Ytest = (Ytest-offset)*sdevChan1/scale+muChan1;
Ytrain = (Ytrain-offset)*sdevChan1/scale+muChan1;

err = sqrt(1/nComp*sum((YOPred-Ytest).^2,2))./...
    sqrt(1/nComp*sum(Ytest.^2,2));
disp(['Max error is ' num2str(max(err))])
disp(['Min error is ' num2str(min(err))])
disp(['Mean error is ' num2str(mean(err))])

netFileName = [fileName '_tstepFCNet_flow.mat'];
errMax = max(err);
errMin = min(err);
errMean = mean(err);
save(netFileName,'net','muChan1','sdevChan1','scale','offset','nComp','errMax','errMin','errMean')
end
end
%save(['/workspace/gokberk/' fileName '_w1step.mat'],'YOPred','Ytest',...
%    'Xtrain','Xtest','YOTrainPred','Ytrain')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lgraph,numParams] = chooseNetwork(netType,M,mu,nComp)



if strcmp(netType,'fcMedium')
numParams = (16*32 + 32*64 + 64*128 + 128*64 + 64*32 + 32* 16 + 32 + 64 + 128+ 64+ 32+ 16);

layers = [imageInputLayer([16 1 1],'Name','input','Normalization','none')

    fullyConnectedLayer(32,'Name','fc1')
    batchNormalizationLayer('Name','bnorm1')
    
    fullyConnectedLayer(64,'Name','fc2')
    batchNormalizationLayer('Name','bnorm2')
    
    fullyConnectedLayer(128,'Name','fc3')
    batchNormalizationLayer('Name','bnorm3')
    
    fullyConnectedLayer(64,'Name','fc4')
    batchNormalizationLayer('Name','bnorm4')
    
    fullyConnectedLayer(32,'Name','fc5')
    batchNormalizationLayer('Name','bnorm5')
    
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(16,'Name','fc6')
    regressionLayer('Name','reg')
    
];

lgraph = layerGraph(layers); 

elseif strcmp(netType,'fcLarge2') % 69,472
    
numParams = (16*48+48*96+96*128+128*192+192*128+...
  128*96+96*16+48+96+128+192+128+96+16);

layers = [imageInputLayer([16 1 1],'Name','input','Normalization','none')

    fullyConnectedLayer(48,'Name','fc1')
    batchNormalizationLayer('Name','bnorm1')

    additionLayer(2,'Name','add1')
    
    fullyConnectedLayer(96,'Name','fc2')
    batchNormalizationLayer('Name','bnorm2')

    additionLayer(2,'Name','add2')
    
    fullyConnectedLayer(128,'Name','fc3')
    batchNormalizationLayer('Name','bnorm3')

    additionLayer(2,'Name','add3')
    
    fullyConnectedLayer(192,'Name','fc4')
    batchNormalizationLayer('Name','bnorm4')

    additionLayer(2,'Name','add4')
    
    fullyConnectedLayer(128,'Name','fc5')
    batchNormalizationLayer('Name','bnorm5')

    additionLayer(2,'Name','add5')
    
    fullyConnectedLayer(96,'Name','fc6')
    batchNormalizationLayer('Name','bnorm6')

    additionLayer(2,'Name','add6')

    fullyConnectedLayer(16,'Name','fc7')
    regressionLayer('Name','reg')
];

lgraph = layerGraph(layers); 
lgraph = connectLayers(lgraph,'bnorm1','add1/in2');
lgraph = connectLayers(lgraph,'bnorm2','add2/in2');
lgraph = connectLayers(lgraph,'bnorm3','add3/in2');
lgraph = connectLayers(lgraph,'bnorm4','add4/in2');
lgraph = connectLayers(lgraph,'bnorm5','add5/in2');
lgraph = connectLayers(lgraph,'bnorm6','add6/in2');



end % if strcmp(netType,...)

end % function chooseNetwork
