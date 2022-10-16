clear; clc;
addpath ./necessaryMatFiles/

%load /workspace/gokberk/only1timeStepPCAData.mat
%load /workspace/gokberk/only1timeStep_inputPCACoeffs.mat
%load ./necessaryMatFiles/pcaCoeffsBasis1step.mat
load ./necessaryMatFiles/pcaBasisNewest.mat
load /workspace/gokberk/n256Dt1E4Kb1_inputPCAs.mat
load /workspace/gokberk/NEW_tenMatOnBendN256AllData.mat
%load /workspace/gokberk/tenMatOnBendN256AllData.mat

% predict inv(Div*G*Ten)*Div*G*(-Ben)*X, 8 modes (real and imag. comps)
nComp = 16; % use 16 principal components
nmodes = 32;
predModes = [(1:nmodes/2)';(N-nmodes/2+1:N)'];

netType = 'fcXlargeTen' 

fileName = 'NEW_fcPCAtenMaTonBendN256_xlargeTen32modes';

nves = nInstances;

nTrain = floor(0.95*nves);
nTest = nves-nTrain;
Xtrain = zeros(nComp,1,1,nTrain);
Ytrain = zeros(nTrain,nmodes*2);
Xtest = zeros(nComp,1,1,nTest);
Ytest = zeros(nTest,nmodes*2);


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
  save setOfTrainTestIdcsN96Clean testIdcs trainIdcs
else
  %load setOfTrainTestIdcsN96
  load setOfTrainTestIdcsStand1step
end

% PREPARE INPUT AND OUTPUT

for i = 1 : nTrain
  Xtrain(:,:,:,i) = reshape(XoldCoeffs(trainIdcs(i),1:nComp),...
      nComp,1,1);
  Ytrain(i,:) = [benRealStore(predModes,trainIdcs(i));...
      benImagStore(predModes,trainIdcs(i))];    
end
for i = 1:nTest
  Xtest(:,:,:,i) = reshape(XoldCoeffs(testIdcs(i),1:nComp),...
      nComp,1,1);
  Ytest(i,:) = [benRealStore(predModes,testIdcs(i));...
      benImagStore(predModes,testIdcs(i))];    
end


% Normalize or scale inputs to a certain range
% Compute mean and stdev of channels
chan1vals = Xtrain(:,1,1,:); 
muChan1 = mean(chan1vals(:)); sdevChan1 = std(chan1vals(:));   
% Normalize inputs only
scale = 1; offset = 0;
Xtrain(:,1,1,:) = scale*(Xtrain(:,1,1,:)-muChan1)/sdevChan1+offset;
Xtest(:,1,1,:) = scale*(Xtest(:,1,1,:)-muChan1)/sdevChan1+offset;

% LOAD NETWORK
[lgraph,numParams] = chooseNetwork(netType,evects,colMeans',nComp);

disp(['There are ' num2str(numParams) ' # of parameters.'])
disp(['There are ' num2str(nTrain) ' # of training instances.'])

miniBatchSize  = 256;
validationFrequency = floor(numel(Ytrain)/miniBatchSize);
options = trainingOptions('sgdm', ...
    'Momentum',0.9,...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',20, ...
    'InitialLearnRate',1e-3, ...
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
YOPred = predict(net,Xtest); % compare with Ytest
YOTrainPred = predict(net,Xtrain);

err = sqrt(1/nComp*sum((YOPred-Ytest).^2,2))./...
    sqrt(1/nComp*sum(Ytest.^2,2));
disp(['Max error is ' num2str(max(err))])
disp(['Min error is ' num2str(min(err))])
disp(['Mean error is ' num2str(mean(err))])

netFileName = [fileName '_FCNet_w1step.mat'];
save(netFileName,'net','muChan1','sdevChan1','scale','offset','nComp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lgraph,numParams] = chooseNetwork(netType,M,mu,nComp)

if strcmp(netType,'fcXlarge') % 85,920 params
    
numParams = (16*48+48*96+96*128+128*256+256*128+128*64+48+96+128+256+128+64);

layers = [imageInputLayer([16 1 1],'Name','input','Normalization','none')

    fullyConnectedLayer(48,'Name','fc1')
    batchNormalizationLayer('Name','bnorm1')
    leakyReluLayer(0.1,'Name','relu1')

    additionLayer(2,'Name','add1')
    
    fullyConnectedLayer(96,'Name','fc2')
    batchNormalizationLayer('Name','bnorm2')
    leakyReluLayer(0.1,'Name','relu2')

    additionLayer(2,'Name','add2')
    
    fullyConnectedLayer(128,'Name','fc3')
    batchNormalizationLayer('Name','bnorm3')
    leakyReluLayer(0.1,'Name','relu3')

    additionLayer(2,'Name','add3')
    
    fullyConnectedLayer(256,'Name','fc4')
    batchNormalizationLayer('Name','bnorm4')
    leakyReluLayer(0.1,'Name','relu4')

    additionLayer(2,'Name','add4')
    
    fullyConnectedLayer(128,'Name','fc5')
    batchNormalizationLayer('Name','bnorm5')
    leakyReluLayer(0.1,'Name','relu5')

    additionLayer(2,'Name','add5')
    
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(64,'Name','fc6')
    regressionLayer('Name','reg')
];

lgraph = layerGraph(layers); 
lgraph = connectLayers(lgraph,'bnorm1','add1/in2');
lgraph = connectLayers(lgraph,'bnorm2','add2/in2');
lgraph = connectLayers(lgraph,'bnorm3','add3/in2');
lgraph = connectLayers(lgraph,'bnorm4','add4/in2');
lgraph = connectLayers(lgraph,'bnorm5','add5/in2');

elseif strcmp(netType,'fcXlargeTen') % 85,920 params
    
numParams = (16*64+64*128+128*256+256*128+128*64+64+128+256+128+64);

layers = [imageInputLayer([16 1 1],'Name','input','Normalization','none')

    fullyConnectedLayer(64,'Name','fc1')
    batchNormalizationLayer('Name','bnorm1')
    leakyReluLayer(0.1,'Name','relu1')

    additionLayer(2,'Name','add1')
    
    fullyConnectedLayer(128,'Name','fc2')
    batchNormalizationLayer('Name','bnorm2')
    leakyReluLayer(0.1,'Name','relu2')

    additionLayer(2,'Name','add2')
    
    fullyConnectedLayer(256,'Name','fc3')
    batchNormalizationLayer('Name','bnorm3')
    leakyReluLayer(0.1,'Name','relu3')

    additionLayer(2,'Name','add3')
    
    fullyConnectedLayer(128,'Name','fc4')
    batchNormalizationLayer('Name','bnorm4')
    leakyReluLayer(0.1,'Name','relu4')

    additionLayer(2,'Name','add4')
    
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(64,'Name','fc5')
    regressionLayer('Name','reg')
];

lgraph = layerGraph(layers); 
lgraph = connectLayers(lgraph,'bnorm1','add1/in2');
lgraph = connectLayers(lgraph,'bnorm2','add2/in2');
lgraph = connectLayers(lgraph,'bnorm3','add3/in2');
lgraph = connectLayers(lgraph,'bnorm4','add4/in2');


elseif strcmp(netType,'fcLarge') % 72,624 parameters
    
numParams = (16*48+48*96+96*128+128*256+256*64+64*48+48*32+32*16 ...
  +48+96+128+256+64+48+32+16);

layers = [imageInputLayer([16 1 1],'Name','input','Normalization','none')

    fullyConnectedLayer(48,'Name','fc1')
    batchNormalizationLayer('Name','bnorm1')
    leakyReluLayer(0.1,'Name','relu1')

    additionLayer(2,'Name','add1')
    
    fullyConnectedLayer(96,'Name','fc2')
    batchNormalizationLayer('Name','bnorm2')
    leakyReluLayer(0.1,'Name','relu2')

    additionLayer(2,'Name','add2')
    
    fullyConnectedLayer(128,'Name','fc3')
    batchNormalizationLayer('Name','bnorm3')
    leakyReluLayer(0.1,'Name','relu3')

    additionLayer(2,'Name','add3')
    
    fullyConnectedLayer(256,'Name','fc4')
    batchNormalizationLayer('Name','bnorm4')
    leakyReluLayer(0.1,'Name','relu4')

    additionLayer(2,'Name','add4')
    
    fullyConnectedLayer(64,'Name','fc5')
    batchNormalizationLayer('Name','bnorm5')
    leakyReluLayer(0.1,'Name','relu5')

    additionLayer(2,'Name','add5')
    
    fullyConnectedLayer(48,'Name','fc6')
    batchNormalizationLayer('Name','bnorm6')
    leakyReluLayer(0.1,'Name','relu6')

    additionLayer(2,'Name','add6')
    
    fullyConnectedLayer(32,'Name','fc7')
    batchNormalizationLayer('Name','bnorm7')
    leakyReluLayer(0.1,'Name','relu7')

    additionLayer(2,'Name','add7')
    
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(16,'Name','fc8')
    regressionLayer('Name','reg')
];

lgraph = layerGraph(layers); 
lgraph = connectLayers(lgraph,'bnorm1','add1/in2');
lgraph = connectLayers(lgraph,'bnorm2','add2/in2');
lgraph = connectLayers(lgraph,'bnorm3','add3/in2');
lgraph = connectLayers(lgraph,'bnorm4','add4/in2');
lgraph = connectLayers(lgraph,'bnorm5','add5/in2');
lgraph = connectLayers(lgraph,'bnorm6','add6/in2');
lgraph = connectLayers(lgraph,'bnorm7','add7/in2');

elseif strcmp(netType,'fcLargeNoAdd')
numParams = (16*48+48*96+96*128+128*256+256*64+64*48+48*32+32*16 ...
  +48+96+128+256+64+48+32+16);

layers = [imageInputLayer([16 1 1],'Name','input','Normalization','none')

    fullyConnectedLayer(48,'Name','fc1')
    batchNormalizationLayer('Name','bnorm1')
    leakyReluLayer(0.1,'Name','relu1')

    fullyConnectedLayer(96,'Name','fc2')
    batchNormalizationLayer('Name','bnorm2')
    leakyReluLayer(0.1,'Name','relu2')

    fullyConnectedLayer(128,'Name','fc3')
    batchNormalizationLayer('Name','bnorm3')
    leakyReluLayer(0.1,'Name','relu3')
    
    fullyConnectedLayer(256,'Name','fc4')
    batchNormalizationLayer('Name','bnorm4')
    leakyReluLayer(0.1,'Name','relu4')

    fullyConnectedLayer(64,'Name','fc5')
    batchNormalizationLayer('Name','bnorm5')
    leakyReluLayer(0.1,'Name','relu5')

    fullyConnectedLayer(48,'Name','fc6')
    batchNormalizationLayer('Name','bnorm6')
    leakyReluLayer(0.1,'Name','relu6')

    fullyConnectedLayer(32,'Name','fc7')
    batchNormalizationLayer('Name','bnorm7')
    leakyReluLayer(0.1,'Name','relu7')

    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(16,'Name','fc8')
    regressionLayer('Name','reg')
];

lgraph = layerGraph(layers); 

elseif strcmp(netType,'fcLarge2') % 69,472
    
numParams = (16*48+48*96+96*128+128*192+192*128+...
  128*96+96*16+48+96+128+192+128+96+16);

layers = [imageInputLayer([16 1 1],'Name','input','Normalization','none')

    fullyConnectedLayer(48,'Name','fc1')
    batchNormalizationLayer('Name','bnorm1')
    leakyReluLayer(0.1,'Name','relu1')

    additionLayer(2,'Name','add1')
    
    fullyConnectedLayer(96,'Name','fc2')
    batchNormalizationLayer('Name','bnorm2')
    leakyReluLayer(0.1,'Name','relu2')

    additionLayer(2,'Name','add2')
    
    fullyConnectedLayer(128,'Name','fc3')
    batchNormalizationLayer('Name','bnorm3')
    leakyReluLayer(0.1,'Name','relu3')

    additionLayer(2,'Name','add3')
    
    fullyConnectedLayer(192,'Name','fc4')
    batchNormalizationLayer('Name','bnorm4')
    leakyReluLayer(0.1,'Name','relu4')

    additionLayer(2,'Name','add4')
    
    fullyConnectedLayer(128,'Name','fc5')
    batchNormalizationLayer('Name','bnorm5')
    leakyReluLayer(0.1,'Name','relu5')

    additionLayer(2,'Name','add5')
    
    fullyConnectedLayer(96,'Name','fc6')
    batchNormalizationLayer('Name','bnorm6')
    leakyReluLayer(0.1,'Name','relu6')

    additionLayer(2,'Name','add6')

    dropoutLayer(0.5,'Name','dropout')
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

elseif strcmp(netType,'fcMedium') % 60,080 parameters
    
numParams = (16*32*2+32*48*2+48*64*2+64*96*2+96*192*2)+...
    (2*48+64*2+96*2+2*32+192+16);

layers = [imageInputLayer([16 1 1],'Name','input','Normalization','none')

    fullyConnectedLayer(32,'Name','fc1')
    batchNormalizationLayer('Name','bnorm1')
    leakyReluLayer(0.1,'Name','relu1')

    additionLayer(2,'Name','add1')
    
    fullyConnectedLayer(48,'Name','fc2')
    batchNormalizationLayer('Name','bnorm2')
    leakyReluLayer(0.1,'Name','relu2')

    additionLayer(2,'Name','add2')
    
    fullyConnectedLayer(64,'Name','fc3')
    batchNormalizationLayer('Name','bnorm3')
    leakyReluLayer(0.1,'Name','relu3')

    additionLayer(2,'Name','add3')
    
    fullyConnectedLayer(96,'Name','fc4')
    batchNormalizationLayer('Name','bnorm4')
    leakyReluLayer(0.1,'Name','relu4')

    additionLayer(2,'Name','add4')
    
    fullyConnectedLayer(192,'Name','fc5')
    batchNormalizationLayer('Name','bnorm5')
    leakyReluLayer(0.1,'Name','relu5')

    additionLayer(2,'Name','add5')
    
    fullyConnectedLayer(96,'Name','fc6')
    batchNormalizationLayer('Name','bnorm6')
    leakyReluLayer(0.1,'Name','relu6')

    additionLayer(2,'Name','add6')
    
    fullyConnectedLayer(64,'Name','fc7')
    batchNormalizationLayer('Name','bnorm7')
    leakyReluLayer(0.1,'Name','relu7')

    additionLayer(2,'Name','add7')
    
    fullyConnectedLayer(48,'Name','fc8')
    batchNormalizationLayer('Name','bnorm8')
    leakyReluLayer(0.1,'Name','relu8')

    additionLayer(2,'Name','add8')
    
    fullyConnectedLayer(32,'Name','fc9')
    batchNormalizationLayer('Name','bnorm9')
    leakyReluLayer(0.1,'Name','relu9')

    additionLayer(2,'Name','add9')
    
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(16,'Name','fc10')
    regressionLayer('Name','reg')
];

lgraph = layerGraph(layers); 
lgraph = connectLayers(lgraph,'bnorm1','add1/in2');
lgraph = connectLayers(lgraph,'bnorm2','add2/in2');
lgraph = connectLayers(lgraph,'bnorm3','add3/in2');
lgraph = connectLayers(lgraph,'bnorm4','add4/in2');
lgraph = connectLayers(lgraph,'bnorm5','add5/in2');
lgraph = connectLayers(lgraph,'bnorm6','add6/in2');
lgraph = connectLayers(lgraph,'bnorm7','add7/in2');  
lgraph = connectLayers(lgraph,'bnorm8','add8/in2');  
lgraph = connectLayers(lgraph,'bnorm9','add9/in2');  

end % if strcmp(netType,...)

end % function chooseNetwork
