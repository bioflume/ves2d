clear; clc;
addpath ./necessaryMatFiles/

% which modes to be trained (range: modeStart:modeEnd)
modeStart = 19;
modeEnd = 21;
disp(['Training for the modes in [' num2str(modeStart) ...
 ', ' num2str(modeEnd) ']' ])

% input PCA or FFT coeffs
iFFT = ~true;
% important modes for the input shape (FFT modes)
%shapeModes = [2;256;4;255;6;3;5;254]; 
shapeModes = [2;96;4;95;6;3;5;94];

% number of PCA components
nComp = 16;

% NETWORK ARCHITECTURE
netType = 'fcLarge' 
% or 'xlarge','large' ---> ResNet
% or 'fcXlarge','fcLarge','fcMedium' ---> Fully-Connected Layers

% OUTPUT INFORMATION
N = 256;
%N = 96;
nmodes = 24;
activeModes = [(1:nmodes/2)';(N-nmodes/2+1:N)'];
outputSize = 96; % 1-48: real, 49-96: imag. components 
                 % 1-24: realx, 25-48: realy, 49-72: imagx, 72-96: imagy

% load MVinf data and FFT coeffs
%load /workspace/gokberk/velocityTrainFFTAllData.mat
load /workspace/gokberk/NEW_MtimesFFTBasisN256All24modesData.mat
if ~iFFT
  %clear shapeModesStore;
  %load /workspace/gokberk/only1timeStep_inputPCACoeffs.mat
  load /workspace/gokberk/n256Dt1E4Kb1_inputPCAs
  %load /workspace/gokberk/only1timeStepN96PCAData.mat
  %clear XnewCoeffs;
end

% number of instances
nves = nInstances;
nTrain = floor(0.95*nves);
nTest = nves-nTrain;
% randomly choose idcs or load already randomly chosen one
if 0
  testIdcs = randperm(nves,nTest)';
  trainIdcs = zeros(nTrain,1);
  idx = 1;
  for k = 1 : nves
    if ~any(k==testIdcs)
      trainIdcs(idx) = k;
      idx = idx+1;
    end
  end
  save setOfTrainTestIdcsMtimesFFTN96 testIdcs trainIdcs
else
  %load setOfTrainTestIdcsN96
  load setOfTrainTestIdcsStand1step
%   load setOfTrainTestIdcsRelaxEvolveStand
end

% INITIALIZE INPUTS
if iFFT
  % if FFT coeffs are inputs, than we may have two channels
  if strcmp(netType,'xlarge') || strcmp(netType,'large')
    Xtrain = zeros(numel(shapeModes),1,2,nTrain);
    Xtest = zeros(numel(shapeModes),1,2,nTest);
  else
    Xtrain = zeros(2*numel(shapeModes),1,1,nTrain);
    Xtest = zeros(2*numel(shapeModes),1,1,nTest);
  end
else
  % if PCA coeffs are inputs, then we have one channel for 16 PCA coeffs
  Xtrain = zeros(nComp,1,1,nTrain);
  Xtest = zeros(nComp,1,1,nTest);        
end % if iFFT

% PREPARE INPUT (SAME FOR ALL MODES pmode)
% Training set:
for i = 1 : nTrain
  % INPUT:  
  if iFFT  
    if strcmp(netType,'xlarge') || strcmp(netType,'large')
    Xtrain(:,:,:,i) = reshape(shapeModesStore(shapeModes,:,trainIdcs(i)),...
        numel(shapeModes),1,2);
    else
    Xtrain(:,:,:,i) = reshape(shapeModesStore(shapeModes,:,trainIdcs(i)),...
      2*numel(shapeModes),1,1);
    end % if strcmp(netType,'xlarge') || strcmp(netType,'large')
  else
    Xtrain(:,:,:,i) = reshape(XoldCoeffs(trainIdcs(i),1:nComp),nComp,1,1);     
  end % if iFFT
end

% Testing set:
for i = 1:nTest
  % INPUT:
  if iFFT
    if strcmp(netType,'xlarge') || strcmp(netType,'large')
    Xtest(:,:,:,i) = reshape(shapeModesStore(shapeModes,:,testIdcs(i)),...
        numel(shapeModes),1,2);
    else
    Xtest(:,:,:,i) = reshape(shapeModesStore(shapeModes,:,testIdcs(i)),...
      2*numel(shapeModes),1,1);
    end % if strcmp(netType,'xlarge') || strcmp(netType,'large')
  else
    Xtest(:,:,:,i) = reshape(XoldCoeffs(testIdcs(i),1:nComp),nComp,1,1);         
  end % if iFFT
end
% Clear these matrices:
if iFFT; clear shapeModesStore; end
if ~iFFT; clear XoldCoeffs; end

% Normalize inputs to a certain range
% Compute mean and stdev of channels
if numel(Xtrain(1,1,:,1))==2 % if there are two channels
  chan1vals = Xtrain(:,1,1,:); 
  muChan1 = mean(chan1vals(:)); sdevChan1 = std(chan1vals(:));
  chan2vals = Xtrain(:,1,2,:); 
  muChan2 = mean(chan2vals(:)); sdevChan2 = std(chan2vals(:));

  % Normalize
  scale = 1; offset = 0;

  Xtrain(:,1,1,:) = scale*(Xtrain(:,1,1,:)-muChan1)/sdevChan1+offset;
  Xtrain(:,1,2,:) = scale*(Xtrain(:,1,2,:)-muChan2)/sdevChan2+offset;

  Xtest(:,1,1,:) = scale*(Xtest(:,1,1,:)-muChan1)/sdevChan1+offset;
  Xtest(:,1,2,:) = scale*(Xtest(:,1,2,:)-muChan2)/sdevChan2+offset;
else
  chan1vals = Xtrain(:,1,1,:);
  muChan1 = mean(chan1vals(:)); sdevChan1 = std(chan1vals(:));
  muChan2 = []; sdevChan2 = [];
  
  scale = 1; offset = 0;

  Xtrain(:,1,1,:) = scale*(Xtrain(:,1,1,:)-muChan1)/sdevChan1+offset;
  Xtest(:,1,1,:) = scale*(Xtest(:,1,1,:)-muChan1)/sdevChan1+offset;
end

% LOAD NETWORK
[lgraph,numParams] = chooseNetwork(netType);

disp(['There are ' num2str(numParams) ' # of parameters.'])
disp(['There are ' num2str(nTrain) ' # of training instances.'])

for imode = modeStart : modeEnd

pmode = activeModes(imode) % predict pmode of the output

if iFFT
fileName = ['./NEWn256Mtimes24modesFFTNets/velPredFFTin_mode' num2str(pmode)];
else
fileName = ['./NEWn256Mtimes24modesFFTNets/velPredPCAin_mode' num2str(pmode)];    
end

% INITIALIZE OUTPUTS
Ytrain = zeros(nTrain,outputSize);
Ytest = zeros(nTest,outputSize);

% PREPARE OUTPUTS
% Training set:
% DOWNSAMPLE 512 to 24 points,stack x-y of real into multReal, 
% x-y of imag into multImag
for i = 1 : nTrain    
  multReal = [interpft(zRealStore(1:end/2,imode,trainIdcs(i)),outputSize/4);...
      interpft(zRealStore(end/2+1:end,imode,trainIdcs(i)),outputSize/4)];
  multImag = [interpft(zImagStore(1:end/2,imode,trainIdcs(i)),outputSize/4);...
      interpft(zImagStore(end/2+1:end,imode,trainIdcs(i)),outputSize/4)];
  Ytrain(i,:) = [multReal;multImag];  
end

% Testing set:
for i = 1:nTest
  multReal = [interpft(zRealStore(1:end/2,imode,testIdcs(i)),outputSize/4);...
      interpft(zRealStore(end/2+1:end,imode,testIdcs(i)),outputSize/4)];
  multImag = [interpft(zImagStore(1:end/2,imode,testIdcs(i)),outputSize/4);...
      interpft(zImagStore(end/2+1:end,imode,testIdcs(i)),outputSize/4)];
  Ytest(i,:) = [multReal;multImag];  
end

miniBatchSize  = 256;
validationFrequency = floor(numel(Ytrain)/miniBatchSize);
options = trainingOptions('sgdm', ...
    'Momentum',0.9,...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',20, ...
    'InitialLearnRate',5e-4, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',10, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{Xtest,Ytest},...
    'ValidationFrequency',validationFrequency,...
    'Verbose',true,...
    'VerboseFrequency',5);

% TRAIN NETWORK
net = trainNetwork(Xtrain,Ytrain,lgraph,options);

% PREDICT TEST DATA
YOPred = predict(net,Xtest);
YTrainPred = predict(net,Xtrain);

% COMPUTE ERROR IN PREDICTED RESULTS
err = sqrt(2/outputSize*sum((YOPred(:,1:end/2)-Ytest(:,1:end/2)).^2+...
    (YOPred(:,end/2+1:end)-Ytest(:,end/2+1:end)).^2,2))./...
    sqrt(2/outputSize*sum(Ytest.^2,2));
disp(['Max error is ' num2str(max(err))])
disp(['Min error is ' num2str(min(err))])
disp(['Mean error is ' num2str(mean(err))])

netFileName = [fileName '_net_w1step.mat'];
save(netFileName,'net','muChan1','sdevChan1','muChan2',...
'outputSize','sdevChan2','scale','offset')

% save([fileName '_w1step.mat'],'YOPred','YOTest','Ytest',...
%   'YPredicted','Xtrain','Xtest','pmode','channel','YOTrainPred','Ytrain')
end % for imode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lgraph,numParams] = chooseNetwork(netType)

filterSizes = [];

if strcmp(netType,'xlarge') % 76,368 params
    
filterSizes(1,:) = [1 1 8];
filterSizes(2,:) = [3 1 8];
filterSizes(3,:) = [1 1 16];
filterSizes(4,:) = [3 1 16];
filterSizes(5,:) = [1 1 32];
filterSizes(6,:) = [3 1 32];
filterSizes(7,:) = [1 1 64];
filterSizes(8,:) = [3 1 64];
filterSizes(9,:) = [1 1 128];

layers = [imageInputLayer([8 1 2],'Name','input','Normalization','none')

    convolution2dLayer(filterSizes(1,1:2),filterSizes(1,3),'Name','conv1')
    batchNormalizationLayer('Name','bnorm1')
    leakyReluLayer(0.1,'Name','relu1')
     
    additionLayer(2,'Name','add1')

    convolution2dLayer(filterSizes(2,1:2),filterSizes(2,3),'Name','conv2','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm2')
    leakyReluLayer(0.1,'Name','relu2')
   
    additionLayer(2,'Name','add2')

    convolution2dLayer(filterSizes(3,1:2),filterSizes(3,3),'Name','conv3')
    batchNormalizationLayer('Name','bnorm3')
    leakyReluLayer(0.1,'Name','relu3')
 
    additionLayer(2,'Name','add3')
   
    convolution2dLayer(filterSizes(4,1:2),filterSizes(4,3),'Name','conv4','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm4')
    leakyReluLayer(0.1,'Name','relu4')
   
    additionLayer(2,'Name','add4')

    convolution2dLayer(filterSizes(5,1:2),filterSizes(5,3),'Name','conv5')
    batchNormalizationLayer('Name','bnorm5')
    leakyReluLayer(0.1,'Name','relu5')
    
    additionLayer(2,'Name','add5')

    convolution2dLayer(filterSizes(6,1:2),filterSizes(6,3),'Name','conv6')
    batchNormalizationLayer('Name','bnorm6')
    leakyReluLayer(0.1,'Name','relu6')
    
    additionLayer(2,'Name','add6')

    convolution2dLayer(filterSizes(7,1:2),filterSizes(7,3),'Name','conv7')
    batchNormalizationLayer('Name','bnorm7')
    leakyReluLayer(0.1,'Name','relu7')
    
    additionLayer(2,'Name','add7')
    
    convolution2dLayer(filterSizes(8,1:2),filterSizes(8,3),'Name','conv8')
    batchNormalizationLayer('Name','bnorm8')
    leakyReluLayer(0.1,'Name','relu8')
    
    additionLayer(2,'Name','add8')

    convolution2dLayer(filterSizes(9,1:2),filterSizes(9,3),'Name','conv9')
    batchNormalizationLayer('Name','bnorm9')
    leakyReluLayer(0.1,'Name','relu9')
    
    additionLayer(2,'Name','add9')

    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(96,'Name','fc')
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

elseif strcmp(netType,'large') % 62,032 parameters
filterSizes = [];
filterSizes(1,:) = [1 1 8];
filterSizes(2,:) = [3 1 8];
filterSizes(3,:) = [1 1 16];
filterSizes(4,:) = [3 1 16];
filterSizes(5,:) = [1 1 32];
filterSizes(6,:) = [3 1 32];
filterSizes(7,:) = [1 1 64];
filterSizes(8,:) = [3 1 64];
filterSizes(9,:) = [1 1 96];

layers = [imageInputLayer([8 1 2],'Name','input','Normalization','none')

    convolution2dLayer(filterSizes(1,1:2),filterSizes(1,3),'Name','conv1')
    batchNormalizationLayer('Name','bnorm1')
    leakyReluLayer(0.1,'Name','relu1')
     
    additionLayer(2,'Name','add1')

    convolution2dLayer(filterSizes(2,1:2),filterSizes(2,3),'Name','conv2','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm2')
    leakyReluLayer(0.1,'Name','relu2')
   
    additionLayer(2,'Name','add2')

    convolution2dLayer(filterSizes(3,1:2),filterSizes(3,3),'Name','conv3')
    batchNormalizationLayer('Name','bnorm3')
    leakyReluLayer(0.1,'Name','relu3')
 
    additionLayer(2,'Name','add3')
   
    convolution2dLayer(filterSizes(4,1:2),filterSizes(4,3),'Name','conv4','Padding',[1 0])
    batchNormalizationLayer('Name','bnorm4')
    leakyReluLayer(0.1,'Name','relu4')
   
    additionLayer(2,'Name','add4')

    convolution2dLayer(filterSizes(5,1:2),filterSizes(5,3),'Name','conv5')
    batchNormalizationLayer('Name','bnorm5')
    leakyReluLayer(0.1,'Name','relu5')
    
    additionLayer(2,'Name','add5')

    convolution2dLayer(filterSizes(6,1:2),filterSizes(6,3),'Name','conv6')
    batchNormalizationLayer('Name','bnorm6')
    leakyReluLayer(0.1,'Name','relu6')
    
    additionLayer(2,'Name','add6')

    convolution2dLayer(filterSizes(7,1:2),filterSizes(7,3),'Name','conv7')
    batchNormalizationLayer('Name','bnorm7')
    leakyReluLayer(0.1,'Name','relu7')
    
    additionLayer(2,'Name','add7')
    
    convolution2dLayer(filterSizes(8,1:2),filterSizes(8,3),'Name','conv8')
    batchNormalizationLayer('Name','bnorm8')
    leakyReluLayer(0.1,'Name','relu8')
    
    additionLayer(2,'Name','add8')

    convolution2dLayer(filterSizes(9,1:2),filterSizes(9,3),'Name','conv9')
    batchNormalizationLayer('Name','bnorm9')
    leakyReluLayer(0.1,'Name','relu9')
    
    additionLayer(2,'Name','add9')
 
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(96,'Name','fc')
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

elseif strcmp(netType,'fcXlarge') % 80,608 parameters

numParams = (16*128+4*128*128+128*96)+5*128+96;

layers = [imageInputLayer([16 1 1],'Name','input','Normalization','none')

    fullyConnectedLayer(128,'Name','fc1')
    batchNormalizationLayer('Name','bnorm1')
    leakyReluLayer(0.1,'Name','relu1')

    additionLayer(2,'Name','add1')
    
    fullyConnectedLayer(128,'Name','fc2')
    batchNormalizationLayer('Name','bnorm2')
    leakyReluLayer(0.1,'Name','relu2')

    additionLayer(2,'Name','add2')
    
    fullyConnectedLayer(128,'Name','fc3')
    batchNormalizationLayer('Name','bnorm3')
    leakyReluLayer(0.1,'Name','relu3')

    additionLayer(2,'Name','add3')
    
    fullyConnectedLayer(128,'Name','fc4')
    batchNormalizationLayer('Name','bnorm4')
    leakyReluLayer(0.1,'Name','relu4')

    additionLayer(2,'Name','add4')
    
    fullyConnectedLayer(128,'Name','fc5')
    batchNormalizationLayer('Name','bnorm5')
    leakyReluLayer(0.1,'Name','relu5')

    additionLayer(2,'Name','add5')
    
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(96,'Name','fc6')
    regressionLayer('Name','reg')
];

lgraph = layerGraph(layers); 
lgraph = connectLayers(lgraph,'bnorm1','add1/in2');
lgraph = connectLayers(lgraph,'bnorm2','add2/in2');
lgraph = connectLayers(lgraph,'bnorm3','add3/in2');
lgraph = connectLayers(lgraph,'bnorm4','add4/in2');
lgraph = connectLayers(lgraph,'bnorm5','add5/in2');

elseif strcmp(netType,'fcLarge') % 75,632 parameters
    
numParams = (16*48+48*96+96*128+128*256+256*96+48+96+128+256+96);

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
    
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(96,'Name','fc5')
    regressionLayer('Name','reg')
];

lgraph = layerGraph(layers); 
lgraph = connectLayers(lgraph,'bnorm1','add1/in2');
lgraph = connectLayers(lgraph,'bnorm2','add2/in2');
lgraph = connectLayers(lgraph,'bnorm3','add3/in2');
lgraph = connectLayers(lgraph,'bnorm4','add4/in2');

    
elseif strcmp(netType,'fcMedium') % numParams = 66,560
numParams = 16*32+32*48+48*64+64*96+96*128+128*192+192*96;

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
    
    fullyConnectedLayer(128,'Name','fc5')
    batchNormalizationLayer('Name','bnorm5')
    leakyReluLayer(0.1,'Name','relu5')

    additionLayer(2,'Name','add5')
    
    fullyConnectedLayer(192,'Name','fc6')
    batchNormalizationLayer('Name','bnorm6')
    leakyReluLayer(0.1,'Name','relu6')

    additionLayer(2,'Name','add6')
    
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(96,'Name','fc7')
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

if strcmp(netType,'large') || strcmp(netType,'xlarge')
% CALCULATE NUMBER OF PARAMETERS
nConvLayers = numel(filterSizes(:,1));
numParams = 0;
for k = 1 : nConvLayers
  if k == 1
    numParams = numParams + 2*prod(filterSizes(k,:));    
  else
    numParams = numParams + filterSizes(k-1,3)*prod(filterSizes(k,:));  
  end  
end % end for
end % end if

end % function chooseNetwork
