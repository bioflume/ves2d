clear all;
clc;
addpath ../src/
addpath ../examples/
addpath ./shannets/
addpath ./shannets/ves_fft_models/
pathofDocument = fileparts(which('Net_ves_relax_midfat.py'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pe = pyenv('Version', '/Users/gokberk/opt/anaconda3/envs/mattorch/bin/python');

load('./necessaryMatFiles/n128Dt1e-05RelaxDataSet.mat')

driftNet = zeros(2,nInstances);
driftData = zeros(2,nInstances);

prams.N = 128; % num. points for true solve in DNN scheme
prams.nv = 1; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.dt = 1E-5; % time step size
prams.dtRelax = prams.dt;
prams.Nbd = 0;
prams.nvbd = 0;
prams.interpOrder = 1;
oc = curve;
dnn = dnnToolsSingleVes([],prams);

XnewNet = zeros(size(XstandStore));

%%

for ives = 1 : nInstances
  disp(ives)
  X = XstandStore(:,ives);
  XnewNet(:,ives) = dnn.relaxWTorchNet(X);    

  driftData(1,ives) = mean(XnewStandStore(1:end/2,ives)) - mean(XstandStore(1:end/2,ives));
  driftData(2,ives) = mean(XnewStandStore(end/2+1:end,ives)) - mean(XstandStore(end/2+1:end,ives));

  driftNet(1,ives) = mean(XnewNet(1:end/2,ives)) - mean(XstandStore(1:end/2,ives));
  driftNet(2,ives) = mean(XnewNet(end/2+1:end,ives)) - mean(XstandStore(end/2+1:end,ives));
end

