clear; clc;
addpath ../src/
addpath ../examples/
addpath ./shannets/
addpath ./shannets/ves_fft_models/

pathofDocument = fileparts(which('Net_ves_relax_midfat.py'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pathofDocument = fileparts(which('Net_ves_adv_fft.py'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pathofDocument = fileparts(which('ves_fft_mode2.pth'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pe = pyenv('Version', '/Users/gokberk/opt/anaconda3/envs/mattorch/bin/python');

load VF30_70VesIC.mat
dnn = dnnTools;
N = 128;
oc = curve;
% Here X is loaded. 75 vesicles discretized witht 64 points


% Standardize a vesicle
k = 1;
Xinit = [interpft(X(1:end/2,k),128); interpft(X(end/2+1:end,k),128)];
[XinitStand, scaling, rotate, trans, sortIdx] = dnn.standardizationStep(Xinit,oc);

theta = (0:127)'/128*2*pi;
bases = 1/128*exp(1i*theta*(0:127));
bases_rr = real(bases);
bases_ii = imag(bases);
basesConcat = [bases_rr;bases_ii];

load ./shannets/ves_fft_in_param.mat
load ./shannets/ves_fft_out_param.mat

Nves = 1;

% Normalize input
modes = [(0:N/2-1) (-N/2:-1)];
modesInUse = 32;
modeList = find(abs(modes)<=modesInUse);


input_list = [];
cnt = 1;
for imode = modeList
  if imode ~= 1
  input_net = zeros(Nves,2,128);
  x_mean = in_param(imode-1,1);
  x_std = in_param(imode-1,2);
  y_mean = in_param(imode-1,3);
  y_std = in_param(imode-1,4);
  input_net(:,1,:) = (Xinit(1:end/2)-x_mean)/x_std;
  input_net(:,2,:) = (Xinit(end/2+1:end)-y_mean)/y_std;
  input_list{cnt} = py.numpy.array(input_net);
  cnt = cnt + 1;
  end
end % imode

% input_conv = py.numpy.array(input_net);
% [XpredictStand] = pyrunfile("relax_predict.py", "predicted_shape", input_shape=input_conv);
[XpredictStand] = pyrunfile("advect_predict.py", "output_list", input_shape=input_list,num_ves=py.int(Nves));
% [XpredictStand] = pyrunfile("advTest.py", "output_list", input_shape=input_conv,num_ves=py.int(Nves));
Nnet = N;
nv = Nves;
Z11r = zeros(Nnet,Nnet,nv); Z12r = Z11r;
Z21r = Z11r; Z22r = Z11r;

for ij = 1 : numel(modeList)-1  
  
  imode = modeList(ij+1); % mode index
  pred = double(XpredictStand{ij}); % size(pred) = [1 2 256]


  % denormalize output
  real_mean = out_param(imode-1,1);
  real_std = out_param(imode-1,2);
  imag_mean = out_param(imode-1,3);
  imag_std = out_param(imode-1,4);
  
  % first channel is real
  pred(1,1,:) = (pred(1,1,:)*real_std) + real_mean;
  % second channel is imaginary
  pred(1,2,:) = (pred(1,2,:)*imag_std) + imag_mean;

  Z11r(:,imode,1) = pred(1,1,1:end/2);
  Z21r(:,imode,1) = pred(1,1,end/2+1:end);
  Z12r(:,imode,1) = pred(1,2,1:end/2);
  Z22r(:,imode,1) = pred(1,2,end/2+1:end);
end