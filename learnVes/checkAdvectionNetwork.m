clear; 
clc;
set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')


addpath ../src/
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

oc = curve;
N = 128;
op = poten(N);

% from shear
load('finalShearXclose.mat')
X = [Xf(1:end/2,1)-mean(Xf(1:end/2,1)); Xf(end/2+1:end,1)-mean(Xf(end/2+1:end,1))];
for it = 1 : 5
  X = oc.redistributeArcLength(X);
end

% Smooth the shape
X = oc.upsThenFilterShape(X,512,16);
vesicle = capsules(X,[],[],1,1,0);

% background velocity
speed = 100;
vinf = @(X) speed*[X(end/2+1:end,:);zeros(size(X(1:end/2,:)))]; 
vback = vinf(X);

% Stokeslet
G = op.stokesSLmatrix(vesicle);
[Ben,Ten,Div] = vesicle.computeDerivs;

M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
MVinf = M*vback;


% 

theta = (0:N-1)'/N*2*pi;
ks = (0:N-1)';
basis = 1/N*exp(1i*theta*ks');

M11 = M(1:end/2,1:end/2); M12 = M(1:end/2,end/2+1:end);
M21 = M(end/2+1:end,1:end/2); M22 = M(end/2+1:end,end/2+1:end);
B1 = real(basis); B2 = imag(basis);
Z11 = M11*B1+M12*B2; Z12 = M12*B1-M11*B2;
Z21 = M21*B1+M22*B2; Z22 = M22*B1-M21*B2;

z = vback(1:end/2) + 1i*vback(end/2+1:end);
zh = fft(z);
V1 = real(zh); V2 = imag(zh);
MVinfRec = [Z11*V1 + Z12*V2; Z21*V1 + Z22*V2];


%%

% load ./shannets/ves_fft_in_param.mat % in_param
% load ./shannets/ves_fft_out_param.mat % out_param
load ./shannets/mergedAdv_NormParams.mat
Nnet = 128;
nv = 1;

modes = [(0:Nnet/2-1) (-Nnet/2:-1)];
modesInUse = 128;
modeList = find(abs(modes)<=modesInUse);

[Xstand,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(X,Nnet);
vesicleSt = capsules(Xstand,[],[],1,1,0);
% Stokeslet
Gst = op.stokesSLmatrix(vesicleSt);
[BenSt,TenSt,DivSt] = vesicleSt.computeDerivs;

Mst = Gst*TenSt*((DivSt*Gst*TenSt)\eye(vesicleSt.N))*DivSt;
M11st = Mst(1:end/2,1:end/2); M12st = Mst(1:end/2,end/2+1:end);
M21st = Mst(end/2+1:end,1:end/2); M22st = Mst(end/2+1:end,end/2+1:end);
B1 = real(basis); B2 = imag(basis);
Z11true = M11st*B1+M12st*B2; Z12true = M12st*B1-M11st*B2;
Z21true = M21st*B1+M22st*B2; Z22true = M22st*B1-M21st*B2;


% Normalize input
input_net = zeros(1,2*127,2*Nnet);  
for ij = 1 : 127
  x_mean = in_param(1,1);
  x_std = in_param(1,2);
  y_mean = in_param(1,3);
  y_std = in_param(1,4);
  input_net(1,2*(ij-1)+1,1:128) = (Xstand(1:end/2)-x_mean)/x_std;
  input_net(1,2*(ij-1)+1,129:256) = (Xstand(end/2+1:end)-y_mean)/y_std;

  rr = real(basis(:,ij+1));
  ii = imag(basis(:,ij+1));

  input_net(1,2*ij,1:128) = rr;
  input_net(1,2*ij,129:256) = ii;
end

input_conv = py.numpy.array(input_net);
[Xpredict] = pyrunfile("advect_predict_merged.py","output_list",input_shape=input_conv,num_ves=py.int(nv));

allmodes_pred = double(Xpredict);

Z11r = zeros(Nnet,Nnet,nv); Z12r = Z11r;
Z21r = Z11r; Z22r = Z11r;

for ij = 1 : 127
  
  imode = modeList(ij+1); % mode index # skipping the first mode
  pred = allmodes_pred(:,2*(ij-1)+1:2*ij,:); % size(pred) = [1 2 256]


  % denormalize output
  real_mean = out_param(ij,1);
  real_std = out_param(ij,2);
  imag_mean = out_param(ij,3);
  imag_std = out_param(ij,4);
  
  % first channel is real
  pred(:,1,:) = (pred(:,1,:)*real_std) + real_mean;
  % second channel is imaginary
  pred(:,2,:) = (pred(:,2,:)*imag_std) + imag_mean;

  Z11r(:,imode,:) = reshape(pred(:,1,1:end/2),[Nnet,1,nv]);
  Z21r(:,imode,:) = reshape(pred(:,1,end/2+1:end),[Nnet,1,nv]);
  Z12r(:,imode,:) = reshape(pred(:,2,1:end/2),[Nnet,1,nv]);
  Z22r(:,imode,:) = reshape(pred(:,2,end/2+1:end),[Nnet,1,nv]);
end


vinfStand = standardize(vback,[0;0],rotate,[0;0],1,sortIdx);
z = vinfStand(1:end/2)+1i*vinfStand(end/2+1:end);


zh = fft(z);
V1 = real(zh); V2 = imag(zh);
% Compute the approximate value of the term M*vinf
MVinfStand = [Z11r*V1+Z12r*V2; Z21r*V1+Z22r*V2];
% Need to destandardize MVinf (take sorting and rotation back)
MVinfNN = zeros(size(MVinfStand));
MVinfNN([sortIdx;sortIdx+Nnet]) = MVinfStand;
MVinfNN = rotationOperator(MVinfNN,-rotate,[0;0]);


save advectionNetCheck MVinf vback MVinfNN X Xstand vinfStand Z11r Z21r Z12r Z22r Z11true Z12true Z21true Z22true

if 0
X = XstandStore(:,120);

addpath ../src/
oc = curve;
N = 128;
op = poten(N);
nmodes = 128;

theta = (0:N-1)'/N*2*pi;
ks = (0:N-1)';
basis = 1/N*exp(1i*theta*ks');
activeModes = [(1:nmodes/2)';(N-nmodes/2+1:N)'];

vesicle = capsules(X,[],[],1,1,0);
[Ben,Ten,Div] = vesicle.computeDerivs;

G = op.stokesSLmatrix(vesicle);

M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
M11 = M(1:end/2,1:end/2); M12 = M(1:end/2,end/2+1:end);
M21 = M(end/2+1:end,1:end/2); M22 = M(end/2+1:end,end/2+1:end);
B1 = real(basis(:,activeModes)); B2 = imag(basis(:,activeModes));
Z11 = M11*B1+M12*B2; Z12 = M12*B1-M11*B2;
Z21 = M21*B1+M22*B2; Z22 = M22*B1-M21*B2;

zReal = [Z11;Z21];
zImag = [Z12;Z22];




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(Xin,Nnet)
oc = curve;
N = numel(Xin)/2;
if Nnet ~= N
  Xin = [interpft(Xin(1:end/2),Nnet);interpft(Xin(end/2+1:end),Nnet)];    
end

% Equally distribute points in arc-length
for iter = 1 : 10
  [Xin,~,~] = oc.redistributeArcLength(Xin);
end


X = Xin;
[trans,rotate,rotCent,scaling,sortIdx] = referenceValues(X);

% Fix misalignment in center and angle due to reparametrization
% X = oc.alignCenterAngle(Xin,X);

% standardize angle, center, scaling and point order

X = standardize(X,trans,rotate,rotCent,scaling,sortIdx);
end % standardizationStep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XrotSort = standardize(X,translation,rotation,rotCent,scaling,sortIdx)
N = numel(sortIdx);

% translate, rotate and scale configuration

Xrotated = rotationOperator(X,rotation,rotCent);   
Xrotated = translateOp(Xrotated,translation);

% now order the points
XrotSort = [Xrotated(sortIdx);Xrotated(sortIdx+N)];

XrotSort = scaling*XrotSort;

end % standardize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = destandardize(XrotSort,translation,rotation,rotCent,scaling,sortIdx)

N = numel(sortIdx);    
    
% scaling back
XrotSort = XrotSort/scaling;

% change ordering back 
X = zeros(size(XrotSort));
X([sortIdx;sortIdx+N]) = XrotSort;

% take translation back
X = translateOp(X,-translation);

% take rotation back
X = rotationOperator(X,-rotation,rotCent);


end % destandardize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [translation,rotation,rotCent,scaling,sortIdx] = referenceValues(Xref)
oc = curve;

N = numel(Xref)/2;

% find translation, rotation and scaling
center = oc.getPhysicalCenterShan(Xref);
V = oc.getPrincAxesGivenCentroid(Xref,center);
% % find rotation angle
w = [0;1]; % y-axis
rotation = atan2(w(2)*V(1)-w(1)*V(2), w(1)*V(1)+w(2)*V(2));


% translation = [-mean(Xref(1:end/2));-mean(Xref(end/2+1:end))];
% rotation = pi/2-oc.getIncAngle2(Xref);
       
% find the ordering of the points
rotCent = center;
Xref = rotationOperator(Xref, rotation, center);
center = oc.getPhysicalCenterShan(Xref);
translation = -center;

Xref = translateOp(Xref, translation);

firstQuad = find(Xref(1:end/2)>=0 & Xref(end/2+1:end)>=0);
theta = atan2(Xref(end/2+1:end),Xref(1:end/2));
[~,idx]= min(theta(firstQuad));
sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];

% amount of scaling
[~,~,length] = oc.geomProp(Xref);
scaling = 1/length;
end % referenceValues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xrot = rotationOperator(X,theta, rotCent)
% Get x-y coordinates
Xrot = zeros(size(X));
x = X(1:end/2)-rotCent(1); y = X(end/2+1:end)-rotCent(2);

% Rotated shape
xrot = (x)*cos(theta) - (y)*sin(theta);
yrot = (x)*sin(theta) + (y)*cos(theta);

Xrot(1:end/2) = xrot+rotCent(1);
Xrot(end/2+1:end) = yrot+rotCent(2);
end % rotationOperator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateOp(X,transXY)
Xnew = zeros(size(X));
Xnew(1:end/2) = X(1:end/2)+transXY(1);
Xnew(end/2+1:end) = X(end/2+1:end)+transXY(2);
end  % translateOp  

