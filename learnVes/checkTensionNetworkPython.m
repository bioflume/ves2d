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
op = poten(N,4,0);

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
% bending force
fBend = vesicle.tracJump(X,zeros(N,1));

% Stokeslet
G = op.stokesSLmatrix(vesicle);
[Ben,Ten,Div] = vesicle.computeDerivs;

LHS = (Div*G*Ten);
LHSinv = LHS\eye(size(LHS));
selfBend = G*fBend; % or G*Ben
selfBend2 = G*(-Ben)*X; % or G*Ben

RHSSelf = -Div*selfBend;
RHSSelf2 = -Div*selfBend2;

RHSVel = -Div*vback;

tenSelf = LHSinv*RHSSelf;
tenSelf2 = LHSinv*RHSSelf2;

tenVel = LHSinv * RHSVel;

% NOW CHECK IF RECONSTRUCTION FROM FOURIER COEFFICIENTS MAKE SENSE
theta = (0:N-1)'/N*2*pi;
ks = (0:N-1)';
basis = 1/N*exp(1i*theta*ks');
nmodes = 128;
activeModes = [(1:nmodes/2)';(N-nmodes/2+1:N)'];
M = ((Div*G*Ten)\eye(vesicle.N))*Div;
M11 = M(:,1:end/2); M12 = M(:,end/2+1:end);
B1 = real(basis(:,activeModes)); B2 = imag(basis(:,activeModes));
Z1true = M11*B1+M12*B2; Z2true = M12*B1-M11*B2;

z = vback(1:end/2) + 1i*vback(end/2+1:end);
zh = fft(z);
V1 = real(zh(activeModes)); V2 = imag(zh(activeModes));
tenVelrecon = -(Z1true*V1 + Z2true*V2);


%% 
load ./shannets/tensionAdv_NormParams.mat
torchTenAdvInNorm = in_param;
torchTenAdvOutNorm = out_param;

Nnet = 128; nv = 1;
modes = [(0:Nnet/2-1) (-Nnet/2:-1)];
modesInUse = 128;
modeList = find(abs(modes)<=modesInUse);

[Xstand,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(X,Nnet);

in_param = torchTenAdvInNorm;
out_param = torchTenAdvOutNorm;

% Normalize input
% Normalize input
input_net = zeros(nv,2*127,Nnet);

for imode = 1 : 127 
  for k = 1 : nv
    input_net(k,2*(imode-1)+1,:) = (Xstand(1:end/2,k)-in_param(1,1))/in_param(1,2);
    input_net(k,2*imode,:) = (Xstand(end/2+1:end,k)-in_param(1,3))/in_param(1,4);
  end
end

input_conv = py.numpy.array(input_net);


[Xpredict] = pyrunfile("tension_advect_allModes_predict.py","output_list",input_shape=input_conv,num_ves=py.int(nv),modesInUse=py.int(modesInUse));
Z1 = zeros(Nnet,Nnet,nv); Z2 = Z1;

pred = double(Xpredict); % size (nv, 2*127, 128)

for ij = 1 : numel(modeList)-1
  
  imode = modeList(ij+1); % mode index # skipping the first mode
  
  % denormalize output
  real_mean = out_param(ij,1);
  real_std = out_param(ij,2);
  imag_mean = out_param(ij,3);
  imag_std = out_param(ij,4);
  
  for k = 1 : nv
    % first channel is real
    Z1(:,imode,k) = (pred(k,2*(ij-1)+1,:)*real_std) + real_mean;
    % second channel is imaginary
    Z2(:,imode,k) = (pred(k,2*ij,:)*imag_std) + imag_mean;
  end
  
end

vinfStand = standardize(vback,[0;0],rotate,[0;0],1,sortIdx);
z = vinfStand(1:end/2)+1i*vinfStand(end/2+1:end);
zh = fft(z);
V1 = real(zh); V2 = imag(zh);

MVinfStand = Z1*V1 + Z2*V2;
MVinf = zeros(size(MVinfStand));
MVinf(sortIdx) = MVinfStand;

%%
x_mean = 2.980232033378272e-11; 
x_std = 0.06010082736611366;
y_mean = -1.0086939616904544e-10; 
y_std = 0.13698545098304749;

% Output normalizing parameters
vx_mean = 327.26141357421875; 
vx_std = 375.0673828125;
Xin = zeros(size(Xstand));
Xin(1:end/2,:) = (Xstand(1:end/2,:)-x_mean)/x_std;
Xin(end/2+1:end,:) = (Xstand(end/2+1:end,:)-y_mean)/y_std;
XinitShape = zeros(nv,2,N);
for k = 1 : nv
XinitShape(k,1,:) = Xin(1:end/2,k)'; 
XinitShape(k,2,:) = Xin(end/2+1:end,k)';
end
XinitConv = py.numpy.array(XinitShape);

[Xpredict] = pyrunfile("self_tension_solve.py","predicted_shape",input_shape=XinitConv);


tenStand = double(Xpredict);
tension = zeros(N,nv);

for k = 1 : nv
  tenOut = tenStand(k,1,:)*vx_std + vx_mean;
  tension(sortIdx(:,k),k) = tenOut/scaling(k)^2;
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
