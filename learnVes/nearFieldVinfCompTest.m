clear; clc;

addpath ../src/

oc = curve;
dnn = dnnTools;

Nup = 128;
N = 128;
nlayers = 5;

load ./necessaryMatFiles/X106KinitShapes.mat
k = 1;
Xinit =  [interpft(Xstore(1:end/2,k),Nup); interpft(Xstore(end/2+1:end,k),Nup)];

maxLayerDist = @(h) sqrt(h);

nmodes = 128;
theta = (0:N-1)'/N*2*pi;
ks = (0:N-1)';
basis = 1/N * exp(1i*theta*ks');
activeModes = [(1:nmodes/2)';(N-nmodes/2+1:N)'];

op = poten(Nup);

[Xinit,scaling,rotate,trans,sortIdx] = dnn.standardizationStep(Xinit,oc);
[~,area,len] = oc.geomProp(Xinit);


% Generate the grid for velocity
[~,tang] = oc.diffProp(Xinit);
% get x and y components of normal vector at each point
nx = tang(Nup+1:2*Nup);
ny = -tang(1:Nup);

% Points where velocity is calculated involve the points on vesicle
tracersX = zeros(2*Nup, nlayers-1);
tracersX(:,1) = Xinit;

% initialize vesicle
vesicle = capsules(Xinit, [], [], 1, 1, 1);
vesicle.setUpRate();

% Generate tracers
h = vesicle.length/vesicle.N;  % arc-length spacing
dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist(h);
for il = 2 : nlayers
tracersX(:,il) = [Xinit(1:end/2)+nx*dlayer(il);Xinit(end/2+1:end)+ny*dlayer(il)];
end

tracers.N = Nup;
tracers.nv = nlayers-1;
tracers.X = tracersX(:,2:nlayers);

Vinf = @(x,y) [y;zeros(size(x))];

% SLP 
G = op.stokesSLmatrix(vesicle);

% derivatives
[Ben,Ten,Div] = vesicle.computeDerivs;

SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);

% Get the near zone
[~,NearV2T] = vesicle.getZone(tracers,2);

% CALCULATE VELOCITY ON THE FOURIER MODES FIRST
VelOnGridModesReal = zeros(2*Nup,nlayers-1,nmodes);
VelOnGridModesImag = zeros(2*Nup,nlayers-1,nmodes);
selfVelModesReal = zeros(2*Nup,nmodes);
selfVelModesImag = zeros(2*Nup,nmodes);

M = Ten*((Div*G*Ten)\eye(Nup))*Div;
M11 = M(1:end/2,1:end/2); M12 = M(1:end/2,end/2+1:end);
M21 = M(end/2+1:end,1:end/2); M22 = M(end/2+1:end,end/2+1:end);
B1 = real(basis(:,activeModes)); B2 = imag(basis(:,activeModes));
Z11 = M11*B1+M12*B2; Z12 = M12*B1-M11*B2;
Z21 = M21*B1+M22*B2; Z22 = M22*B1-M21*B2;
Zmat = [Z11 Z12; Z21 Z22];


kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;
SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);

for imode = 1 : nmodes
  BvecR = Zmat(:,imode);
  BvecI = Zmat(:,imode+nmodes);
  
  VelOnGridModesReal(:,:,imode) = op.nearSingInt(vesicle,BvecR,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);
  VelOnGridModesImag(:,:,imode) = op.nearSingInt(vesicle,BvecI,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);
  selfVelModesReal(:,imode) = G*BvecR;
  selfVelModesImag(:,imode) = G*BvecI;
end

% Now reconstruct from the fourier coeff of velocity
% Choose vinf
vinf = [Xinit(end/2+1:end); 0*Xinit(1:end/2)];
 
vz = vinf(1:end/2)+1i*vinf(end/2+1:end);
coeffs = fft(vz);
Vre = real(coeffs); Vim = imag(coeffs);

% Reconstruct the self-velocity
selfVelX_recon = selfVelModesReal*Vre - selfVelModesImag*Vim;
selfVelY_recon = selfVelModesReal*Vim + selfVelModesImag*Vre;

% Actual function eval
selfVelactual = G*M*vinf;

errSelfVel = norm(selfVelactual-[selfVelX_recon;selfVelY_recon]);



