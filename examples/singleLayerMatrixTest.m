% function errors = singleLayerMatrixTest(N)
% 
N = 16;
addpath ../src/
viscosity = 1;
radius = 1;
Npoints = N;

theta = [0:Npoints-1]'/Npoints * 2 * pi;
X = [radius * cos(theta); radius*sin(theta)];

ker = kernels(Npoints);
wall = walls(X, [0 0], zeros(size(X)));
alpertSLP = ker.stokesSLmatrixAlpert(wall,viscosity);

regulSLP = ker.stokesSLmatrixRegulWeightless(wall, viscosity);
H = [wall.sa;wall.sa]*2*pi/wall.N;


alpertSymmSLP = ker.stokesSLmatrixAlpertWeightless(wall,viscosity);
alpertSymmSLP = 0.5 * (alpertSymmSLP+alpertSymmSLP');

normal = wall.normal;
normalH = wall.normal .* H;

alpertRes = alpertSLP * normal;
regulRes = regulSLP * normalH;
symmAlpert = alpertSymmSLP * normalH;

errors(1) = 1/N * sum(alpertRes(1:end/2).^2 + alpertRes(end/2+1:end).^2);
errors(2) = 1/N * sum(regulRes(1:end/2).^2 + regulRes(end/2+1:end).^2);
errors(3) = 1/N * sum(symmAlpert(1:end/2).^2 + symmAlpert(end/2+1:end).^2);
% end
%


