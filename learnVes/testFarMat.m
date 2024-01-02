%%
clear; clc;

addpath ../src/

N = 128;

thet = (0:N-1)'*2*pi/N;
Xwalls = [1 * cos(-thet); 1*sin(-thet)];

op = poten(N);
walls = capsules(Xwalls,[],[],zeros(1,1),zeros(1,1),1);

% Xtracers(:,1) = [1*cos(thet);1*sin(thet)];
% Xtracers(:,2) = [1.5*cos(thet); 1.5*sin(thet)];
% tracers.X = Xtracers;
% tracers.nv = 2;
% tracers.N = size(Xtracers,1)/2;

% f = rand(2*N,1);
% D = op.stokesDLmatrixFar(walls,tracers);
% matVec = D*f; matVecx = matVec(1:end/2); matVecy = matVec(end/2+1:end);
% stokesDLPmatVec = zeros(2*tracers.N,tracers.nv);
% for k = 1 : tracers.nv
%   stokesDLPmatVec(1:end/2,k) = matVecx((k-1)*tracers.N+1:k*tracers.N);
%   stokesDLPmatVec(end/2+1:end,k) = matVecy((k-1)*tracers.N+1:k*tracers.N);
% end
% 
% 
% [~,stokesDLPmatFree] = op.exactStokesDL(walls,f,[],Xtracers,[1]);
% 


%% precompute near field
nLayers = 3;
oc = curve;

[~,tang] = oc.diffProp(Xwalls);
nx = -tang(walls.N+1:2*walls.N);
ny = tang(1:walls.N);
h = walls.length/walls.N;

Ntra = walls.N*(nLayers-1);
tracersXwall= zeros(2*Ntra,1);
dlayer = (0:nLayers-1)'/(nLayers-1)' * sqrt(h);
dlayer = dlayer(2:end);
for il = 1 : nLayers-1
  tracersXwall((il-1)*walls.N+1:il*walls.N,1) = Xwalls(1:end/2)+nx*dlayer(il);
  tracersXwall(Ntra+(il-1)*walls.N+1:Ntra+il*walls.N,1) = Xwalls(end/2+1:end) + ny*dlayer(il);
end
wallDLP = op.stokesDLmatrix(walls);
matOnSelf = -1/2*eye(2*walls.N) + wallDLP;
matOnLayers = op.stokesDLmatrixFar(walls,tracersXwall);

% Given density calculate velocity at the layers
density = [cos(thet); sin(thet)];

xxInput = [Xwalls(1:end/2);tracersXwall(1:end/2)];
yyInput = [Xwalls(end/2+1:end);tracersXwall(end/2+1:end)];
selfDensity = matOnSelf*density;
selfDensityX = selfDensity(1:end/2);
selfDensityY = selfDensity(end/2+1:end);
layersDensity = matOnLayers * density;
layersDensityX = layersDensity(1:end/2);
layersDensityY = layersDensity(end/2+1:end);


checkThet = 0.2;
checkRad = linspace(1.001,1.4,100)';

velPointsX = checkRad*cos(checkThet);
velPointsY = checkRad*sin(checkThet);

% Find near points for interpolation, use direct evaluation for far points
dr = sqrt(velPointsX.^2 + velPointsY.^2)-1;
nearIDs = find(dr<=sqrt(h));
farIDs = find(dr>sqrt(h));
tMAT = tic;
Fx = scatteredInterpolant(xxInput, yyInput, [selfDensityX;layersDensityX],'natural');
Fy = scatteredInterpolant(xxInput, yyInput, [selfDensityY;layersDensityY],'natural');

interpVelX = Fx(velPointsX(nearIDs), velPointsY(nearIDs));
interpVelY = Fy(velPointsX(nearIDs), velPointsY(nearIDs));
timeMatInterp = toc(tMAT);

[~,farDirectVel] = op.exactStokesDL(walls,density,[],[velPointsX(farIDs);velPointsY(farIDs)],[1]);
farDirectVelX = farDirectVel(1:end/2);
farDirectVelY = farDirectVel(end/2+1:end);

interpMethVelX = zeros(size(velPointsX));
interpMethVelY = zeros(size(velPointsY));
interpMethVelX(nearIDs) = interpVelX;
interpMethVelY(nearIDs) = interpVelY;
interpMethVelX(farIDs) = farDirectVelX;
interpMethVelY(farIDs) = farDirectVelY;

% 
% figure(1); clf;
% plot(Xwalls(1:end/2),Xwalls(end/2+1:end),'k','linewidth',2)
% hold on
% plot(xxInput, yyInput, 'k.','markersize',10)
% plot(velPointsX,velPointsY,'rs','markersize',6)
% axis equal

figure(2); clf;
plot(checkRad,interpMethVelX, 'r', 'linewidth',2)
hold on

figure(3); clf;
plot(checkRad, interpMethVelY, 'r','linewidth',2)
hold on

% Compare with near-singular integration
tracers.X = [velPointsX; velPointsY];
tracers.nv = 1;
tracers.N = numel(velPointsX);


%% Near singular integration

% Get the near structure
[~,NearW2T] = walls.getZone(tracers,2);
kernel = @op.exactStokesDL;

jump = -1/2;
DLP = wallDLP + jump*eye(2*walls.N);
DLPfun = @(X) op.exactStokesDLdiag(walls, DLP, X);
velTracers = op.nearSingInt(walls,density,DLPfun,[],NearW2T,kernel,kernel,tracers,false,false);

figure(2); 
plot(checkRad,velTracers(1:end/2), 'k--', 'linewidth',2)

figure(3); 
plot(checkRad, velTracers(end/2+1:end), 'k--','linewidth',2)


errs(1) = mean(sqrt((interpMethVelX-velTracers(1:end/2)).^2 + (interpMethVelY-velTracers(end/2+1:end)).^2)./sqrt(velTracers(1:end/2).^2+velTracers(end/2+1:end).^2));
%% Compare with regularized calculation
if 0
DLPnoCorr = op.stokesDLmatrixNoCorr(walls);
[~,nearField,velIgnoreNear] = op.divideNearFarSLP(walls,...
        density,DLPnoCorr,NearW2T,kernel,kernel,tracers,false);

figure(2); 
plot(checkRad,velIgnoreNear(1:end/2), 'm--', 'linewidth',2)

figure(3); 
plot(checkRad, velIgnoreNear(end/2+1:end), 'm--','linewidth',2)
end

if 0
[~,inaccNearVel] = op.exactStokesDL(walls,density,[],tracers.X,1);

figure(2); 
plot(checkRad,inaccNearVel(1:end/2), 'b-.', 'linewidth',2)
legend('Interpolation','True','no regul')
legend boxoff

figure(3); 
plot(checkRad, inaccNearVel(end/2+1:end), 'b-.','linewidth',2)
legend('Interpolation','True','no regul')
legend boxoff
end

%% Alternative interpolation
tRBF = tic;
opX = rbfcreate([xxInput';yyInput'],[selfDensityX;layersDensityX]','RBFFunction','cubic');
opY = rbfcreate([xxInput';yyInput'],[selfDensityY;layersDensityY]','RBFFunction','cubic');
rbfVelX = rbfinterp([velPointsX(nearIDs)'; velPointsY(nearIDs)'], opX);
rbfVelY = rbfinterp([velPointsX(nearIDs)'; velPointsY(nearIDs)'], opY);
timeRBF = toc(tRBF);

RBFinterpMethVelX = zeros(size(velPointsX));
RBFinterpMethVelY = zeros(size(velPointsY));
RBFinterpMethVelX(nearIDs) = rbfVelX;
RBFinterpMethVelY(nearIDs) = rbfVelY;
RBFinterpMethVelX(farIDs) = farDirectVelX;
RBFinterpMethVelY(farIDs) = farDirectVelY;

figure(2);
plot(checkRad,RBFinterpMethVelX, 'g', 'linewidth',2)
hold on

figure(3); 
plot(checkRad, RBFinterpMethVelY, 'g','linewidth',2)
hold on
errs(2) = mean(sqrt((RBFinterpMethVelX-velTracers(1:end/2)).^2 + (RBFinterpMethVelY-velTracers(end/2+1:end)).^2)./sqrt(velTracers(1:end/2).^2+velTracers(end/2+1:end).^2));
