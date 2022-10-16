%clear all
addpath ../src;

density = @(omega) [exp(cos(omega));exp(cos(sin(omega)))];
%density = @(omega) [ones(size(omega));ones(size(omega))];
geomName = 'geom3';

Nup = 64; thetaUp = (0:Nup-1)'*2*pi/Nup;
Xup = geoms(Nup,geomName);
vesicleUp = capsules(Xup,[],[],1,0,false);
denUp = density(thetaUp);

[xtar,ytar] = meshgrid(linspace(-4.0,-3.01,20),linspace(-1.1,1.1,20));
Xtar = [xtar;ytar];
targets = capsules(Xtar,[],[],1,0,false);
[~,NearV2T] = vesicleUp.getZone(targets,3);

op = poten(Nup,vesicleUp.antiAlias);

kernel = @op.exactStokesSL;
SLP = @(X) op.exactStokesSLdiag(vesicleUp,[],X);
Fslp = op.nearSingInt(vesicleUp,denUp,SLP,SLP,...
    NearV2T,kernel,kernel,targets,false);
uExact = Fslp(1:end/2,:);
vExact = Fslp(end/2+1:end,:);

N = 8; theta = (0:N-1)'*2*pi/N;
den = density(theta);
X = geoms(N,geomName);

% No anti-aliasing
vesicle = capsules(X,[],[],1,0,false);
op = poten(N,vesicle.antiAlias);
[~,NearV2T] = vesicle.getZone(targets,3);

kernel = @op.exactStokesSL;
SLP = @(X) op.exactStokesSLdiag(vesicle,[],X);
Fslp = op.nearSingInt(vesicle,den,SLP,SLP,...
    NearV2T,kernel,kernel,targets,false);
f1 = SLP(den);
u1 = Fslp(1:end/2,:);
v1 = Fslp(end/2+1:end,:);


% Use anti-aliasing
vesicle = capsules(X,[],[],1,0,true);
op = poten(N,vesicle.antiAlias);
kernel = @op.exactStokesSL;

SLP = @(X) op.exactStokesSLdiag(vesicle,[],X);
Fslp = op.nearSingInt(vesicle,den,SLP,SLP,...
    NearV2T,kernel,kernel,targets,false);
f2 = SLP(den);
u2 = Fslp(1:end/2,:);
v2 = Fslp(end/2+1:end,:);













%DOUBLE LAYER POTENTIAL
%kernel = @op.exactStokesDL;
%DLP = @(X) 0.5*X + op.exactStokesDLdiag(vesicleUp,[],X);
%Fdlp = op.nearSingInt(vesicleUp,denUp,DLP,DLP,...
%    NearV2T,kernel,kernel,targets,false);
%uDup = Fdlp(1:end/2,:);
%vDup = Fdlp(end/2+1:end,:);
%
%
