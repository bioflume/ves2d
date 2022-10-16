clear all;
set(0,'DefaultAxesFontSize',25)
set(0,'DefaultAxesFontName', 'Helvetica')

addpath ../../src
addpath ../

file = './gtruthsNew/shear2VesRA99VC10Data.bin';  
viscCont = 10;
farField = 'shear';
  
Nup = 96;
plotEffvsTime = ~false;

% Load data
[posx,posy,~,~,~,~,~,time,N,nv] = loadFile(file); 

% Parameters 
skip = 100;
prams.N = N; 
prams.nv = nv;
prams.T = time(end);
prams.m = numel(time);
prams.viscCont = viscCont*ones(nv,1);
prams.kappa = 1e-1;

options.farField = farField;
options.inextens = 'method1';
options.order = 1;
options.vesves = 'implicit';
options.near = true;
options.antiAlias = true;
options.nsdc = 1;
options.orderGL = 2;

[options,prams] = initVes2D(options,prams);
ntime = numel(time);


oc = curve;

% Compute the effective viscosity 
effViscGT = zeros(size(time(1:skip:end)));

op = poten(N);
tt = tstep(options,prams);
idx = 1;
for t = 1 : skip : ntime
    
    sig12 = zeros(Nup,nv);
    % put points together
    X = [interpft(posx(:,:,t),Nup);interpft(posy(:,:,t),Nup)];
    
    % USE computeSigAndEta to compute tension and velocity given config.
    vesicle = capsules([posx(:,:,t);posy(:,:,t)],[],[],prams.kappa,prams.viscCont,true);
    vesicle.setUpRate(op);
    [ten,~,~,u,~] = vesicle.computeSigAndEta(tt,[]);

    % get tangent, curvature and jacobian
    [sa,tan,cur] = oc.diffProp(X);
    % get x and y components of normal and tangent
    nx = tan(Nup+1:2*Nup,:);
    ny = -tan(1:Nup,:);
    tx = tan(1:Nup,:);
    ty = tan(Nup+1:2*Nup,:);
    % Upsample tension
    tenUp = interpft(ten(:,:),Nup);
    % Upsample velocity
    ux = interpft(u(1:N,:),Nup); uy = interpft(u(N+1:2*N,:),Nup);
    for k = 1 : nv 
        term1 = -(prams.kappa*cur(:,k).^2.*nx(:,k).*ny(:,k)+...
            tenUp(:,k).*tx(:,k).*ty(:,k));
        term2 = (viscCont-1)*(ux(:,k).*ny(:,k) + uy(:,k).*nx(:,k));
        sig12(:,k) = term1 + term2;
    end
    effViscGT(idx) = 2*pi/(Nup*nv)*sum(sum(sig12.*sa));
    idx = idx + 1;
end

        
if plotEffvsTime
    figure(1);clf;hold on;
    plot(time(1:skip:end),effViscGT,'Color',[.5,.5,.5],'linewidth',4);
    xlabel('Time')
    ylabel('Intrinsic viscosity')
    axis square
    grid on
end
    
