clear all;
set(0,'DefaultAxesFontSize',25)
set(0,'DefaultAxesFontName', 'Helvetica')

addpath ../../src
addpath ../

file = './couette75VesN96GTcontData2.bin';
%file = '/scratch/gokberk/00recentCouetteRuns/couette75VesN96GTData.bin';
viscCont = 1;
farField = 'couette';
options.confined = true;

% Log File
logFile = './statistics/couette75VesGTt80_100_stats.log';
fid = fopen(logFile,'w');
fprintf(fid,'%s\n',file);
fclose(fid);



Nup = 96;
plotEffvsTime = ~false;

% Load data
[posx,posy,~,wallx,wally,~,~,time,N,nv] = loadFile(file); 
%time = time(1:10001);

% Parameters 
skip = 160;
prams.N = N; 
prams.nv = nv;
prams.T = time(end);
prams.m = numel(time);
prams.viscCont = viscCont*ones(nv,1);
prams.kappa = 1e-1;

options.fmm = true;
options.fmmDLP = true;
options.farField = farField;
options.inextens = 'method1';
options.order = 1;
options.vesves = 'implicit';
options.near = true;
options.antiAlias = true;
options.nsdc = 0;
options.orderGL = 2;

if options.confined
  prams.nvbd = size(wallx,2);
  prams.Nbd = size(wallx,1);
  prams.Nbdcoarse = prams.Nbd;
end

[options,prams] = initVes2D(options,prams);
ntime = numel(time);

oc = curve;
op = poten(N);
tt = tstep(options,prams);

if options.confined
  [walls,~] = tt.initialConfined(prams,[wallx;wally]);
end

% Viscosity File
%--------------------------------------------------------------------------
viscFile = './statistics/couette75VesGTt80_100_viscosity.bin';
fid = fopen(viscFile,'w');
fwrite(fid,numel(time(1:skip:ntime)),'double');
fclose(fid);
%--------------------------------------------------------------------------

% Compute the effective viscosity 
effVisc = zeros(size(time(1:skip:end)));

idx = 1;


for t = 1 : skip : ntime
    tic
    sig12 = zeros(Nup,nv);
    % put points together
    X = [interpft(posx(:,:,t),Nup);interpft(posy(:,:,t),Nup)];
    
    % USE computeSigAndEta to compute tension and velocity given config.
    vesicle = capsules([posx(:,:,t);posy(:,:,t)],[],[],prams.kappa,prams.viscCont,true);
    vesicle.setUpRate(op);
    
    if options.confined
      [ten,eta,RSt,u,~] = vesicle.computeSigAndEta(tt,walls); 
    else
      [ten,~,~,u,~] = vesicle.computeSigAndEta(tt,[]);
    end

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
        term2NoVC(:,k) = (ux(:,k).*ny(:,k) + uy(:,k).*nx(:,k));
    end
    effVisc(idx) = 2*pi/(Nup*nv)*sum(sum(sig12.*sa));
    effTerm2(idx) = 2*pi/(Nup*nv)*sum(sum(term2NoVC.*sa));
    %write viscosity data
    fid = fopen(viscFile,'a');
    fwrite(fid,[effVisc(idx);effTerm2(idx)],'double');
    fclose(fid);
    
    idx = idx + 1;
    telapsed = toc;
    
    % Keep diary
    fid = fopen(logFile,'a');
    fprintf(fid,'%s\n',['Elapsed Time at tt = ' num2str(t) ' of ' num2str(ntime) ' = ' num2str(telapsed) ' sec']);
    fclose(fid);
end

        
if plotEffvsTime
    figure(1);clf;hold on;
    plot(time(1:skip:end),effVisc,'Color',[.5,.5,.5],'linewidth',4);
    xlabel('Time')
    ylabel('Intrinsic viscosity')
    axis square
    grid on
end
    
