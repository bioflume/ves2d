clear all;
set(0,'DefaultAxesFontSize',25)
set(0,'DefaultAxesFontName', 'Helvetica')

addpath ../../src
addpath ../

file = './lowRes/taylorGreen4VesVC1N48NoTricksData.bin';  

if strncmpi(file,'./lowRes/shear2VesRA99VC1',25)
  viscCont = 1;
  farField = 'shear';
  load ./GTeffectViscosities/shear2VesVC1RA99.mat
end
if strncmpi(file,'./lowRes/shear2VesRA99VC10',26)
  viscCont = 10;
  farField = 'shear';
  load ./GTeffectViscosities/shear2VesVC10RA99.mat
end
if strncmpi(file,'./lowRes/shear2VesRA65VC1',25)  
  viscCont = 1;
  farField = 'shear';
  load ./GTeffectViscosities/shear2VesVC1RA65.mat
end
if strncmpi(file,'./lowRes/shear2VesRA65VC10',26)
  viscCont = 10;
  farField = 'shear';
  load ./GTeffectViscosities/shear2VesVC10RA99.mat
end
if strncmpi(file,'./lowRes/taylorGreen4VesVC1',27)  
  viscCont = 1;
  farField = 'taylorGreen';
  load ./GTeffectViscosities/taylorGreen4VesVC1.mat
end
if strncmpi(file,'./lowRes/taylorGreen4VesVC10',28)
  viscCont = 10;
  farField = 'taylorGreen';
  load ./GTeffectViscosities/taylorGreen4VesVC10.mat
end
  
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
prams.gmresTol = 1e-10;

options.farField = farField;
options.inextens = 'method1';
options.order = 1;
options.vesves = 'implicit';
options.near = true;
options.antiAlias = true;
options.nsdc = 0;
options.orderGL = 2;

[options,prams] = initVes2D(options,prams);
ntime = numel(time);


oc = curve;

% Compute the effective viscosity 
effViscLR = zeros(size(time));

op = poten(N);
tt = tstep(options,prams);

for t = 1 : ntime
    
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
    effViscLR(t) = 2*pi/(Nup*nv)*sum(sum(sig12.*sa));
end


% Take time integral
time1 = linspace(0,prams.T,numel(effViscGT));
avgViscGT = trapz(time1,effViscGT);
avgViscLR = trapz(time,effViscLR);
% Compute the error
viscErr = abs(avgViscGT-avgViscLR)/abs(avgViscGT) * 100;

disp('*****************************************************')
disp(['Effective viscosity from Ground Truth: ' num2str(avgViscGT,'%4.2e')])
disp(['Effective viscosity from Low Resolutn: ' num2str(avgViscLR,'%4.2e')])
disp('---------------')
disp(['Error in effective viscosity: ' num2str(viscErr,'%4.2f') '%'])
disp('*****************************************************')
        
if plotEffvsTime
    figure(1);clf;hold on;
    plot(time1,effViscGT,'Color',[.5,.5,.5],'linewidth',4);
    plot(time,effViscLR,'r','linewidth',2)
    xlabel('Time')
    ylabel('Intrinsic viscosity')
    legend('Ground truth','Low resolution')
    legend boxoff
    axis square
    grid on
end
    
