function Couette2Transport

%% Options
% -------------------------------------------------------------------------
PPoptions.saveData          = 1; % save data on-the-fly
PPoptions.compute_vel       = 0; % compute velocity at the grid points
PPoptions.compute_sigEta    = 0; % compute sigma, eta, RSt (compute velocity on small grid)
PPoptions.interpolate       = 0; % interpolate velocity to finer deltaT
PPoptions.pp_IC             = 1; % compute vesicle IC
% -------------------------------------------------------------------------

%% Information 
% -------------------------------------------------------------------------
info.deltaTnew   = 1e-2; % deltaT for transport problem
info.viscCont    = 1;    % Viscosity Contrast
info.FMM         = true; % FMM on for pp?
info.bgFlow      = 'couette';  %'couette' or 'couette10' or 'couette100' 
info.innWallVel  = 1;    % inner wall radial velocity
                         % 'couette'    = 1
                         % 'couette10'  = 10
                         % 'couette100' = 100
info.frqncy      = 10;

% -------------------------------------------------------------------------

%% File Names
% -------------------------------------------------------------------------
file.VesicleData        = '/scratch/gokberk/Case411/Case411TA_couette.bin';
SmallFile.VelData       = '/scratch/gokberk/Case411/Case411TA_Oscil_Vel8x8_TA.bin';
LargeFile.VelData       = '/scratch/gokberk/Case411/Case411TA_Oscil_Vel256x1024_TA.bin';
file.VelLargeInterpData = '/scratch/gokberk/Case411/Case411TA_Oscil_Vel256x1024_intp.bin';
file.BeforeBounds       = '/scratch/gokberk/Case411/timeVels/file';
file.SigmaData          = '/scratch/gokberk/Case411/Case411TA_Sigma.bin';
file.EtaData            = '/scratch/gokberk/Case411/Case411TA_Eta.bin';
file.RStData            = '/scratch/gokberk/Case411/Case411TA_RSt.bin';
file.uData              = '/scratch/gokberk/Case411/Case411TA_u.bin';
file.IC                 = '/scratch/gokberk/Case411/Case411TA_IC256x1024.bin';

% Log Files
SmallFile.VelLog        = '/scratch/gokberk/Case411/Case411TA_Oscil_Vel8x8_TA.log';
LargeFile.VelLog        = '/scratch/gokberk/Case411/Case411TA_Oscil_Vel256x1024_TA.log';
file.VelLargeInterpLog  = '/scratch/gokberk/Case411/Case411TA_Oscil_Vel256x1024_intp.log';
% -------------------------------------------------------------------------

%% Parameters
% -------------------------------------------------------------------------
% Transport problem parameters
% -------------------------------------------------------------------------
radius = [10; 20]; 

% Small Grid Size
small.N_theta = 8; %    : # of points in [0,2*pi)
small.N_radii = 7;  % +1 : # of points in radial direction
[small.x,small.y,~,~,~,~,~] = generateGrid(small.N_theta,small.N_radii,radius);

% Postprocessing points
% Postprocess velocity only for interior nodes, (N_radii-1)xN_theta
small.xInt = small.x(2:end-1,1:end-1); 
small.yInt = small.y(2:end-1,1:end-1); 
small.Xtar = [small.xInt(:);small.yInt(:)];

% Large Grid Size
large.N_theta =  1024;  %    : # of points in [0,2*pi)
large.N_radii =   255;  % +1 : # of points in radial direction
[large.x,large.y,~,large.theta,~,~,~] = generateGrid(large.N_theta,large.N_radii,radius);

% Postprocessing points
% Postprocess velocity only for interior nodes, (N_radii-1)xN_theta
large.xInt = large.x(2:end-1,1:end-1); 
large.yInt = large.y(2:end-1,1:end-1); 
large.Xtar = [large.xInt(:);large.yInt(:)];


%% Postprocessing
% -------------------------------------------------------------------------

% Compute sigma, eta, RSt and velocity on the small grid
if PPoptions.compute_sigEta
    
    postprocessVel(PPoptions,file,SmallFile,small,info,1);
    % 1: compute sigma,eta and RSt together with velocity
end

% Compute velocity on the large grid
if PPoptions.compute_vel
    
    postprocessVel(PPoptions,file,LargeFile,large,info,0);
    % 0: do not compute sigma,eta and RSt, use the ones computed
end

% Interpolate for finer deltaT
if PPoptions.interpolate
   
    velInterpInTime(file,LargeFile,large,radius,info);
    
end

% Postprocess for IC
if PPoptions.pp_IC 
    
    postprocessIC(file,large,info)
    
end

% -------------------------------------------------------------------------
end
%% End of function

function postprocessVel(PPoptions,file,sizeFile,sizePrams,info,computeSigs)

% Compute Velocity, sigma, eta, RSt parameters
% -------------------------------------------------------------------------
[posx,posy,~,wallx,wally,~,~,time,N,nv] = loadFile(file.VesicleData);

addpath ../../src
addpath ../

prams.m = numel(time); % # of time steps
prams.T = time(end)  ; % time horizon
prams.N = N          ; % # of points per vesicle
prams.nv = nv        ; % # of vesicles
prams.Nbd = size(wallx,1); % # of disc points on the walls
prams.Nbdcoarse = size(wallx,1);
prams.nvbd = size(wallx,2); % # of boundaries (2 walls)
prams.kappa = 1e-1;
prams.frqncy = info.frqncy;
prams.viscCont = info.viscCont*ones(1,nv); % viscosity contrast

options.confined = true;
options.fmm = info.FMM;
options.fmmDLP = info.FMM;
options.near = true;
options.farField = info.bgFlow;
options.farFieldSpeed = @(t) sin(2*pi*prams.frqncy*t/prams.T);
[options,prams] = initVes2D(options,prams);

% set all other options and parameters to default values.  None of them
% play a role at this point
tt = tstep(options,prams);

walls = tt.initialConfined(prams,[wallx;wally],0);

u = zeros(2*N,nv);

% Compute WallGetZone once and save for asked grid size
% -------------------------------------------------------------------------
tracers.N = numel(sizePrams.Xtar)/2;
tracers.nv = 1;
tracers.X = sizePrams.Xtar;
[~,NearW2T] = walls.getZone(tracers,2);

% Start Postprocessing Velocity 
% -------------------------------------------------------------------------
if ~computeSigs
    sigmaT = loadDataFile(file.SigmaData); 
    etaT   = loadDataFile(file.EtaData);  
    RStT   = loadDataFile(file.RStData);
%     uT     = loadDataFile(file.uData);
end

% Initialize Log File
fid = fopen(sizeFile.VelLog,'w');
fprintf(fid,'%s\n',['# of vesicles = ' num2str(nv) ', Grid Size = ', num2str(sizePrams.N_radii+1) 'x' num2str(sizePrams.N_theta) ', nTime = ' num2str(prams.m) ', Th = ' num2str(time(end))]);
fclose(fid);
    
% Initialize Files
if PPoptions.saveData && computeSigs
    
    sigmaSize = [prams.N*prams.nv;prams.m];
    writeData([],sigmaSize,file.SigmaData,'initialize')
    
    EtaSize = [2*prams.Nbd*prams.nvbd;prams.m];
    writeData([],EtaSize,file.EtaData,'initialize')
    
    RStSize = [3*prams.nvbd;prams.m];
    writeData([],RStSize,file.RStData,'initialize');
    
%     uSize = [2*prams.N, prams.nv];
%     writeData([],uSize,file.uData,'initialize');
    
end

if PPoptions.saveData
    velSize = [(sizePrams.N_radii-1)*(sizePrams.N_theta)*2;prams.m];
    writeData([],velSize,sizeFile.VelData,'initialize');   
end
    
for tst = 1 : prams.m
    
    tic 
    
    % Since the boundary is time dependent
    walls = tt.initialConfined(prams,[wallx;wally],time(tst));
    [~,NearW2T] = walls.getZone(tracers,2);
    %--------------------------------------------------------------------
    
    X = [posx(:,:,tst);posy(:,:,tst)];    % get ves positions
    
    if computeSigs
        vesicle = capsules(X,[],[],prams.kappa,prams.viscCont,0);
        [sigma,eta,RSt,~,~] = vesicle.computeSigAndEta(tt,walls,time(tst));
    else 
        sigma = reshape(sigmaT(:,tst),prams.N,prams.nv);
        eta   = reshape(etaT(:,tst),2*prams.Nbd,prams.nvbd);
        RSt   = reshape(RStT(:,tst),3,prams.nvbd);
%         u     = reshape(uT(:,tst),uSize(1),uSize(2));
    end
    
    vel = tt.tracersVel(X,sigma,u,prams.kappa,prams.viscCont,walls,eta,RSt,sizePrams.Xtar,NearW2T,time(tst));
    tim = toc; 


    % Save a log file
    keepDiary(sizeFile.VelLog,prams.m,tst,tim,[]);
  

    if PPoptions.saveData && computeSigs
    
        writeData(sigma,[],file.SigmaData,'write')
    
        writeData(eta,[],file.EtaData,'write')
  
        writeData(RSt,[],file.RStData,'write');
        
%         writeData(u,[],file.uData,'write');
    
    end
    
    if PPoptions.saveData
    
        writeData(vel,[],sizeFile.VelData,'write');
    
    end

end
end % Postprocessing is done

%% Interpolate velocity in time
function velInterpInTime(file,sizeFile,sizePrams,radius,info)
% Load vesicle simulation data
% -------------------------------------------------------------------------
[~,~,~,~,~,~,~,time,~,~] = loadFile(file.VesicleData);
Th = time(end);

% Interpolate velocity parameters
% -------------------------------------------------------------------------
ntime_new  = Th/info.deltaTnew + 1    ; % # of time steps      
time_new   = linspace(0,Th,ntime_new)'; % time discrete
% -------------------------------------------------------------------------

% log file
fid = fopen(file.VelLargeInterpLog,'w');
fprintf(fid,'%s\n', ['Time horizon = ' num2str(Th)]);
fprintf(fid,'%s\n', ['deltaTnew = ' num2str(info.deltaTnew)]);
fprintf(fid,'%s\n', ['Grid = ' num2str(sizePrams.N_radii+1) 'x' num2str(sizePrams.N_theta+1)]);
fclose(fid);

% Load postprocessed velocity
vel = loadDataFile(sizeFile.VelData);

% Start interpolating
NoP = size(vel,1); % # of points to be interp'ed(2* due to x and y)

% Save interpolated velocity to a binary file
% Interpolation in time is done for each point
% but velocity field at each time step has to be saved:
for tst = 1 : ntime_new
  fileNames{tst} = [file.BeforeBounds num2str(tst) '.bin'];
  timeFiles(tst) = fopen(fileNames{tst},'w');
end

for point = 1 : NoP
    tic 
    % Interpolate velocity in time:
    interpedVel = interp1(time,vel(point,:),time_new,'spline');
    % Save the velocity to the specific time step:
    for tst = 1 : ntime_new
      fwrite(timeFiles(tst),interpedVel(tst),'double');
    end
    tim = toc;
    keepDiary(file.VelLargeInterpLog,NoP,point,tim,[]);
end

for tst = 1 : ntime_new
  fclose(timeFiles(tst));
end
%% Enforce boundaries

fidFin = fopen(file.VelLargeInterpData,'w');
fwrite(fidFin,[(sizePrams.N_radii+1)*(sizePrams.N_theta+1)*2;ntime_new],'double');
%innerVel = info.innWallVel * radius(1);

innerVel = @(t) radius(1)*sin(2*pi*info.frqncy*t/Th);

for tst = 1 : ntime_new
    
    % Load the interpolated velocity
    fid = fopen(fileNames{tst},'r');
    Vt = fread(fid,'double');
    fclose(fid);

    Vx = reshape(Vt(1:end/2),sizePrams.N_radii-1,sizePrams.N_theta);
    Vy = reshape(Vt(end/2+1:end),sizePrams.N_radii-1,sizePrams.N_theta);

    Vxt = zeros(sizePrams.N_radii+1,sizePrams.N_theta+1); Vyt = Vxt;

    Vxt(2:end-1,1:sizePrams.N_theta) = Vx;
    Vxt(1,1:sizePrams.N_theta) = -innerVel(time_new(tst))*sin(sizePrams.theta);
    Vxt(:,end) = Vxt(:,1);

    Vyt(2:end-1,1:sizePrams.N_theta) = Vy;
    Vyt(1,1:sizePrams.N_theta) = innerVel(time_new(tst))*cos(sizePrams.theta);
    Vyt(:,end) = Vyt(:,1);
    
    disp([num2str(tst) ' of ' num2str(ntime_new)])
    fwrite(fidFin,[Vxt(:);Vyt(:)],'double');
end
fclose(fidFin);

end % Interpolation is done

%% Postprocess IC
function postprocessIC(file,sizePrams,info)

% Get Vesicle info
[posx,posy,~,wallx,wally,~,~,time,N,nv] = loadFile(file.VesicleData);

% Target Points
x = sizePrams.x(:,1:end-1); y = sizePrams.y(:,1:end-1);
Xtar = [x(:);y(:)];

% Choose time step in vesicle simulation as a IC
tst = 30; 
X = [posx(:,:,tst);posy(:,:,tst)];

addpath ../../src
addpath ../

prams.m = numel(time); % # of time steps
prams.T = time(end)  ; % time horizon
prams.N = N          ; % # of points per vesicle
prams.nv = nv        ; % # of vesicles
prams.Nbd = size(wallx,1); % # of disc points on the walls
prams.Nbdcoarse = size(wallx,1);
prams.nvbd = size(wallx,2); % # of boundaries (2 walls)
prams.kappa = 1e-1;
prams.viscCont = ones(1,nv); % viscosity contrast

options.confined = true;
options.fmm = true;
options.fmmDLP = true;
options.near = true;
options.farField = 'couette';
options.farFieldSpeed = @(t) sin(2*pi*info.frqncy*t/prams.T);
[options,prams] = initVes2D(options,prams);


% set all other options and parameters to default values.  None of them
% play a role at this point
tt = tstep(options,prams);

walls = tt.initialConfined(prams,[wallx;wally]);

u = zeros(2*N,nv);

vesicle = capsules(X,[],[],prams.kappa,prams.viscCont,0);
[sigma,~,~,~] = vesicle.computeSigAndEta(tt,walls);

vesicles = capsules(X,sigma,u,prams.kappa,prams.viscCont,0);
op = poten(N);
[~,id] = op.exactLaplaceDL(vesicles,[ones(N,nv);zeros(N,nv)],[],Xtar,1:nv);
id = round(id(1:(sizePrams.N_radii+1)*(sizePrams.N_theta)));

% To avoid any concentration values except 0&1 due to integration errors:
id(id<0) = 0;
id(id>1) = 1;

C_ic = zeros(sizePrams.N_radii+1,sizePrams.N_theta+1);
C_ic(:,1:sizePrams.N_theta) = reshape(id,sizePrams.N_radii+1,sizePrams.N_theta);
C_ic(:,end) = C_ic(:,1);

fid = fopen(file.IC,'w');
fwrite(fid,[sizePrams.N_radii+1;sizePrams.N_theta+1],'double');
fwrite(fid,C_ic(:),'double');
fclose(fid);
end % postprocessing is done

%% Save logs
function keepDiary(fileName,ntime,tst,tim,vel)
fid = fopen(fileName,'a');
fprintf(fid,'%s\n',['Elapsed Time at tt = ' num2str(tst) ' of ' num2str(ntime) ' = ' num2str(tim) ' sec']);
fclose(fid);

% If there is NaN in vel, warn:
if size(find(isnan(vel)==1)) > 0
    fid = fopen(fileName,'a');
    fprintf(fid,'%s\n',['NaN velocity at time step ' num2str(tst)]);
    fclose(fid);
end
end

%% Load from bin file
function dat = loadDataFile(fileName)
fid = fopen(fileName,'r');
val = fread(fid,'double');
nrow = val(1);
ncol = val(2);
dat = reshape(val(3:end),nrow,ncol);
end

%% Write data to bin file
function writeData(dat,sizeDat,fileName,str)
if strcmp(str,'initialize')
    fid = fopen(fileName,'w');
    fwrite(fid,sizeDat,'double');
    fclose(fid);
elseif strcmp(str,'write')
    fid = fopen(fileName,'a');
    fwrite(fid,dat,'double');
    fclose(fid);
end
end
