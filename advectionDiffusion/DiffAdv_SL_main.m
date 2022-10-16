clear all; clc;
%% ------------------------------------------------------------------------
% PARAMETERS
% -------------------------------------------------------------------------

% Grid Size
% ++++++++++++++++++
params.Ntheta = 512; % + 1 = # of points in azimuthal direction
params.Nradii = 127; % + 1 = # of points in radial direction
% ------------------

% Geometry
% ++++++++++++++++++++++
params.radii = [10; 20]; % Radii of cylinders
% ----------------------

% Time Discretization
% +++++++++++++++++++++++
params.Th = 90; % Time Horizon 
params.deltaT_diff = 0.04; % Time step size for diffusion
params.deltaT_adv = params.deltaT_diff/2; % deltaT_adv = deltaT_diff/2 assumed
% -----------------------

% Diffusion Constant
% ++++++++++++++++++
Case111_diffCoeffs = [5e-3;1e-2;2e-2;4e-2;1e-1;2e-1;4e-1;1]*4.44/5;
% Case111_diffCoeffs = [5e-3;1e-2;2e-2;4e-2;1e-1;2e-1;4e-1;1]*2.83/5;
params.D = Case111_diffCoeffs(1);
% ------------------

%% ------------------------------------------------------------------------
% OPTIONS
% -------------------------------------------------------------------------

% Choose velocity profile:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 1: Analytic Couette with r = [10, 20]
% 2: Both have the same angular velocity (1 rad/s)
% 3: Cylinders with angular velocity of 1 rad/s inside, 0.5 rad/s outside
% 4: Zero Velocity 
% 5: Vesicle Velocity

options.Vel = 5;
% -------------------------------------------------------------------------

% Choose initial condition:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 1: Gaussian bump for r = [10, 20]
% 2: Gaussian bump for r = [1,  2 ]
% 3: IC satisfying homogeneous Neuman BCs
% 4: Bessel function
% 5: Dye IC
% 6: Discontinuous 4 circles for r = [10, 20]
% 7: Vesicle IC
% 8: Layer for r = [10, 20]
% 9: Random IC
%10: Moon IC
%11: Alternative IC (Same interaction surface and VF as 8)

options.IC = 5;
% -------------------------------------------------------------------------


% Interpolation Schemes
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options.conInterpScheme = 'spline';
options.velInterpScheme = 'spline';
% -------------------------------------------------------------------------

% Choose BC type
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% options.BC_type = 'Dirichlet';
options.BC_type = 'Neumann';
% -------------------------------------------------------------------------

% Measure mixing
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options.measure = 1;
% -------------------------------------------------------------------------

% Save data
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options.saveData      = 1;
params.ConDataFile    = ...
'/scratch/gokberk/mixNormTest/DyeIC/Case111_VesVel_Con_Pe1.bin';
params.mixingDataFile = ...
'/scratch/gokberk/mixNormTest/DyeIC/Case111_VesVel_Mix_Pe1.bin';
params.LogFile        = ...
'/scratch/gokberk/mixNormTest/DyeIC/Case111_VesVel_Pe1.log';
% -------------------------------------------------------------------------

% Visualize
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options.visualize = 0;
options.saveFig   = 0;
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% SPACE DISCRETIZATION
% -------------------------------------------------------------------------
[x,y,r1D,theta,~,~,~] = generateGrid(params.Ntheta,params.Nradii,params.radii);

[THETA,R] = meshgrid([theta;2*pi],r1D);

%% ------------------------------------------------------------------------
% TIME DISCRETIZATION
% -------------------------------------------------------------------------
ntime   = params.Th/params.deltaT_diff; % Number of time intervals for diffusion

% Initialize OC Matrix to keep concentration 
% OC = zeros(params.Nradii+1,params.Ntheta+1,ntime+1);



%% ------------------------------------------------------------------------
% FLOW PROPERTIES
% -------------------------------------------------------------------------

if options.Vel == 1

    % 1.1) Stokes flow for Couette Apparatus with r=[10,20]
    V1 = @(x,y)  y/3 .* (1-400./(x.^2+y.^2)); % x-velocity
    V2 = @(x,y) -x/3 .* (1-400./(x.^2+y.^2)); % y-velocity
    
elseif options.Vel == 2

    % 1.2) Both cylinders have the same angular velocity of 1 rad/s 
    V1 = @(x,y) -y; % x-velocity
    V2 = @(x,y)  x; % y-velocity
    
elseif options.Vel == 3

    % 1.3) Angular velocity for inner = 1 rad/s, for outer = 0.5 rad/s
    V1 = @(x,y) -y./sqrt(x.^2+y.^2)*2*pi;
    V2 = @(x,y)  x./sqrt(x.^2+y.^2)*2*pi;
    
elseif options.Vel == 4
    
    % 1.4) Zero velocity
    V1 = @(x,y) zeros(params.Nradii+1,params.Ntheta+1);
    V2 = @(x,y) zeros(params.Nradii+1,params.Ntheta+1);

elseif options.Vel == 5

    % 2) If VESICLE velocity field is used, uncomment this section:
%     fid = fopen('ves40velocity256x1025_p01.bin','r');
%     fid = fopen('Case1Velocity128x512_p01.bin','r');
%     fid = fopen('Case1Velocity128x512_TAdap_p01.bin','r');
%     nrow = (params.Nradii+1)*(params.Ntheta+1)*2;
%     ncol = (params.Th/(params.deltaT_adv/2))+1;
%     vel = fread(fid,[nrow ncol],'double');

    % 256x1024 grid with ntime = 9001 is too large to store in workspace
%     fid = fopen('Case211TA_Vel128x512_p01.bin','r');

    fid = fopen('/scratch/gokberk/Case111/Case111TA_48N_Vel128x512_intp.bin','r');
   % fid = fopen('/scratch/gokberk/VF5_Oscil/AD/Case111TA_Oscil_Vel128x512_intp.bin','r');
    val = fread(fid,2,'double');
%  Read and remove velocity (not necessary for the largest simulation)
   % for i = 1 : 2999
    %    vel = fread(fid,val(1),'double');
    %end
    %clear vel;
%     fclose(fid);
%     vel = reshape(val(3:end),val(1),val(2));
end

if options.Vel ~= 5
    
    % To make it compatible with vesicle velocity field 
    ntime_mid = ntime * params.deltaT_diff/params.deltaT_adv*2; % Number of time intervals for delta_t/2
    Vx        = V1(x,y);
    Vy        = V2(x,y);
    vel       = repmat([Vx(:);Vy(:)],1,ntime_mid+1);
    
end

% ------------------------------------------------------------------------
% FLUID PROPERTIES
% -------------------------------------------------------------------------

if options.IC == 1
    
    % 1.1) Gaussian Bump

    % 1.1.1) If radius = [10,20]
    centers = [15;0]; % Center coordinates of the bumps [x1 x2;y1 y2]
    OC(:,:,1) = exp(-0.1*(((x-centers(1)).^2 + (y-centers(2)).^2)));
    
elseif options.IC == 2
    
    % 1.1.2) If radius = [1,2]
    centers = [1.5;0];
    OC(:,:,1) = exp(-200*(((x-centers(1)).^2 + (y-centers(2)).^2)));

elseif options.IC == 3

    % 1.2) IC which satisfies Neumann BCs
    C = 1;
    A = 12*C;
    B = -7.5*C;
    r2 = x.^2 + y.^2;
    OC(:,:,1) = (A*r2 + B*r2.^2 + C*r2.^3).*cos(THETA);

elseif options.IC == 4

    % 1.3) Bessel Function
    mu0     = 0.6773360049;
    coeff1  = bessely(0,mu0)-bessely(1,mu0)/mu0; coeff2 = -(besselj(0,mu0)-besselj(1,mu0)/mu0);
    C_exact = @(time,theta) exp(-mu0^2*time)*(coeff1*besselj(1,mu0*r1D)+coeff2*bessely(1,mu0*r1D))*cos(theta);
    OC(:,:,1) = C_exact(0,[theta;2*pi]'); % at t = 0

elseif options.IC == 5
    
    % 1.4) Circle with 1s and 0s
    C_ic = zeros(params.Nradii+1,params.Ntheta+1);

    % If radius = [1, 2]
    C_ic(((x-15).^2 + y.^2)<16) = 1;


    OC(:,:,1) = C_ic;
    
elseif options.IC == 6
    
    % 1.4) Circle with 1s and 0s
    % If radius = [10, 20]

    C_ic = zeros(params.Nradii+1,params.Ntheta+1);

    % Equally distributed 4 circles
    C_ic((x.^2 + (y+15).^2)<2) = 1;
    C_ic((x.^2 + (y-15).^2)<2) = 1;
    C_ic(((x-15).^2 + y.^2)<2) = 1;
    C_ic(((x+15).^2 + y.^2)<2) = 1;
    
    OC(:,:,1) = C_ic;

elseif options.IC == 7

    % 2) If VESICLE IC is used: 
%     load IC40_256x1025  
%     load IC_Case1_128x513
    load CouetteSin_IC128x513
    
elseif options.IC == 8
    
    % 7) A layer
    C_ic = zeros(params.Nradii+1,params.Ntheta+1);
    C_ic(R>=14&R<=16) = 1;
%     OC(:,:,1) = C_ic;

elseif options.IC == 9
    
    % 9) Load random IC
    load RandIC_r5
    
elseif options.IC == 10
    
    % 7) Moon
    C_ic = zeros(params.Nradii+1,params.Ntheta+1);
    C_ic(x<=0) = 1;
%     OC(:,:,1) = C_ic;

elseif options.IC == 11
    
    % 8) Alternative IC
    C_ic = zeros(params.Nradii+1,params.Ntheta+1);
    C_ic(R>=12&R<=14& x>=0) = 1;
    C_ic(R>=16&R<=18& x>=0) = 1;
end

% To view the initial condition
% surf(x,y,C_ic,'EdgeColor','none');view(2);shading interp;axis equal;
% pause
%% ------------------------------------------------------------------------
%  SOLVING
% -------------------------------------------------------------------------

% Call the solver
SemiLag_Solver(C_ic,fid,val,params,options);


%% ------------------------------------------------------------------------
%  POSTPROCESSING
% -------------------------------------------------------------------------

if options.visualize 
    postprocessing(OC,params.Th,params.deltaT_diff,x,y,10,options.saveFig)
%     postprocessingMultiple(AnalytLargeDiffOC,AnalytMedDiffOC,AnalytLessDiffOC,LargeDiffOC,MedDiffOC,LessDiffOC,params.Th,params.deltaT_diff,x,y,1,options.saveFig)
end
