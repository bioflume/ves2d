clear functions; clc
%dbstop in Ves2D at 103

% Choosing the number of vesicles
n = 64;                                 % Number of vesicles
nv = 3;                                 % Number of discretization points
% Setting modeling parameters
prams.T = 10;                           % Simulation time
prams.m = 400;                          % Number of time steps
prams.kappa = 1e-1;                     % Bending modulus
prams.viscCont = ones(nv,1);            % Viscosity contrast 
prams.rhoIn = ones(nv,1);               % Density of the fluid inside
prams.rhoOut = 1;                       % Density of the fluid outside
prams.order = 2;                        % order of the slover
prams.ts = prams.T/prams.m;             % step size (time)
prams.Incompressibility = 1;            % Incompressibility
prams.g = [0;-1];                       % Acceleration of gravity
% Far field velocity parameters
prams.flowType = 'confined';            % Bounded or unbounded
prams.M = [128 64];                     % Number of discretization points
prams.bd = @(ind,m) sampleBd(ind,m,1,'couette');
                                        % Function handle of the fixed
                                        % boundaries
prams.bc = @(x) forcing(x,'couette');   % The function for the boundary condition 
prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,prams);
                                        % The far field velocity
% Getting the initial marker point location
X = boundary(n,'couette','nv',nv);
% Setting Options
options.GmresTol = 10^-14;              % GMRES slover tolernce
options.Lp =  2*pi;                     % Arclength parameter range for the
                                        % quadrature integration (fictitious)
options.usePlot = 1;                    % To plot the intermidate states
options.progressBar = 1;
options.AxesHandle = gca;               % The axes for the  plot
options.showAngle = 0;                  % Showing the inclination angle on the plot
options.saveFig = 0;                    % Saving the plot as *.jpg file
options.saveStride = [];                % Save the plot every other saveStride
options.figNameInit =[]; 
options.axisOn = 0;                     % To show or hide plot axis
options.axis = [];                      % Force the range of the axes
options.track = 0;                      % Color marker points differently to be


% Calling update function with all of the above parameters 
[XF status] = Ves2D(X,prams,options,@monitor);