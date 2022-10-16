clear functions global;
clc;

n = 128;                                 % Number of discretization points

% Setting modeling parameters
prams.T = 6;                            % Simulation time
prams.m = 60;                           % Number of time steps

prams.kappa = 1e0;                      % Bending modulus
prams.viscCont = 1;                     % Viscosity contrast 
prams.rhoIn = 1;                        % Density of the fluid inside
prams.rhoOut = 1;                       % Density of the fluid outside
prams.order = 3;                        % order of the slover
prams.ts = prams.T/prams.m;             % step size (time)
prams.m = prams.T/prams.ts;             % Number of time steps
prams.Incompressibility = 1;            % Incompressibility
prams.vInf = @(X) farFieldVel(X,'parabolic','R',5,'Umax',15);
                                        % The far field velocity

% Get the initial marker point locations -- an ellipse
g = (0:n-1)'*2*pi/n; a = .8; b = 1.8; c =0;
X =  [a*cos(g);c+b*sin(g)];

% Setting Options
options.scheme = 1;
options.GmresTol = 10^-12;              % GMRES slover tolernce
options.usePlot = 1;                    % To plot the intermidate states   
options.AxesHandle = gca;               % The axes for the  plot
options.axisOn = 1;                     % To show or hide plot axis
options.axis = [-3 3 -2 2];             % Force the range of the axes

% Calling update function with all of the above parameters 
[Xfinal status] = Ves2D(X,prams,options);
