clear all; clc

% Setting the correct path 
P = path; i = find(pwd==filesep); i = i(end);
subPath = pwd; subPath = [subPath(1:i) 'src' filesep 'base'];
if isempty(strfind(P, subPath)),addpath(subPath);end

% Setting modeling parameters
prams.T = 20;                           % Simulation time
prams.kappa = 1e-1;                     % Bending modulus
prams.rhoOut = 1;                       % Density of the fluid outside
prams.Incompressibility = 1;            % Incompressibility
prams.g = [0;-1];                       % Acceleration of gravity
% Far field velocity parameters
prams.flowType = 'confined';            % Bounded or unbounded
prams.bd = @(ind,m) sampleBd(ind,m,1,'couette');
                                        % Function handle of the fixed
                                        % boundaries
prams.bc = @(x) forcing(x,'couette');   % The function for the boundary condition
                                        % The far field velocity
% Setting Options
options.GmresTol = 10^-14;              % GMRES slover tolernce
options.Lp =  2*pi;                     % Arclength parameter range for the
                                        % quadrature integration (fictitious)
options.usePlot = 0;                    % To plot the intermidate states
options.AxesHandle = [];                % The axes for the  plot
options.showAngle = 0;                  % Showing the inclination angle on the plot
options.saveFig = 0;                    % Saving the plot as *.jpg file
options.saveStride = [];                % Save the plot every other saveStride
options.figNameInit =[];
options.velVec = 0;                     % To show the velocity vector of the
                                        % marker points on the plot
options.axisOn = 0;                     % To show or hide plot axis
options.axis = [];                      % Force the range of the axes
options.track = 0;                      % Color marker points differently to be
                                        % able to track them
options.velProfile = 0;                 % To show the velocity profile on the
                                        % inset of the figure
options.showWaitbar = 0;

% Calling update function with all of the above parameters
for nv = [1 3 4 5]
    prams.viscCont = ones(nv,1);
    prams.rhoIn = ones(nv,1);
    for O = 1:4
        prams.order = O;
        for a = 5:8
            n = 2^a;
            prams.m = n;
            prams.ts = prams.T/prams.m;
            prams.M = [n n];
            prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,prams);

            X = boundary(n,'couette','nv',nv);

            Results = Ves2D(X,prams,options,@monitor);

            fileName =  ['..' filesep 'docs' filesep 'testRuns' filesep ...
                'res' num2str(nv) '_' num2str(O) '_' num2str(n)];
            save(fileName,'Results');
        end
    end
end

