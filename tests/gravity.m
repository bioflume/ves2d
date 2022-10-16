clear all;clc
% Setting the correct path
P = path; i = find(pwd==filesep); i = i(end);
subPath = pwd; subPath = [subPath(1:i) 'src' filesep 'base'];
if isempty(strfind(P, subPath)),addpath(subPath);end

% Setting Options
options.GmresTol = 10^-14;              % GMRES slover tolernce
options.Lp =  2*pi;                     % Arclength parameter range for the
% quadrature integration (fictitious)
options.usePlot = 1;                    % To plot the intermidate states
options.AxesHandle = gca;               % The axes for the  plot
options.showAngle = 0;                  % Showing the inclination angle on the plot
options.saveFig = 0;                    % Saving the plot as *.jpg file
options.velVec = 0;                     % To show the velocity vector of the
% marker points on the plot
options.axisOn = 1;                     % To show or hide plot axis
options.axis = [-3 3 -3 3];             % Force the range of the axes
options.track = 0;                      % Color marker points differently to be
% able to track them
options.velProfile = 0 ;                % To show the velocity profile on the
% inset of the figure

% Choosing the number of vesicles
nv = 1;                                 % Number of vesicles
% Setting modeling parameters
prams.T = 8;                            % Simulation time
prams.kappa = 1e-1;                     % Bending modulus
prams.viscCont = ones(nv,1);            % Viscosity contrast
prams.rhoIn = 2*ones(nv,1);             % Density of the fluid inside
prams.rhoOut = 1;                       % Density of the fluid outside
prams.order = 0;                        % order of the slover
prams.Incompressibility = 1;            % Incompressibility
prams.vInf = @(X) farFieldVel(X,'shear',0);
% The far field velocity
prams.g = [0;-1];                       % Acceleration of gravity

% Get the initial marker point location
fprintf('  order \t  n  \t  |X-X0|/|X0|  \t\t  |A-A0|/A0   \t\t   |L-L0|/L0\n');
fprintf('----------------------------------------------------------------------------\n');
fd = [' %3.0f   \t %3.0f \t' repmat('  %2.4e\t\t',1,3) '\n'];

for order = 1:4
    prams.order = order;
    for a = [8 5 6 7]
        prams.n = 2^a;                   % Number of discretization points
                 prams.m = prams.n;               % Number of time steps
                 prams.ts = prams.T/prams.m;      % step size (time)
                 X = boundary(prams.n,'nv',nv);
        
                 % Calling update function with all of the above parameters
                 name = genvarname(['Results_' num2str(prams.order) '_' num2str(prams.n)],who);
                 eval([name '= update(X,prams,options);']);
        if(a~=8)
            X0 = eval(['Results_' num2str(order) '_256{end}.X;']);
            [t A0 L0] = reducedVolume(X0);

            eval(['X = Results_' num2str(order) '_' num2str(prams.n) '{end}.X;']); 
            X = reshape(X,[],2);
            Xref = X0(1:256/prams.n:end); Xref = reshape(Xref,[],2);

            e = max(sqrt(dot(X-Xref,X-Xref,2)))/max(sqrt(dot(Xref,Xref,2)));
            [t A L] = reducedVolume(X(:));
            fprintf(fd,order,prams.n,e,abs(A-A0)/A0,abs(L-L0)/L0);
        end
    end
end
save('gravity', '-regexp', 'Res*');





