clc;
% Setting the correct path
P = path; i = find(pwd==filesep); i = i(end);
subPath = pwd; subPath = [subPath(1:i) 'src' filesep 'base'];
if isempty(strfind(P, subPath)),addpath(subPath);end

% Setting modeling parameters
prams.kappa = 1e0;
prams.Incompressibility = 1;

% Initial marker points
n = 64;
%g = (0:n-1)'*2*pi/n; a = 1; b = .7; c=.7;
%X =  [a*cos(g);c+b*sin(g)];
X = boundary(n);

% To get the time the method diverges
% global state;

prams.order = 1;
prams.m = 1;
options = [];
options.usePlot = 1;
options.scheme = 1;
options.verbose = 0;

%options.axisOn = 1;
%options.axis = [-5 5 -2 2];

dbstop error
Ts = []; Prams = [];
fprintf('order       shear          T            m             ts \n');
fprintf('------------------------------------------------------------\n');

for order = 1:4
    prams.order = order;
    FLAG = 0;
    for ind= [5 10 15 20 30 40 50]
        prams.ts = 32;
        s = .1*ind;
        if(prams.m<=2^12), FLAG = 0;end
        while(~FLAG)
            clear functions global
            prams.T = 32;
            prams.m = prams.T/prams.ts;
            if(prams.m<=32)
                prams.T = 32*prams.ts;
                prams.m = 32;
            end
            prams.vInf = @(X) farFieldVel(X,'shear',s);
            [Results FLAG]= Ves2D(X,prams,options);
            if(prams.m>2^12), FLAG = 1; end
            prams.ts = prams.ts/1.175;
        end
        Prams{end+1} = prams;
        fprintf('  %1.0f\t\t\t % 3.1f\t\t % 4.0f\t\t % 4.0f\t\t %8.7e\n',...
             prams.order,s,prams.T,prams.m,1.175*prams.ts);
    end
end

% Scheme 1:
% order          shear         T           m              ts 
% -----------------------------------------------------------------
%   1			  0.5		  128		   32		 4.0000000e+000
%   1			  1.0		   32		   32		 1.0000000e+000
%   1			  1.5		   32		  256		 1.2500000e-001
%   1			  2.0		   32		  256		 1.2500000e-001
%   1			  3.0		   32		  1024		 3.1250000e-002
%   1			  4.0		   32		  2048		 1.5625000e-002
%   1			  5.0		   32		  4096		 7.8125000e-003
%   2			  0.5		   32		   64		 5.0000000e-001
%   2			  1.0		   32		  128		 2.5000000e-001
%   2			  1.5		   32		  256		 1.2500000e-001
%   2			  2.0		   32		  256		 1.2500000e-001
%   2			  3.0		   32		  512		 6.2500000e-002
%   2			  4.0		   32		  1024		 3.1250000e-002
%   2			  5.0		   32		  1024		 3.1250000e-002
%   3			  0.5		   32		   32		 1.0000000e+000
%   3			  1.0		   32		  128		 2.5000000e-001
%   3			  1.5		   32		  128		 2.5000000e-001
%   3			  2.0		   32		  256		 1.2500000e-001
%   3			  3.0		   32		  256		 1.2500000e-001
%   3			  4.0		   32		  512		 6.2500000e-002
%   3			  5.0		   32		  1024		 3.1250000e-002
%   4			  0.5		   32		  256		 1.2500000e-001
%   4			  1.0		   32		  512		 6.2500000e-002
%   4			  1.5		   32		  512		 6.2500000e-002
%   4			  2.0		   32		  1024		 3.1250000e-002
%   4			  3.0		   32		  1024		 3.1250000e-002
%   4			  4.0		   32		  2048		 1.5625000e-002
%   4			  5.0		   32		  2048		 1.5625000e-002
% 
% Scheme: 2
% order       shear          T            m             ts 
% ------------------------------------------------------------
%   1			  0.5		   32		   32		 1.0000000e+000
%   1			  1.0		   32		  512		 6.2500000e-002
%   1			  1.5		   32		  512		 6.2500000e-002
%   1			  2.0		   32		  512		 6.2500000e-002
%   1			  3.0		   32		  512		 6.2500000e-002
%   1			  4.0		   32		  512		 6.2500000e-002
%   1			  5.0		   32		  1024		 3.1250000e-002
%   2			  0.5		   32		  8192		 3.9062500e-003
%   2			  1.0		   32		  8192		 6.4000000e+001
%   2			  1.5		   32		  8192		 6.4000000e+001
%   2			  2.0		   32		  8192		 6.4000000e+001
%   2			  3.0		   32		  8192		 6.4000000e+001
%   2			  4.0		   32		  8192		 6.4000000e+001
%   2			  5.0		   32		  8192		 6.4000000e+001
%   3			  0.5		   32		  8192		 3.9062500e-003
%   3			  1.0		   32		  8192		 6.4000000e+001
%   3			  1.5		   32		  8192		 6.4000000e+001
%   3			  2.0		   32		  8192		 6.4000000e+001
%   3			  3.0		   32		  8192		 6.4000000e+001
%   3			  4.0		   32		  8192		 6.4000000e+001
%   3			  5.0		   32		  8192		 6.4000000e+001
%   4			  0.5		   32		  8192		 3.9062500e-003
%   4			  1.0		   32		  8192		 6.4000000e+001
%   4			  1.5		   32		  8192		 6.4000000e+001
%   4			  2.0		   32		  8192		 6.4000000e+001
%   4			  3.0		   32		  8192		 6.4000000e+001
%   4			  4.0		   32		  8192		 6.4000000e+001
%   4			  5.0		   32		  8192		 6.4000000e+001
% 
%
% Scheme 1:
% order       shear          T            m             ts 
% ------------------------------------------------------------
%   1			  0.5		  174		   32		 5.4291701e+000
%   1			  1.0		   41		   32		 1.2717234e+000
%   1			  1.5		   32		   66		 4.8324286e-001
%   1			  2.0		   32		  205		 1.5627891e-001
%   1			  3.0		   32		  633		 5.0540009e-002
%   1			  4.0		   32		  1418		 2.2565553e-002
%   1			  5.0		   32		  2300		 1.3910168e-002
%   2			  0.5		  126		   32		 3.9324002e+000
%   2			  1.0		   48		   32		 1.4942749e+000
%   2			  1.5		   32		   48		 6.6717717e-001
%   2			  2.0		   32		   78		 4.1127052e-001
%   2			  3.0		   32		  205		 1.5627891e-001
%   2			  4.0		   32		  332		 9.6335594e-002
%   2			  5.0		   32		  539		 5.9384510e-002
%   3			  0.5		   32		   35		 9.2112149e-001
%   3			  1.0		   32		   66		 4.8324286e-001
%   3			  1.5		   32		   91		 3.5001746e-001
%   3			  2.0		   32		  126		 2.5352102e-001
%   3			  3.0		   32		  174		 1.8362772e-001
%   3			  4.0		   32		  283		 1.1319432e-001
%   3			  5.0		   32		  390		 8.1987740e-002
%   4			  0.5		   32		  205		 1.5627891e-001
%   4			  1.0		   32		  283		 1.1319432e-001
%   4			  1.5		   32		  390		 8.1987740e-002
%   4			  2.0		   32		  539		 5.9384510e-002
%   4			  3.0		   32		  744		 4.3012774e-002
%   4			  4.0		   32		  1027		 3.1154567e-002
%   4			  5.0		   32		  1418		 2.2565553e-002




