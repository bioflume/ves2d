function SemiLag_Solver(C_ic,fidVel,val,params,options)
% -------------------------------------------------------------------------
%% LOAD PARAMETERS and OPTIONS
% -------------------------------------------------------------------------
% Parameters
N_theta = params.Ntheta; N_radii = params.Nradii; radius = params.radii;
Th = params.Th; deltaT_diff = params.deltaT_diff; deltaT_adv = params.deltaT_adv;
D = params.D;

% Options 
velInterpScheme = options.velInterpScheme;
conInterpScheme = options.conInterpScheme;

% Boundary Condition Type
if strcmp(options.BC_type,'Neumann')
    BC_ND = 1;
end
if strcmp(options.BC_type,'Dirichlet')
    BC_ND = 2;
end

% Error tolerance for interpolations
etol = 0;

% -------------------------------------------------------------------------
%% SPACE DISCRETIZATION
% -------------------------------------------------------------------------
M = floor((N_theta+1)/2); % number of Fourier coefficients for each (cos and sin) (~k)

[x,y,r1D,theta,s,R1D,l] = generateGrid(N_theta,N_radii,radius);

% Uniform Mesh for SL advection
% [thetg,sg] = meshgrid([theta;2*pi],s);   % Uniform theta-s space
thetaEx = [theta(end-9:end)-2*pi ;theta; theta(1:10)+2*pi]; % extend theta due to periodicity

[thetg,sg] = meshgrid(thetaEx,s); % for new Cheb grid (s does not contain 0 and pi)

% ------------------------------------------------------------------------
%% TIME DISCRETIZATION
% -------------------------------------------------------------------------
ntime     = Th/deltaT_diff;         % Number of time intervals for diffusion
ntime_adv = deltaT_diff/deltaT_adv; % Number of time intervals for advection
I = 1/deltaT_diff;

% Initialize mixing measure 
if options.measure
    H1Norm = zeros(ntime+1,1);
end

% ---------------------------------------------------------------------
%% Setting the stage for the calculations
% ---------------------------------------------------------------------
% Trig. functions
bcos  = @(l,s) cos(s*l');
bsin  = @(l,s) sin(s*l');

% Derivatives of "s" w.r.t. Chebyshev points (R)
dsdR = @(R) 1./sqrt(1-R.^2);
d2sdR2 = @(R) R./((1-R.^2).^(3/2));

% Derivative of R w.r.t. exact geometry r (due to linear mapping)
% dRdr = 2/(radius(2)-radius(1));
dRdr = (R1D(end)-R1D(1))/(radius(2)-radius(1));

% ---------------------------------------------------------------------
%% Forming the system of equations
% ---------------------------------------------------------------------
% Multipliers of Chebyshev Coefficients (b^c_kl, b_ol) (LHS) = (I-0.5*D*Laplacian)
BC_Multip = @(l,s,dsdR,d2sdR2,r,sgn) sgn*0.5*D*(repmat(dsdR.^2,1,length(l)).*bcos(l,s).*repmat(l'.^2,length(r),1)*dRdr^2 ...
                            + repmat(d2sdR2,1,length(l)).*bsin(l,s).*repmat(l',length(r),1)*dRdr^2 ...
                            + repmat(dsdR,1,length(l)).*repmat(1./r,1,length(l)).*bsin(l,s).*repmat(l',length(r),1)*dRdr) ...
                            + bcos(l,s)*I;

% Compute LHS matrices (excluding BC nodes, i.e. r1 and r2) for Cheb. pts.
% Therefore, constant sections of the matrices are computed once
% (LHS_S = LHS_C, thus, is not created)

LHS_C = zeros(N_radii+1,N_radii+1);
RHS_C = zeros(N_radii+1,N_radii+1);

LHS_C(2:end-1,:) = BC_Multip(l,s(2:end-1),dsdR(R1D(2:end-1)),d2sdR2(R1D(2:end-1)),r1D(2:end-1), 1);
RHS_C(2:end-1,:) = BC_Multip(l,s(2:end-1),dsdR(R1D(2:end-1)),d2sdR2(R1D(2:end-1)),r1D(2:end-1),-1);


% Boundary Nodes
if BC_ND == 1
    % New Cheb Grid (Neumann BCs)
    LHS_C(1,:) = (sin(s(1)*l)).*l;
    LHS_C(N_radii+1,:) = (sin(s(end)*l)).*l;
    RHS_C(1,:) = LHS_C(1,:);
    RHS_C(end,:) = LHS_C(end,:);
else
    % Dirichlet BCs
    LHS_C(1,:) = cos(s(1)*l);
    LHS_C(end,:) = cos(s(end)*l);
    RHS_C(1,:) = LHS_C(1,:);
    RHS_C(end,:) = LHS_C(end,:);
end

w1 = 1/sqrt(N_radii+1);
wL = sqrt(2/(N_radii+1));

LHS_C(:,1) = w1 * LHS_C(:,1); RHS_C(:,1) = w1 * RHS_C(:,1); 
LHS_C(:,2:end) = wL * LHS_C(:,2:end); RHS_C(:,2:end) = wL * RHS_C(:,2:end); 

% Fourier Mode Matrix
FourModeMatrix = repmat(0.5*D./(r1D(2:end-1).^2),1,length(l)).*bcos(l,s(2:end-1)); 
FourModeMatrix(:,1) = w1 * FourModeMatrix(:,1); 
FourModeMatrix(:,2:end) = wL * FourModeMatrix(:,2:end);

% ---------------------------------------------------------------------
%% Initialize Data Files
% ---------------------------------------------------------------------
if options.saveData
    fid = fopen(params.ConDataFile,'w');
    fwrite(fid,[N_radii+1;N_theta+1;ntime+1],'double');
    fwrite(fid,Th,'double');
    fwrite(fid,radius,'double');
    fclose(fid);

    fid = fopen(params.LogFile,'w');
    dateTime = datestr(clock);
    fprintf(fid,'%s\n',['Advection-Diffusion solver is started on: ' dateTime]);
    fprintf(fid,'%s\n',['Initial condition: ' num2str(options.IC) '. Velocity Profile: ' num2str(options.Vel)]);
    fprintf(fid,'%s\n',['Boundary Condition: ' options.BC_type]);
    fprintf(fid,'%s\n',['Nradial x Ntheta = ' num2str(N_radii+1) 'x' num2str(N_theta) '. radius = [' num2str(radius(1)) ', ' num2str(radius(2)) '].']);
    fprintf(fid,'%s\n',['Time horizon = ' num2str(Th) '. Time step for diffusion = ' num2str(deltaT_diff) '. Time step for advection = ' num2str(deltaT_adv)]);
    fclose(fid);

    if options.measure
        fid = fopen(params.mixingDataFile,'w');
        fwrite(fid,numel(H1Norm),'double');
        fwrite(fid,Th,'double');
        fclose(fid);
    end
end
% ---------------------------------------------------------------------

%% ########################################################################
% SOLVING w/ STRANG SPLITTING + BACKWARD EULER for SMOOTHING
% #########################################################################

% Initialize the Chebyshev Coeff Matrices for COS and SIN
B_C = zeros(N_radii+1,M-1);
B_S = zeros(N_radii+1,M-1);

% -------------------------------------------------------------------------
% INITIALLY DO FEW BACKWARD EULER TO LET IC DIFFUSE and SMOOTH OUT
% -------------------------------------------------------------------------

% disp('Few Backward Euler steps for smoothing')
% disp('--------------------------------------------------------------')
% tic % Start timing
% 
delT_BE = deltaT_diff;
BE_Lap = @(l,s,dsdR,d2sdR2,r) delT_BE*D*(repmat(dsdR.^2,1,length(l)).*bcos(l,s).*repmat(l'.^2,length(r),1)*dRdr^2 ...
                            + repmat(d2sdR2,1,length(l)).*bsin(l,s).*repmat(l',length(r),1)*dRdr^2 ...
                            + repmat(dsdR,1,length(l)).*repmat(1./r,1,length(l)).*bsin(l,s).*repmat(l',length(r),1)*dRdr) ...
                            + bcos(l,s);
BE_FourModeMatrix = repmat(delT_BE*D./(r1D(2:end-1).^2),1,length(l)).*bcos(l,s(2:end-1));  
BE_FourModeMatrix(:,1) = w1 * BE_FourModeMatrix(:,1); 
BE_FourModeMatrix(:,2:end) = wL * BE_FourModeMatrix(:,2:end);

% Compute LHS matrices (excluding BC nodes, i.e. r1 and r2) for Cheb. pts.
% Therefore, constant sections of the matrices are computed once

BE_LHS = zeros(N_radii+1,N_radii+1);
BE_LHS(2:end-1,:) = BE_Lap(l,s(2:end-1),dsdR(R1D(2:end-1)),d2sdR2(R1D(2:end-1)),r1D(2:end-1));

% BE_RHS = bcos(l,s);

% Boundary Nodes
if BC_ND == 1
    % New Cheb Grid (Neumann BCs)
    BE_LHS(1,:) = (sin(s(1)*l)).*l;
    BE_LHS(N_radii+1,:) = (sin(s(end)*l)).*l;
else
    % Dirichlet BCs
    BE_LHS(1,:) = cos(s(1)*l);
    BE_LHS(end,:) = cos(s(end)*l);
end

BE_LHS(:,1) = w1 * BE_LHS(:,1); 
BE_LHS(:,2:end) = wL * BE_LHS(:,2:end);
% 
% ntime_smooth = ceil(2/deltaT_diff);
% BE_OC = zeros(N_radii+1,N_theta+1,ntime_smooth+1);
BE_OC(:,:,1) = C_ic;

% for tt = 1 : ntime_smooth
%     C_ic = BE_OC(:,:,tt);
%     C_ic = BE_Smoothing(N_radii,N_theta,M,BE_LHS,BE_FourModeMatrix,B_C,B_S,C_ic);
%     BE_OC(:,:,tt+1) = C_ic;
% end
C_ic = BE_OC(:,:,end);

% tim = toc;  %report timing for single loop
% disp(['Elapsed time is ' num2str(tim) 'sec'])
% disp('--------------------------------------------------------------')

% Save Data
if options.saveData
%    keepDiary(params.LogFile,0,0,tim,0);
   writeData(params.ConDataFile,C_ic);
end
% -------------------------------------------------------------------------

% # of steps of velocity to read at every tt_diff
nsteps_read = deltaT_diff/deltaT_adv * 2;

for tt = 1 : ntime
    
    disp([num2str(tt*deltaT_diff) ' unit of ' num2str(Th) ' :'])
    disp('--------------------------------------------------------------')
        
    % Read velocity at every time step (deltaT_diff/deltaT_adv/2 steps)
    vel = fread(fidVel,nsteps_read*val(1),'double');
    vel = reshape(vel,val(1),nsteps_read);
    
    tic % time for single loop
    
    
    % ---------------------------------------------------------------------
    % 2nd-ORDER STRANG SPLITTING
    % ---------------------------------------------------------------------
    disp('Step 1: ADVECTION')    
    for tt_adv = 1 : ntime_adv/2 
        C_ic = SL_advection(N_radii,N_theta,tt_adv,deltaT_adv,vel,x,y,thetg,sg,radius,velInterpScheme,conInterpScheme,etol,R1D,C_ic);
    end % end first advection
 
    % Using C(tt+deltaT/2) as an IC, solve the diffusion problem:	
    if D ~= 0
        disp('Step 2: DIFFUSION')
        C_ic = BE_Smoothing(N_radii,N_theta,M,BE_LHS,BE_FourModeMatrix,B_C,B_S,C_ic);
%         C_ic = CN_diffusion(N_radii,N_theta,M,LHS_C,RHS_C,FourModeMatrix,B_C,B_S,C_ic);
    end
    
   disp('Step 3: ADVECTION') 
   for tt_adv = ntime_adv/2+1 : ntime_adv
       C_ic = SL_advection(N_radii,N_theta,tt_adv,deltaT_adv,vel,x,y,thetg,sg,radius,velInterpScheme,conInterpScheme,etol,R1D,C_ic);
   end % end 2nd advection
   % ---------------------------------------------------------------------
   
   % ---------------------------------------------------------------------
   % Mixing measurement
   % ---------------------------------------------------------------------
   if (tt == 1 && options.measure)
      H1Norm(1) = MeasureMixingLaplacian(BE_OC(:,1:end-1,end),N_radii,N_theta,radius,BC_ND);
      if options.saveData
          writeData(params.mixingDataFile,H1Norm(1));
      end
   end

   if options.measure
      H1Norm(tt+1) = MeasureMixingLaplacian(C_ic(:,1:end-1),N_radii,N_theta,radius,BC_ND);
      disp(['Normalized Mixing Norm: ' num2str(H1Norm(tt+1)/H1Norm(1))])
   end
   % ---------------------------------------------------------------------
   
   % Check if it blows up
   blowUp = max(max(abs(C_ic))) > 10^3;
   
   tim = toc;  %report timing for single loop
   disp(['Elapsed time is ' num2str(tim) 'sec'])
   if blowUp
       disp('Blowing up')
   end
   disp('--------------------------------------------------------------')
   
   % Save Data
   if options.saveData
       keepDiary(params.LogFile,ntime,tt,tim,blowUp);
       writeData(params.ConDataFile,C_ic);
       if options.measure
           writeData(params.mixingDataFile,H1Norm(tt+1));
       end
   end
   
end % end of time loop

% Close velocity file
fclose(fidVel);

% Conclude log file
if options.saveData
    fid = fopen(params.LogFile,'a');
    dateTime = datestr(clock);
    fprintf(fid,'%s\n',['Advection-Diffusion solver is completed on: ' dateTime]);
    fclose(fid);
end

end%% END OF FUNCTION

% #########################################################################
%% BACKWARD EULER SMOOTHING
% #########################################################################
function C_ic = BE_Smoothing(N_radii,N_theta,M,LHS_C,FourModeMatrix,B_C,B_S,C_ic)
% IC for diffusion problem
CicDiff = C_ic(:,1:end-1);

% Interpolate for spectral coefficients
[~, ~, ~, ~, a0, ac, as, aM] = interpFFT(N_theta, N_radii, CicDiff);

for kk = 1 : M-1
    
    % For ac and as:
    %----------------------------------------------------------------------
    % Form LHS:
    LHS_CK = LHS_C;
    LHS_CK(2:end-1,:) = LHS_CK(2:end-1,:) + kk^2*FourModeMatrix;
    
    % Form RHS:  
    RHS_FC = ac(:,kk);
    RHS_FS = as(:,kk);
    
    RHS_FC(1) = 0; RHS_FC(end) = 0;
    RHS_FS(1) = 0; RHS_FS(end) = 0;
     
    % Solve:
    B_C(:,kk) = LHS_CK\RHS_FC;
    B_S(:,kk) = LHS_CK\RHS_FS;

end

% For am:
%--------------------------------------------------------------------------
% Form LHS:
LHS_CM = LHS_C;

LHS_CM(2:end-1,:) = LHS_CM(2:end-1,:) + M^2*FourModeMatrix;

% Form RHS:
RHS_FM = aM;
RHS_FM(1) = 0; RHS_FM(end) = 0;

% Solve:
B_M = LHS_CM\RHS_FM;


% For a0:
%--------------------------------------------------------------------------
% Form RHS:
RHS_F0 = a0;
RHS_F0(1) = 0; RHS_F0(end) = 0;

% Solve:
B_0 = LHS_C\RHS_F0;

% Reconstruct Concentration
CicDiff = reconst_fun(N_theta, N_radii, B_0, B_C, B_S, B_M);  
C_ic = [CicDiff CicDiff(:,1)];
end

% #########################################################################
%% SEMI-LAGRANGIAN ADVECTION
% #########################################################################
function C_ic = SL_advection(N_radii,N_theta,tt_adv,deltaT_adv,vel,x,y,thetg,sg,radius,velInterpScheme,conInterpScheme,etol,R1D,C_ic)
% Velocity components at t_n+deltaT_diff and X_a^(n+1)
Vx = reshape(vel(1:end/2,tt_adv*2),(N_radii+1),(N_theta+1));
Vy = reshape(vel(end/2+1:end,tt_adv*2),(N_radii+1),(N_theta+1));

% Do not take the velocity on theta = 2*pi
Vx = Vx(:,1:end-1);
Vy = Vy(:,1:end-1);
x = x(:,1:end-1);
y = y(:,1:end-1);

% Compute the MidPoint coordinates
x_mid = x - deltaT_adv/2 * Vx;
y_mid = y - deltaT_adv/2 * Vy;

% if (max(max(x_mid.^2+y_mid.^2)) > 400) || (min(min(x_mid.^2+y_mid.^2)) < 100)
%     disp('SL: point goes outside the domain')
% end

% Map to the uniform theta-s grid to interpolate velocity at midpoints and t_n + deltaT_diff/2
[s_mid,thet_mid] = arbiTOreg(x_mid,y_mid,radius,R1D);

% Velocity on regular grid at t_n+deltaT_diff/2
VGx = reshape(vel(1:end/2,tt_adv*2-1),(N_radii+1),(N_theta+1));
VGy = reshape(vel(end/2+1:end,tt_adv*2-1),(N_radii+1),(N_theta+1));

% Do not take the velocity on theta = 2*pi
VGx = VGx(:,1:end-1);
VGy = VGy(:,1:end-1);

% Enforce periodicity in theta
VGxEx   = [VGx(:,end-9:end) VGx VGx(:,1:10)];
VGyEx  = [VGy(:,end-9:end) VGy VGy(:,1:10)];

C_ic = C_ic(:,1:end-1);
C_icEx = [C_ic(:,end-9:end) C_ic C_ic(:,1:10)];


% Interpolate Velocity (w/ velInterpScheme)
Vx_mid = interp2(thetg,sg,VGxEx,thet_mid,s_mid,velInterpScheme); 
% Vx_mid(abs(Vx_mid-VGx)<=etol) = VGx(abs(Vx_mid-VGx)<=etol);
Vy_mid = interp2(thetg,sg,VGyEx,thet_mid,s_mid,velInterpScheme);
% Vy_mid(abs(Vy_mid-VGy)<=etol) = VGy(abs(Vy_mid-VGy)<=etol);

% Compute the Departure point coordinates at tn
x_d = x - deltaT_adv * Vx_mid;
y_d = y - deltaT_adv * Vy_mid;

% Check if it leaves the domain

% if (max(max(x_d.^2+y_d.^2)) > 400) || (min(min(x_d.^2+y_d.^2)) < 100)
%     disp('SL: point goes outside the domain')
% end

% Map to theta-s grid and interpolate
[s_d,thet_d] = arbiTOreg(x_d,y_d,radius,R1D);
C_ic = interp2(thetg,sg,C_icEx,thet_d,s_d,conInterpScheme); 
% C_new(abs(C_new-C_ic)<=etol) = C_ic(abs(C_new-C_ic)<=etol); % C at tt_adv*deltaT_adv
% C_ic = [zeros(1,N_theta); C_new; zeros(1,N_theta)];   

C_ic = [C_ic C_ic(:,1)];

end

% #########################################################################
%% CRANK-NICOLSON DIFFUSION
% #########################################################################
function C_ic = CN_diffusion(N_radii,N_theta,M,LHS_C,RHS_C,FourModeMatrix,B_C,B_S,C_ic)

% IC for diffusion problem
CicDiff = C_ic(:,1:end-1);

% Interpolate for spectral coefficients
[bc00, bc_c0, bc_s0, bcM0, ~, ~, ~, ~] = interpFFT(N_theta, N_radii, CicDiff);

for kk = 1 : M-1
    
    % For ac and as:
    %----------------------------------------------------------------------
    % Form LHS:
    LHS_CK = LHS_C;
    LHS_CK(2:end-1,:) = LHS_CK(2:end-1,:) + kk^2*FourModeMatrix;
    
    % Form RHS:  
    RHS_CK = RHS_C; 
    RHS_CK(2:end-1,:) = RHS_CK(2:end-1,:) - kk^2*FourModeMatrix;
    RHS_FC = RHS_CK*bc_c0(:,kk);
    RHS_FS = RHS_CK*bc_s0(:,kk);
    
    RHS_FC(1) = 0; RHS_FC(end) = 0;
    RHS_FS(1) = 0; RHS_FS(end) = 0;
   
    % Solve:
    B_C(:,kk) = LHS_CK\RHS_FC;
    B_S(:,kk) = LHS_CK\RHS_FS;

end

% For am:
%--------------------------------------------------------------------------
% Form LHS:
LHS_CM = LHS_C;

LHS_CM(2:end-1,:) = LHS_CM(2:end-1,:) + M^2*FourModeMatrix;

% Form RHS:
RHS_CM = RHS_C; 
RHS_CM(2:end-1,:) = RHS_CM(2:end-1,:) - M^2*FourModeMatrix;
RHS_FM = RHS_CM*bcM0;
RHS_FM(1) = 0; RHS_FM(end) = 0;

% Solve:
B_M = LHS_CM\RHS_FM;

% For a0:
%--------------------------------------------------------------------------
% Form RHS:
RHS_F0 = RHS_C*bc00;
RHS_F0(1) = 0; RHS_F0(end) = 0;

% Solve:
B_0 = LHS_C\RHS_F0;

% Reconstruct Concentration
CicDiff = reconst_fun(N_theta, N_radii, B_0, B_C, B_S, B_M);  
C_ic = [CicDiff CicDiff(:,1)];
end

% #########################################################################
%% MAPPING FROM ARBITRARY DOMAIN (r,theta) TO UNIFORM DOMAIN (s,theta)
% #########################################################################
function [sp,thetp] = arbiTOreg(xp,yp,radius,R1D)
rinr = radius(1);
rout = radius(2);

% Find theta at the arbitrary point
thetp = atan2(yp,xp);
thetp(:,ceil(end/2)+1:end) = thetp(:,ceil(end/2)+1:end) + 2*pi;
thetp(thetp>=2*pi) = mod(thetp(thetp>=2*pi),2*pi);
thetp(thetp<0)    = mod(thetp(thetp<0),2*pi);

% Find r at the arbitrary point [r0, r1]
rp = (xp.^2+yp.^2).^(0.5);

% Find R at the arbitrary point [-1, 1]
% Rp = 2*(rp-rinr)/(rout-rinr) - 1;
Rp = (R1D(end)-R1D(1))*(rp - rinr)/(rout - rinr) + R1D(1);

% Rp(Rp>1) = 1;
% Rp(Rp<-1) = -1;

% Find s at the arbitrary point
sp = acos(-Rp);
end

% #########################################################################
%% WRITE DATA to BIN FILE
% #########################################################################
function writeData(fileName,dataArr)
fid = fopen(fileName,'a');
fwrite(fid,dataArr,'double');
fclose(fid);
end

% #########################################################################
%% KEEP DIARY
% #########################################################################
function keepDiary(fileName,ntime,tst,tim,blowUp)
fid = fopen(fileName,'a');
fprintf(fid,'%s\n',['Elapsed Time at step = ' num2str(tst) ' of ' num2str(ntime) ' = ' num2str(tim) ' sec']);


if blowUp
    fprintf(fid,'%s\n','!BLOWS UP!');
end    

fclose(fid);
end

% #########################################################################
%% END OF FUNCTION
% #########################################################################
