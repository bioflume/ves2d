clear all;
% load SimpleU
addpath ../src
addpath ../examples

prams.N = 256; theta = (0:prams.N-1)'*2*pi/prams.N;
a = 1; b = 0.5;
X = [a*cos(theta);b*sin(theta)];
prams.nv = 1;
prams.viscCont = 4;
prams.kappa = 1;
options.farField = 'relaxation';

[options,prams] = initVes2D(options,prams);

vesicle = capsules(X,[],[],prams.kappa,prams.viscCont,false);

tt = tstep(options,prams);
op = poten(prams.N);

[Ben,Ten,Div] = vesicle.computeDerivs;
% Ben is the positive fourth arclength derivative
G = op.stokesSLmatrix(vesicle);
D = op.stokesDLmatrix(vesicle); 
%%GOKBERK: has to be stokesDLmatrix, (does not change the result)
% BRYAN: THIS MAKES A HUGE DIFFERENCE.  THE SLP DOES NOT HAVE THE JUMP
% IN ACROSS THE INTERFACE OF THE VESICLE THAT THE DLP DOES
alpha = 0.5*(1+prams.viscCont);
%% GOKBERK: I think this is not correct
 Mat = [alpha*eye(2*prams.N)-D -G*Ten; Div zeros(prams.N)];
 rhs = [-G*Ben*vesicle.X; zeros(prams.N,1)];
% BRYAN: MY LINEAR SYSTEM IS CORRECT.  THERE IS NO BENDING TERM IN THE
% MATRIX BECAUSE WE ARE GIVEN THE POSITION AND AM SOLVING FOR THE
% VELOCITY AND TENSION.
% BRYAN.  I SEEM TO HAVE FORGOT G IN rhs
%% GOKBERK: The following linear system is right, it gives smooth u
%Mat = [alpha*eye(2*prams.N)-D+0*G*Ben -G*Ten; Div zeros(prams.N)];
%rhs = [ -G*Ben*vesicle.X; zeros(prams.N,1)];
%% GOKBERK: However, it does not make any difference whether u is smooth or not??
u_sig = Mat\rhs;
u = u_sig(1:2*prams.N);
sig = u_sig(2*prams.N+1:end);

ntra = 200;

% Velocity along horizontal slice
[xtra,ytra] = meshgrid(linspace(-2*a,2*a,ntra),linspace(0.2,0.2,1));
Xtra = [xtra(:);ytra(:)];
tracers = capsules(Xtra,[],[],1,0,false);
vel = tt.tracersVel(X,sig,u,1,prams.viscCont,...
    [],[],[],Xtra,[]);

velx = vel(1:end/2);
vely = vel(end/2+1:end);

s = find((xtra/a).^2 + (ytra/b).^2 < 1);
velx(s) = velx(s)/prams.viscCont;
vely(s) = vely(s)/prams.viscCont;
% BRYAN: NOT SURE WHY YOU COMMENTED THIS

figure(1); clf; hold on
plot(xtra,velx,'b-o')
plot(xtra(s),velx(s),'r.')

figure(2); clf; hold on
plot(xtra,vely,'b-o')
plot(xtra(s),vely(s),'r.')


% Velocity along vertical slice
[xtra,ytra] = meshgrid(linspace(0.2,0.2,1),linspace(-2*b,2*b,ntra));
Xtra = [xtra(:);ytra(:)];
tracers = capsules(Xtra,[],[],1,0,false);
vel = tt.tracersVel(X,sig,u,1,prams.viscCont,...
    [],[],[],Xtra,[]);

velx = vel(1:end/2);
vely = vel(end/2+1:end);

s = find((xtra/a).^2 + (ytra/b).^2 < 1);
velx(s) = velx(s)/prams.viscCont;
vely(s) = vely(s)/prams.viscCont;
% BRYAN: NOT SURE WHY YOU COMMENTED THIS

figure(3); clf; hold on
plot(ytra,velx,'b-o')
plot(ytra(s),velx(s),'r.')
%plot(-b,u(3*prams.N/4+1),'ro')
%plot(b,u(prams.N/4+1),'ro')

figure(4); clf; hold on
plot(ytra,vely,'b-o')
plot(ytra(s),vely(s),'r.')
%plot(-b,u(3*prams.N/4+1+prams.N),'ro')
%plot(b,u(prams.N/4+1+prams.N),'ro')

% Surface plot of the velocity field
[xtra,ytra] = meshgrid(linspace(-2*a,2*a,ntra),linspace(-2*b,2*b,ntra));
Xtra = [xtra(:);ytra(:)];
tracers = capsules(Xtra,[],[],1,0,false);
vel = tt.tracersVel(X,sig,u,1,prams.viscCont,...
    [],[],[],Xtra,[]);

velx = vel(1:end/2);
vely = vel(end/2+1:end);

s = find((xtra/a).^2 + (ytra/b).^2 < 1);
velx(s) = velx(s)/prams.viscCont;
vely(s) = vely(s)/prams.viscCont;

velx = reshape(velx,ntra,ntra);
vely = reshape(vely,ntra,ntra);

v = (velx.^2 + vely.^2).^0.5;

figure(5); clf; hold on;
surf(xtra,ytra,v,'edgecolor','none'); shading interp; colorbar;


