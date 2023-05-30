clear; 

% Original RSlets
X = rand(128,2);

x = X(1:end/2,:);
y = X(end/2+1:end,:);

cx = 2; cy = 1;
stokeslet = rand(2,1);
rotlet = rand(1,1);

% the center of the rotlet/stokeslet terms

rho2 = (x-cx).^2 + (y-cy).^2;
% distance squared

LogTerm = -0.5*log(rho2)*stokeslet(1);
rorTerm = 1./rho2.*((x-cx).*(x-cx)*stokeslet(1) + ...
    (x-cx).*(y-cy)*stokeslet(2));
RotTerm = (y-cy)./rho2*rotlet;
velx = 1/4/pi*(LogTerm + rorTerm) + RotTerm;
% x component of velocity due to the stokeslet and rotlet

LogTerm = -0.5*log(rho2)*stokeslet(2);
rorTerm = 1./rho2.*((y-cy).*(x-cx)*stokeslet(1) + ...
    (y-cy).*(y-cy)*stokeslet(2));
RotTerm = -(x-cx)./rho2*rotlet;
vely = 1/4/pi*(LogTerm + rorTerm) + RotTerm;
% y component of velocity due to the stokeslet and rotlet

vel = [velx;vely];


%% New Code

xfactor = 1/4/pi;

ntarget = 2;
nsource = 1;
Ntarget = 64;

pot = zeros(2*Ntarget,ntarget);

for k = 1 : ntarget
  
  K = [1:nsource];
  
  for ksrc = K
    diffx = x(:,k)-cx(ksrc);
    diffy = y(:,k)-cy(ksrc);
    
    rho2 = diffx.^2 + diffy.^2; 

    LogTerm = -0.5*log(rho2)*stokeslet(1,ksrc);
    rorTerm = 1./rho2.*(diffx.*diffx*stokeslet(1,ksrc) + ...
        diffx.*diffy*stokeslet(2,ksrc));
    RotTerm = diffy./rho2*rotlet(ksrc);
    

    pot(1:Ntarget,k) = pot(1:Ntarget,k) + xfactor*(LogTerm + rorTerm) + RotTerm;

    LogTerm = -0.5*log(rho2)*stokeslet(2,ksrc);
    rorTerm = 1./rho2.*(diffy.*diffx*stokeslet(1,ksrc) + ...
        diffy.*diffy*stokeslet(2,ksrc));
    RotTerm = -diffx./rho2*rotlet(ksrc);

    pot(Ntarget+1:2*Ntarget,k) = pot(Ntarget+1:2*Ntarget,k) + xfactor*(LogTerm + rorTerm) + RotTerm;

  end
end

%%

GR = zeros(2*Ntarget,3,ntarget);
isrc = 1;
for itr = 1 : ntarget

diffx = x(:,itr) - cx(isrc); 
diffy = y(:,itr) - cy(isrc);

rho = sqrt(diffx.^2 + diffy.^2);
rho2 = (diffx.^2 + diffy.^2).^(-1);

logpart = -rho;

GR(1:Ntarget,1,itr) = logpart + diffx.^2.*rho2;
GR(Ntarget+1:2*Ntarget,1,itr) =  diffx.*diffy.*rho2;

GR(1:Ntarget,2,itr) = diffx.*diffy.*rho2;
GR(Ntarget+1:2*Ntarget,2,itr) = logpart + diffy.^2.*rho2;

GR(1:Ntarget,3,itr) = diffy.*rho2;
GR(Ntarget+1:2*Ntarget,3,itr) =  -diffx.*rho2;
  
end

GR = xfactor * GR;


