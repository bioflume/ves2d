clear all;
clc;

Ntheta = 1024;
Nradii = 255;

radii = [10; 20];

[x,y,r1D,theta,~,~,~] = generateGrid(Ntheta,Nradii,radii);

[THETA,R] = meshgrid([theta;2*pi],r1D);

THETA = THETA(:,1:end-1);

% gradUth = 1/300 * (-100-4e4./R.^2);
% Uth = 1/300 * (-100*R + 100*400./R);
load VF40Vel0.mat
%%% GET Uth Ur from UvesX UvesY

% 1) Blob IC
Cblob = zeros(Nradii+1,Ntheta+1);
Cblob(((x-15).^2 + y.^2)<16) = 1;

%2) Ring IC
Cring = zeros(Nradii+1,Ntheta+1);
Cring(R>=14&R<=16) = 1;

%3) Vesicle IC
load IC_t40_Largest
% load RandIC_STDP17
Cves = C_ic;

%4) Moon IC
Cmoon = zeros(Nradii+1,Ntheta+1);
Cmoon(x<=0) = 1;

%5) angular patch
Cang = zeros(Nradii+1,Ntheta+1);
Cang(THETA>=0 & THETA<=pi/18) = 1;

%6) Random IC
load RandIC_STDP17
Crand = C_ic;

%7) Serpentine ring IC
load serpentineRing
Cserp = C_ic;


% To compute derivative of C w.r.t. r

for k = 1 : Ntheta+1
    
    pp = spline(r1D,Cblob(:,k));
    p_der = fnder(pp,1);
    dCblobdr(:,k) = ppval(p_der,r1D); 
    
    pp = spline(r1D,Cring(:,k));
    p_der = fnder(pp,1);
    dCringdr(:,k) = ppval(p_der,r1D);
    
    pp = spline(r1D,Cves(:,k));
    p_der = fnder(pp,1);
    dCvesdr(:,k) = ppval(p_der,r1D);
    
    pp = spline(r1D,Cmoon(:,k));
    p_der = fnder(pp,1);
    dCmoondr(:,k) = ppval(p_der,r1D);
    
    pp = spline(r1D,Cang(:,k));
    p_der = fnder(pp,1);
    dCangdr(:,k) = ppval(p_der,r1D);
    
    pp = spline(r1D,Crand(:,k));
    p_der = fnder(pp,1);
    dCranddr(:,k) = ppval(p_der,r1D);
    
    pp = spline(r1D,Cserp(:,k));
    p_der = fnder(pp,1);
    dCserpdr(:,k) = ppval(p_der,r1D);
    
    pp = spline(r1D,Uth(:,k));
    p_der = fnder(pp,1);
    dUthdr(:,k) = ppval(p_der,r1D);
    
    
end

for k = 1 : Nradii + 1
    pp = spline(theta,Cblob(k,1:end-1));
    p_der = fnder(pp,1);
    dCblobdt(k,:) = ppval(p_der,theta);
    
    pp = spline(theta,Cring(k,1:end-1));
    p_der = fnder(pp,1);
    dCringdt(k,:) = ppval(p_der,theta);
    
    pp = spline(theta,Cves(k,1:end-1));
    p_der = fnder(pp,1);
    dCvesdt(k,:) = ppval(p_der,theta);
    
    pp = spline(theta,Cmoon(k,1:end-1));
    p_der = fnder(pp,1);
    dCmoondt(k,:) = ppval(p_der,theta);
    
    pp = spline(theta,Cang(k,1:end-1));
    p_der = fnder(pp,1);
    dCangdt(k,:) = ppval(p_der,theta);
    
    pp = spline(theta,Crand(k,1:end-1));
    p_der = fnder(pp,1);
    dCranddt(k,:) = ppval(p_der,theta);
    
    pp = spline(theta,Cserp(k,1:end-1));
    p_der = fnder(pp,1);
    dCserpdt(k,:) = ppval(p_der,theta);
    
    pp = spline(theta,Uth(k,1:end-1));
    p_der = fnder(pp,1);
    dUthdt(k,:) = ppval(p_der,theta);
    
end


% Compute L2 norms
gradCblob = dCblobdr(:,1:end-1) + 1./R(:,1:end-1).*dCblobdt;
gradCring = dCringdr(:,1:end-1) + 1./R(:,1:end-1).*dCringdt;
gradCves  = dCvesdr(:,1:end-1)  + 1./R(:,1:end-1).*dCvesdt ;
gradCmoon  = dCmoondr(:,1:end-1)  + 1./R(:,1:end-1).*dCmoondt ;
gradCang  = dCangdr(:,1:end-1)  + 1./R(:,1:end-1).*dCangdt ;
gradCrand  = dCranddr(:,1:end-1)  + 1./R(:,1:end-1).*dCranddt ;
gradCserp  = dCserpdr(:,1:end-1)  + 1./R(:,1:end-1).*dCserpdt ;
gradUth  = dUthdr(:,1:end-1)  + 1./R(:,1:end-1).*dUthdt ;


L2dCblob = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCblob.^2.*R(:,1:end-1))));
L2dCring = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCring.^2.*R(:,1:end-1))));
L2dCves  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCves.^2 .*R(:,1:end-1))));
L2dCmoon  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCmoon.^2 .*R(:,1:end-1))));
L2dCang  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCang.^2 .*R(:,1:end-1))));
L2dCrand  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCrand.^2 .*R(:,1:end-1))));
L2dCserp  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCserp.^2 .*R(:,1:end-1))));

L2Uth = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,Uth(:,1:end-1).^2 .*R(:,1:end-1))));

% Compute the measures


M4Cblob = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCblobdt.*Uth(:,1:end-1))/(L2dCblob*L2dU).*R(:,1:end-1)));

M4Cring = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCringdt.*Uth(:,1:end-1))/(L2dCring*L2dU).*R(:,1:end-1)));

M4Cves = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCvesdt.*Uth(:,1:end-1))/(L2dCves*L2dU).*R(:,1:end-1)));

M4Cmoon = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCmoondt.*Uth(:,1:end-1))/(L2dCmoon*L2dU).*R(:,1:end-1)));

M4Cang = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCangdt.*Uth(:,1:end-1))/(L2dCang*L2dU).*R(:,1:end-1)));

M4Crand = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCranddt.*Uth(:,1:end-1))/(L2dCrand*L2dU).*R(:,1:end-1)));

M4Cserp = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCserpdt.*Uth(:,1:end-1))/(L2dCserp*L2dU).*R(:,1:end-1)));


disp('-------------------')
disp('Measure #4 for ICs')
disp('-------------------')

disp(['Ring IC: ' num2str(M4Cring)])
disp(['Blob IC: ' num2str(M4Cblob)])
disp(['Vesicle IC: ' num2str(M4Cves)])
disp(['Moon IC: ' num2str(M4Cmoon)])
disp(['Angular patch IC: ' num2str(M4Cang)])
disp(['Random IC: ' num2str(M4Crand)])
disp(['Serpentine IC: ' num2str(M4Cserp)])

% disp('-------------------')
% disp('(H1norm(C)-norm(C))/norm(C) for ICs')
% disp('-------------------')
% 
% disp(['Ring IC: ' num2str((H1Cring-L1Cring)/L1Cring)])
% disp(['Blob IC: ' num2str((H1Cblob-L1Cblob)/L1Cblob)])
% disp(['Vesicle IC: ' num2str((H1Cves-L1Cves)/L1Cves)])


