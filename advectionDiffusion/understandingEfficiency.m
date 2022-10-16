clear all;
clc;

Ntheta = 1024;
Nradii = 255;

radii = [10; 20];

[x,y,r1D,theta,~,~,~] = generateGrid(Ntheta,Nradii,radii);

[THETA,R] = meshgrid([theta;2*pi],r1D);

THETA = THETA(:,1:end-1);

gradUth = 1/300 * (-100-4e4./R.^2);
Uth = 1/300 * (-100*R + 100*400./R);

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
    
end


% Compute L2 norms
gradCblob = dCblobdr(:,1:end-1) + 1./R(:,1:end-1).*dCblobdt;
gradCring = dCringdr(:,1:end-1) + 1./R(:,1:end-1).*dCringdt;
gradCves  = dCvesdr(:,1:end-1)  + 1./R(:,1:end-1).*dCvesdt ;
gradCmoon  = dCmoondr(:,1:end-1)  + 1./R(:,1:end-1).*dCmoondt ;
gradCang  = dCangdr(:,1:end-1)  + 1./R(:,1:end-1).*dCangdt ;
gradCrand  = dCranddr(:,1:end-1)  + 1./R(:,1:end-1).*dCranddt ;
gradCserp  = dCserpdr(:,1:end-1)  + 1./R(:,1:end-1).*dCserpdt ;

L2dCblob = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCblob.^2.*R(:,1:end-1))));
L2dCring = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCring.^2.*R(:,1:end-1))));
L2dCves  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCves.^2 .*R(:,1:end-1))));
L2dCmoon  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCmoon.^2 .*R(:,1:end-1))));
L2dCang  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCang.^2 .*R(:,1:end-1))));
L2dCrand  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCrand.^2 .*R(:,1:end-1))));
L2dCserp  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradCserp.^2 .*R(:,1:end-1))));

L2Cblob = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,Cblob(:,1:end-1).^2.*R(:,1:end-1))));
L2Cring = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,Cring(:,1:end-1).^2.*R(:,1:end-1))));
L2Cves  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,Cves(:,1:end-1).^2 .*R(:,1:end-1))));
L2Cmoon  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,Cmoon(:,1:end-1).^2 .*R(:,1:end-1))));
L2Cang  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,Cang(:,1:end-1).^2 .*R(:,1:end-1))));
L2Crand  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,Crand(:,1:end-1).^2 .*R(:,1:end-1))));
L2Cserp  = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,Cserp(:,1:end-1).^2 .*R(:,1:end-1))));

L1dCblob = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(gradCblob).*R(:,1:end-1))));
L1dCring = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(gradCring).*R(:,1:end-1))));
L1dCves  = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(gradCves) .*R(:,1:end-1))));
L1dCmoon  = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(gradCmoon) .*R(:,1:end-1))));
L1dCang  = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(gradCang) .*R(:,1:end-1))));
L1dCrand  = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(gradCrand) .*R(:,1:end-1))));
L1dCserp  = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(gradCserp) .*R(:,1:end-1))));

L1Cblob = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(Cblob(:,1:end-1)).*R(:,1:end-1))));
L1Cring = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(Cring(:,1:end-1)).*R(:,1:end-1))));
L1Cves  = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(Cves(:,1:end-1)) .*R(:,1:end-1))));
L1Cmoon  = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(Cmoon(:,1:end-1)) .*R(:,1:end-1))));
L1Cang  = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(Cang(:,1:end-1)) .*R(:,1:end-1))));
L1Crand  = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(Crand(:,1:end-1)) .*R(:,1:end-1))));
L1Cserp  = 1/(pi*300) * (trapz(theta,trapz(r1D,abs(Cserp(:,1:end-1)) .*R(:,1:end-1))));

H1Cblob = 1/sqrt(pi*300)*MeasureMixingLaplacian(Cblob(:,1:end-1),Nradii,Ntheta,radii,[]);
H1Cring = 1/sqrt(pi*300)*MeasureMixingLaplacian(Cring(:,1:end-1),Nradii,Ntheta,radii,[]);
H1Cves = 1/sqrt(pi*300)*MeasureMixingLaplacian(Cves(:,1:end-1),Nradii,Ntheta,radii,[]);
H1Cmoon = 1/sqrt(pi*300)*MeasureMixingLaplacian(Cmoon(:,1:end-1),Nradii,Ntheta,radii,[]);
H1Cang = 1/sqrt(pi*300)*MeasureMixingLaplacian(Cang(:,1:end-1),Nradii,Ntheta,radii,[]);
H1Crand = 1/sqrt(pi*300)*MeasureMixingLaplacian(Crand(:,1:end-1),Nradii,Ntheta,radii,[]);
H1Cserp = 1/sqrt(pi*300)*MeasureMixingLaplacian(Cserp(:,1:end-1),Nradii,Ntheta,radii,[]);

L2dUthet = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,gradUth(:,1:end-1).^2.*R(:,1:end-1))));
L2dU = 1/sqrt(pi*300) * sqrt(trapz(theta,trapz(r1D,Uth(:,1:end-1).^2.*R(:,1:end-1))));

% Compute the measures
% M2Cblob = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCblobdr(:,1:end-1).*gradUth(:,1:end-1))/(L2dCblob*L2dUthet).*R(:,1:end-1)))/L2dCblob;
% 
% M2Cring = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCringdr(:,1:end-1).*gradUth(:,1:end-1))/(L2dCring*L2dUthet).*R(:,1:end-1)))/L2dCring;
% 
% M2Cves = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCvesdr(:,1:end-1).*gradUth(:,1:end-1))/(L2dCves*L2dUthet).*R(:,1:end-1)))/L2dCves;
% 
% M1Cblob = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCblobdr(:,1:end-1).*gradUth(:,1:end-1))/(L2dCblob*L2dUthet).*R(:,1:end-1)));
% 
% M1Cring = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCringdr(:,1:end-1).*gradUth(:,1:end-1))/(L2dCring*L2dUthet).*R(:,1:end-1)));
% 
% M1Cves = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCvesdr(:,1:end-1).*gradUth(:,1:end-1))/(L2dCves*L2dUthet).*R(:,1:end-1)));
% 
% M3Cblob = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCblobdr(:,1:end-1).*gradUth(:,1:end-1))/(L2dCblob*L2dUthet).*R(:,1:end-1)))*(H1Cblob-L1Cblob)/L1Cblob;
% 
% M3Cring = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCringdr(:,1:end-1).*gradUth(:,1:end-1))/(L2dCring*L2dUthet).*R(:,1:end-1)))*(H1Cring-L1Cring)/L1Cring;
% 
% M3Cves = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCvesdr(:,1:end-1).*gradUth(:,1:end-1))/(L2dCves*L2dUthet).*R(:,1:end-1)))*(H1Cves-L1Cves)/L1Cves;

M4Cblob = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCblobdt.*Uth(:,1:end-1))/(L2dCblob*L2dU).*R(:,1:end-1)));

M4Cring = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCringdt.*Uth(:,1:end-1))/(L2dCring*L2dU).*R(:,1:end-1)));

M4Cves = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCvesdt.*Uth(:,1:end-1))/(L2dCves*L2dU).*R(:,1:end-1)));

M4Cmoon = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCmoondt.*Uth(:,1:end-1))/(L2dCmoon*L2dU).*R(:,1:end-1)));

M4Cang = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCangdt.*Uth(:,1:end-1))/(L2dCang*L2dU).*R(:,1:end-1)));

M4Crand = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCranddt.*Uth(:,1:end-1))/(L2dCrand*L2dU).*R(:,1:end-1)));

M4Cserp = 1/(pi*300) * trapz(theta,trapz(r1D,abs(dCserpdt.*Uth(:,1:end-1))/(L2dCserp*L2dU).*R(:,1:end-1)));

% Report them
% disp('-------------------')
% disp('Measure #1 for ICs')
% disp('-------------------')
% 
% disp(['Ring IC: ' num2str(M1Cring)])
% disp(['Blob IC: ' num2str(M1Cblob)])
% disp(['Vesicle IC: ' num2str(M1Cves)])
% 
% disp('-------------------')
% disp('Measure #2 for ICs')
% disp('-------------------')
% 
% disp(['Ring IC: ' num2str(M2Cring)])
% disp(['Blob IC: ' num2str(M2Cblob)])
% disp(['Vesicle IC: ' num2str(M2Cves)])
% 
% disp('-------------------')
% disp('Measure #3 for ICs')
% disp('-------------------')
% 
% disp(['Ring IC: ' num2str(M3Cring)])
% disp(['Blob IC: ' num2str(M3Cblob)])
% disp(['Vesicle IC: ' num2str(M3Cves)])

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


