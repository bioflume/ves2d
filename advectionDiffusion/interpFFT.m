 function [bc0, bc_c, bc_s,bcM,a0,ac,as,aM] = interpFFT(N_theta, N_radii, C_ic)

%-------------------
% Fourier in theta
%-------------------
M = floor((N_theta+1)/2);

CompCoeff = zeros(N_radii+1,N_theta);
a0 = zeros(N_radii+1,1);
ac = zeros(N_radii+1,M-1);
as = zeros(N_radii+1,M-1);
aM = zeros(N_radii+1,1);

for rr = 1 : N_radii+1
    f = C_ic(rr,:);
    CompCoeff(rr,:) = fft(f);
    
    % Fourier Coefficients
    a0(rr,1) = CompCoeff(rr,1)/N_theta;
    ac(rr,:) =  2*real(CompCoeff(rr,2:M))/N_theta;
    as(rr,:) = -2*imag(CompCoeff(rr,2:M))/N_theta;
    aM(rr,1) = real(CompCoeff(rr,M+1))/N_theta;
end

%-------------------
% Chebyshev in r
%-------------------

bc0 = dct(a0);
bcM = dct(aM);

bc_c = zeros(N_radii+1,M-1);
bc_s = zeros(N_radii+1,M-1);

for kk = 1 : M-1
    bc_c(:,kk) = dct(ac(:,kk));
    bc_s(:,kk) = dct(as(:,kk));
end