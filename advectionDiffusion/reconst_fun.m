function C = reconst_fun(N_theta, N_radii, bc0, bc_c, bc_s,bcM)

M = floor((N_theta+1)/2);
%-------------------
% Chebyshev in r
%-------------------

a0 = idct(bc0);
aM = idct(bcM);

ac = zeros(N_radii+1,M-1);
as = zeros(N_radii+1,M-1);

for kk = 1 : M-1
    ac(:,kk) = idct(bc_c(:,kk));
    as(:,kk) = idct(bc_s(:,kk));
end
    
%-------------------
% Fourier in theta
%-------------------
C = zeros(N_radii+1,N_theta);
CompCoeff = zeros(N_radii+1,N_theta);
for rr = 1 : N_radii+1
    CompCoeff(rr,1) = a0(rr)*N_theta;
    
    CompCoeff(rr,2:M) = (ac(rr,:) - 1i*as(rr,:))*N_theta/2;
    CompCoeff(rr,M+1) = aM(rr)*N_theta;
    CompCoeff(rr,M+2:end) = conj(fliplr(CompCoeff(rr,2:M)));
    C(rr,:) = ifft(CompCoeff(rr,:));
end