% Looks at the effect of aliasing on computing the traction jump
addpath ../src
addpath ../examples

phi = pi/5;

Nup = 1024;
theta = (0:Nup-1)'*2*pi/Nup;
%r = 1 + .2*cos(5*theta);
%x = r.*cos(theta); y = r.*sin(theta);
x = 3*cos(theta)*cos(phi) + sin(theta)*sin(phi); 
y = -3*cos(theta)*sin(phi) + sin(theta)*cos(phi);
vesicleUp = capsules([x;y],[],[],1,1,0);
%sigmaUp = cos(theta);
sigmaUp = zeros(Nup,1);

N = 32;
theta = (0:N-1)'*2*pi/N;
%r = 1 + .2*cos(5*theta);
%x = r.*cos(theta); y = r.*sin(theta);
x = 3*cos(theta)*cos(phi) + sin(theta)*sin(phi); 
y = -3*cos(theta)*sin(phi) + sin(theta)*cos(phi);
vesicle = capsules([x;y],[],[],1,1,0);
%sigma = cos(theta);
sigma = zeros(N,1);

modesUp = (-Nup/2:Nup/2-1)';
fup = vesicleUp.tracJump(vesicleUp.X,sigmaUp);
fzup = fup(1:end/2) + 1i*fup(end/2+1:end);
fzuph = fftshift(fft(fzup))/Nup;
% build Fourier modes of bending of upscaled vesicle shape

modes = (-N/2:N/2-1)';
f1 = vesicle.tracJump(vesicle.X,sigma);
fz1 = f1(1:end/2) + 1i*f1(end/2+1:end);
fz1h = fftshift(fft(fz1))/N;
% build Fourier modes of bending of low reoslution vesicle without
% anti-aliasing

vesicle.antiAlias = true;
op = poten(N,4,6*N);
vesicle.setUpRate(op);
% choose the upsampling rate
f2 = vesicle.tracJump(vesicle.X,sigma);
fz2 = f2(1:end/2) + 1i*f2(end/2+1:end);
fz2h = fftshift(fft(fz2))/N;
% build Fourier modes of bending of low reoslution vesicle with
% anti-aliasing


as = find(modesUp == modes(1));
ae = find(modesUp == modes(end));

figure(1); clf
semilogy(modes,abs(fzuph(as:ae)),'b-')
% plot the "true" Fourier modes
hold on
semilogy(modes,abs(fz1h),'ro')
% plot the aliased Fourier modes
semilogy(modes,abs(fz2h),'go')
% plot the anti-aliased Fourier modes
xlim([-N/2+1 N/2-1])
ylim([1e-4 1e3])
legend('exact','aliased','anti-aliased')

figure(2); clf
semilogy(modes,abs(fzuph(as:ae) - fz1h)/max(abs(fzuph)),'ro')
% plot the error of the aliased Fourier modes
hold on
semilogy(modes,abs(fzuph(as:ae) - fz2h)/max(abs(fzuph)),'go')
% plot the error of the anti-aliased Fourier modes
xlim([-N/2+1 N/2-1])
legend('aliased','anti-aliased')

err1 = abs(fzuph(as:ae) - fz1h)/max(abs(fzuph));
err2 = abs(fzuph(as:ae) - fz2h)/max(abs(fzuph));
disp(vesicle.uprate)










