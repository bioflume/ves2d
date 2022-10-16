clear op;
addpath ../src

geom = 'geom1';
if strcmp(geom,'geom1')
  a = 3; b = 1;
  geom = @(omega) [a*cos(omega);b*sin(omega)];
elseif strcmp(geom,'geom2')
  a = 1; b = 3; c = 0.85;
  r = @(omega) 0.5*sqrt( (a*cos(omega)).^2 + (b*sin(omega)).^2) + ...
        .07*cos(12*(omega));
  geom = @(omega) [c*r(omega).*cos(omega);r(omega).*sin(omega)];
elseif strcmp(geom,'geom3')
  r = @(omega) 1 + 0.2*cos(5*omega);
  geom = @(omega) [r(omega).*cos(omega);r(omega).*sin(omega)];
end

bc = 'bc1';
if strcmp(bc,'bc1')
  f = @(omega) [exp(cos(omega));exp(cos(sin(omega)))]; 
elseif strcmp(bc,'bc2')
  f = @(omega) [sin(omega);cos(omega)];
elseif strcmp(bc,'bc3')
  f = @(omega) [cos(2*omega);sin(2*omega)];
end

lp = 'slp';

Nup = 512; thetaUp = (0:Nup-1)'*2*pi/Nup; 
fup = f(thetaUp);
modesUp = (-Nup/2+1:Nup/2-1)';
Xup = geom(thetaUp);
if strcmp(bc,'bc1')
  vesicleUp = capsules(Xup,[],[],1,0,false,6);
  op = poten(Nup,false);
  if strcmp(lp,'slp')
    Gup = op.stokesSLmatrix(vesicleUp);
  elseif strcmp(lp,'dlp')
    Gup = op.stokesDLmatrix(vesicleUp);
  end
  zup = Gup*fup; zup = zup(1:end/2);
elseif strcmp(bc,'bc2')
  zup = Gup*fup; zup = zup(1:end/2);
elseif strcmp(bc,'bc3')
  zup = Gup*fup; zup = zup(1:end/2);
end

N = 8; theta = (0:N-1)'*2*pi/N;
f1 = f(theta);
f2 = f(theta);
modes = (-N/2+1:N/2-1)';
X = geom(theta);
vesicle1 = capsules(X,[],[],1,0,false,6);
vesicle2 = capsules(X,[],[],1,0,false,6);
vesicle2.antiAlias = true; 
% for now, don't want to use anti-aliasing to form Jacobian, tangent
% vectors, etc.  Just want to use it when computing the SLP
op = poten(N,4,3.0*N);
if strcmp(lp,'slp')
  G1 = op.stokesSLmatrix(vesicle1);
elseif strcmp(lp,'dlp')
  G1 = op.stokesDLmatrix(vesicle1);
end
z1 = G1*f1; z1 = z1(1:end/2);

%vesicle2.setUpRate(op);
%vesicle2.uprate = 6;

tstart = tic;
if strcmp(lp,'slp');
  z2 = op.exactStokesSLdiag(vesicle2,[],f2); 
elseif strcmp(lp,'dlp');
  z2 = op.exactStokesDLdiag(vesicle2,[],f2); 
end
z2 = z2(1:end/2);


zuph = fftshift(fft(zup))/numel(zup); zuph = zuph(2:end);
z1h = fftshift(fft(z1))/numel(z1); z1h = z1h(2:end);
z2h = fftshift(fft(z2))/numel(z2); z2h = z2h(2:end);
% compute fft of each layer potential and remove the first aliased
% frequency

s = find(modesUp >= modes(1) & modesUp <= modes(end));

figure(1); clf;
semilogy(modesUp(s),abs(zuph(s)),'k')
hold on
%semilogy(modes,abs(z1h),'bo')
semilogy(modes,abs(z2h),'ro')
xlim([modes(1)-1 modes(end)+1])

%figure(2); clf;
figure(2);
semilogy(modes,abs(z1h - zuph(s)),'b-o')
hold on
semilogy(modes,abs(z2h - zuph(s)),'r-o')

%totalErr = norm(z2h - zuph(s));
slow = find(abs(modes) < N/4);
shigh = find(abs(modes) >= N/4);

lowErr = norm(z2h(slow) - zuph(s(slow)));
highErr = norm(z2h(shigh) - zuph(s(shigh)));

%fprintf('Total Error in energy is %4.2e\n',totalErr);
fprintf('Error in energy in low frequencies is  %4.2e\n',...
    lowErr);
fprintf('Error in energy in high frequencies is %4.2e\n',...
    highErr);




