set(0,'DefaultAxesFontSize',22)
options.tracers = false; % plot saved tracers
options.pressure = false; % plot saved pressures
options.stress = false; % plot the saved stresses
options.center = false; % move vesicle to x-axis for easier viewing
options.quiver = false;
options.jacobian = false;
options.dist = false;
options.angle = false;

irate = 1; % controls the speed of the visualization

if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/RestartData.bin .
  file = 'RestartData.bin';
  confined = false;
  ax = 4*[-1 1 -1 1];
  options.dist = true;
end

if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/diffuserData.bin .
  file = 'diffuserData.bin';
  confined = true;
  ax = [-7 7 -2 2];
  ax = [-12 12 -2 2];
end

if 0
  !scp quaife@ronaldo.ices.utexas.edu:Ves2D/examples/output/couetteData.bin .
  file = 'couetteData.bin';
  confined = true;
  ax = 21*[-1 1 -1 1];
end

if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/couetteGokberkData.bin .
  file = 'couetteGokberkData.bin';
  confined = true;
  ax = 21*[-1 1 -1 1];
end

if 0
  !scp quaife@ronaldo.ices.utexas.edu:/scratch/quaife/Ves2D/examples/output/extensional2VesData.bin .
  file = 'extensional2VesData.bin';
  confined = false;
  ax = 1/2*[-6 6 -6 6];
  options.dist = true;
end

if 0
%   !scp -q ronaldo.ices.utexas.edu:Ves2D/examples/output/choke1VesData.bin .
  file = 'choke1VesTestN16Data.bin';
  ax = [-11 11 -5 5];
  options.confined = true;
  % Single vesicle shear flow
end

if 1
  file = 'shear1VesData.bin';
  ax = 3*[-1 1 -1 1];
  options.confined = false;
end

if 0
%  !scp quaife@ronaldo.ices.utexas.edu:Ves2D/examples/output/shear1VesData.bin .
  file = 'shear1VesData.bin';
  options.pressure = false;
  options.stress = false;
  nx = 50; ny = 50;
  confined = false;
  ax = [-4 4 -4 4];
  options.angle = true;
  options.jacobian = true;
end


if 0
%  !scp quaife@ronaldo.ices.utexas.edu:Ves2D/examples/output/relaxation1VesData.bin .
  file = 'relaxation1VesData.bin';
  confined = false;
  ax = 2*[-2 2 -2 2];
end


if 0
%  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/shear2VesData.bin .
  file = 'shear2VesData.bin';
  confined = false;
  ax = [-12 12 -4 4];
end

if 0
%  !scp quaife@ronaldo.ices.utexas.edu:/scratch/quaife/Ves2D/examples/output/shear4VesData.bin .
  file = 'shear4VesData.bin';
  confined = false;
  ax = [-8 8 -5 5];
  options.dist = true;
end

if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/shear16VesData.bin .
  file = 'shear16VesData.bin';
  confined = false;
  ax = 1/2*[-10 10 -18 18];
end

if 0
%  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/shear20VesData.bin .
  file = 'shear20VesData.bin';
  confined = false;
%  ax = [-15 15 -4 4];
  ax = [-0 6 -3 3];
end

if 0
  !scp quaife@ronaldo.ices.utexas.edu:/scratch/quaife/Ves2D/examples/output/shear100VesData.bin .
  file = 'shear100VesData.bin';
  confined = false;
  ax = [-2 25 -8 8];
%  ax = [-150 150 -12 12];
end

if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/cShapeData.bin .
  file = 'cShapeData.bin';
  confined = false;
  ax = [-20 6 -3 3];
end

if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/microfluidicData.bin .
  file = 'microfluidicData.bin';
  confined = true;
  ax = [-2 2 -2 2];
end

if 0
%  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/taylorGreen9VesDataRefined.bin .
  file = 'taylorGreen9VesData.bin';
  confined = false;
  ax = [-0.1 pi+0.1 -0.1 pi+0.1];
end


if 0
  !scp quaife@ronaldo.ices.utexas.edu:Ves2D/examples/output/choke1VesLibinData.bin .
  file = 'choke1VesLibinData.bin';
  confined = true;
%  ax = [-40 40 -4 4];
  ax = [-10 10 -4 4];
end

if 0
%  !scp quaife@ronaldo.ices.utexas.edu:Ves2D/examples/output/choke1Ves2Data.bin .
  file = 'choke1Ves2Data.bin';
  confined = true;
  ax = [-40 40 -4 4];
%  ax = [-10 10 -4 4];
end

if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/parabolic15VesData.bin .
  file = 'parabolic15VesData.bin';
  confined = false;
  ax = [20 40 -10 10];
end

if 0
%  !scp quaife@ronaldo.ices.utexas.edu:/scratch/quaife/Ves2D/examples/output/parabolic1VesData.bin .
%   file = 'parabolic1VesParacN64Data.bin';
  file = 'parabolic1VesData.bin';
  confined = false;
  ax = 20*[-1 1 -1 1];
  irate = 1;
  options.center = true;
end

if 0
  file = 'choke1VesAdap2Data.bin';
  command = ['!scp ronaldo.ices.utexas.edu:Ves2D/examples/output/' ...
        file ' .'];
  eval(command);
  confined = true;
  ax = [-10 10 -4 4];
end


if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/figureEightData.bin .
  file = 'figureEightData.bin';
  confined = true;
  ax = [-2.1 2.1 -1.1 1.1];
end


if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/ex6_couetteData.bin .
  file = 'ex6_couetteData.bin';
  confined = true;
  ax = [-24 24 -24 24];
end


if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/tubeData.bin .
  file = 'tubeData.bin';
  confined = true;
  ax = [-10.5 10.5 -3.5 3.5];
end


if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/cylinderData.bin .
  file = 'cylinderData.bin';
  confined = true;
  ax = 11*[-1 1 -1 1];
end


if 0
  !scp ronaldo.ices.utexas.edu:Ves2D/examples/output/couette4plyData.bin .
  file = 'couette4plyData.bin';
  confined = true;
  ax = [-21 21 -21 21];
end


[posx,posy,ten,wallx,wally,ea,el,time,n,nv] = loadFile(file);
% load positions, tension, errors, time, number of points, and number
% of vesicles

ntime = numel(time);

if options.tracers
%  fid = fopen('couetteTracers.bin','r');
  fid = fopen('couetteWithoutVesiclesData.bin','r');
  tracers = fread(fid,'double');
  ntime = numel(tracers)/2/7550;
%  ntracers = numel(tracers)/2/ntime;
%  ntracers = 7269;
%  ntracers = 450;
  ntracers = 7550;
  tracerx = zeros(ntracers,ntime);
  tracery = zeros(ntracers,ntime);
  for k = 1:ntime
    istart = (k-1)*2*ntracers + 1;
    iend = istart + ntracers - 1;
    tracerx(:,k) = tracers(istart:iend);
    istart = istart + ntracers;
    iend = iend + ntracers;
    tracery(:,k) = tracers(istart:iend);
  end
end

if options.pressure
  fileName = ['output/' file(1:end-8) 'Pressure.bin'];
  fid = fopen(fileName,'r');
  pressure = fread(fid,'double');
%  nPress = numel(pressure)/(ntime+2-1);
  nPress = nx*ny;
  pressX = pressure(1:nPress);
  pressY = pressure(nPress+1:2*nPress);
%  press = zeros(nPress,ntime-1);
  press = [];
  for k = 1:ntime - 1
    istart = (k+1)*nPress + 1;
    iend = istart + nPress - 1;
    if iend <= numel(pressure)
      press(:,k) = pressure(istart:iend);
    end
  end
  xmin = min(min(pressX));
  xmax = max(max(pressX));
  ymin = min(min(pressY));
  ymax = max(max(pressY));
%  cmin = min(min(press(:,end-1:end)));
%  cmax = max(max(press(:,end-1:end)));


  if 1
    pressX = reshape(pressX,nx,ny);
    pressY = reshape(pressY,nx,ny);
    presstemp = zeros(nx,ny,ntime-1);
%    presstemp = zeros(nx,ny,size(press,2));
    for k = 1:numel(press)/nx/ny
      presstemp(:,:,k) = reshape(press(:,k),nx,ny);
    end
    press = presstemp;
  end
end

options.savefig = ~false;
if options.savefig
  count = 1;
  !mkdir ./frames
end

if options.stress
  fileName = ['output/' file(1:end-8) 'Stress11.bin'];
  fid = fopen(fileName,'r');
  s11 = fread(fid,'double');
%  nStress = numel(s11)/(ntime-1+2);
  nStress = nx*ny;
  stress11X = s11(1:nStress);
  stress11Y = s11(nStress+1:2*nStress);
  stress11 = zeros(nStress,ntime-1);
  for k = 1:ntime - 1
    istart = (k+1)*nStress + 1;
    iend = istart + nStress - 1;
    if iend <= numel(s11)
      stress11(:,k) = s11(istart:iend);
    end
  end
  if 1
    stress11X = reshape(stress11X,nx,ny);
    stress11Y = reshape(stress11Y,nx,ny);
    stresstemp = zeros(nx,ny,ntime-1);
    for k = 1:numel(stress11)/nx/ny
      stresstemp(:,:,k) = reshape(stress11(:,k),nx,ny);
    end
    stress11 = stresstemp;
  end


  fileName = ['output/' file(1:end-8) 'Stress12.bin'];
  fid = fopen(fileName,'r');
  s12 = fread(fid,'double');
%  nStress = numel(s12)/(ntime-1+2);
  nStress = nx*ny;
  stress12X = s12(1:nStress);
  stress12Y = s12(nStress+1:2*nStress);
  stress12 = zeros(nStress,ntime-1);
  for k = 1:ntime - 1
    istart = (k+1)*nStress + 1;
    iend = istart + nStress - 1;
    if iend <= numel(s12)
      stress12(:,k) = s12(istart:iend);
    end
  end
  if 1
    stress12X = reshape(stress12X,nx,ny);
    stress12Y = reshape(stress12Y,nx,ny);
    stresstemp = zeros(nx,ny,ntime-1);
    for k = 1:numel(stress12)/nx/ny
      stresstemp(:,:,k) = reshape(stress12(:,k),nx,ny);
    end
    stress12 = stresstemp;
  end


  fileName = ['output/' file(1:end-8) 'Stress21.bin'];
  fid = fopen(fileName,'r');
  s21 = fread(fid,'double');
%  nStress = numel(s21)/(ntime-1+2);
  nStress = nx*ny;
  stress21X = s21(1:nStress);
  stress21Y = s21(nStress+1:2*nStress);
  stress21 = zeros(nStress,ntime-1);
  for k = 1:ntime - 1
    istart = (k+1)*nStress + 1;
    iend = istart + nStress - 1;
    if iend <= numel(s21)
      stress21(:,k) = s21(istart:iend);
    end
  end
  if 1
    stress21X = reshape(stress21X,nx,ny);
    stress21Y = reshape(stress21Y,nx,ny);
    stresstemp = zeros(nx,ny,ntime-1);
    for k = 1:numel(stress21)/nx/ny
      stresstemp(:,:,k) = reshape(stress21(:,k),nx,ny);
    end
    stress21 = stresstemp;
  end


  fileName = ['output/' file(1:end-8) 'Stress22.bin'];
  fid = fopen(fileName,'r');
  s22 = fread(fid,'double');
%  nStress = numel(s22)/(ntime-1+2);
  nStress = nx*ny;
  stress22X = s22(1:nStress);
  stress22Y = s22(nStress+1:2*nStress);
  stress22 = zeros(nStress,ntime-1);
  for k = 1:ntime - 1
    istart = (k+1)*nStress + 1;
    iend = istart + nStress - 1;
    if iend <= numel(s22)
      stress22(:,k) = s22(istart:iend);
    end
  end
  if 1
    stress22X = reshape(stress22X,nx,ny);
    stress22Y = reshape(stress22Y,nx,ny);
    stresstemp = zeros(nx,ny,ntime-1);
    for k = 1:numel(stress22)/nx/ny
      stresstemp(:,:,k) = reshape(stress22(:,k),nx,ny);
    end
    stress22 = stresstemp;
  end
end


if options.dist
  dist = [];
  ind = [];
end

figure(1); clf
for k = 1:irate:ntime
  xx = interpft(posx(:,:,k),96); yy = interpft(posy(:,:,k),96);  
  vec1 = [xx(:,:);xx(1,:)];
  vec2 = [yy(:,:);yy(1,:)];
  if options.center
    vec1 = vec1 - 1/2*(min(vec1) + max(vec1));
  end
  figure(1); clf
  if options.tracers
    hold on
%    plot(tracerx(1:1101,k),tracery(1:1101,k),'bo')
%    plot(tracerx(1102:2386,k),tracery(1102:2386,k),'go')
%    plot(tracerx(2387:3855,k),tracery(2387:3855,k),'co')
%    plot(tracerx(3856:5492,k),tracery(3856:5492,k),'yo')
%    plot(tracerx(5493:7269,k),tracery(5493:7269,k),'mo')
    plot(tracerx(:,k),tracery(:,k),'bo')
    plot(tracerx(1:1110,k),tracery(1:1110,k),'bo')
    plot(tracerx(1111:2420,k),tracery(1111:2420,k),'go')
    plot(tracerx(2421:3930,k),tracery(2421:3930,k),'co')
    plot(tracerx(3931:5640,k),tracery(3931:5640,k),'yo')
    plot(tracerx(5641:7550,k),tracery(5641:7550,k),'mo')
  end
  plot(vec1,vec2,'r-','linewidth',2)

  if options.center
    hold on
    plot(50*[-2 2],[0 0],'k')
    plot([0 0],50*[-2 2],'k')
  end
  hold on
%  plot(vec1(1,:),vec2(1,:),'r.','markersize',30)
  plot(vec1(:,:),vec2(:,:),'r','markersize',30)
  axis equal
  axis(ax);
%  set(gca,'visible','off')
  titleStr = ['t = ' num2str(time(k),'%4.2e') ...
      ' eA = ' num2str(ea(k),'%4.2e') ...
      ' eL = ' num2str(el(k),'%4.2e')];
  title(titleStr)

  if options.dist
    vec1 = posx(:,2,k) + 1i*posy(:,2,k);
    vec2 = posx(:,3,k) + 1i*posy(:,3,k);

    ddist = 10000;
    for j = 1:numel(vec1)
      ddist = min(ddist,min(abs(vec1(j)-vec2)));
    end
    dist = [dist ddist];
  end

  if (options.quiver && k+1 < ntime)
    velx = (posx(:,:,k+1) - posx(:,:,k))/(time(k+1)-time(k));
    vely = (posy(:,:,k+1) - posy(:,:,k))/(time(k+1)-time(k));
    hold on
%    quiver(posx(:,:,k),posy(:,:,k),velx,vely,'k')
    quiver(vec1(1:end-1,:),vec2(1:end-1,:),velx,vely,'k')
    hold off
  end
  % Quiver plot of a finite difference applied to position
  
  if options.confined
    hold on
    vec1 = [wallx;wallx(1,:)];
    vec2 = [wally;wally(1,:)];
    plot(vec1,vec2,'k','linewidth',3);
    hold off
  end
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');

  if (options.pressure && k >= 2)
%    figure(1); hold on
%    plot([xmin xmin],[ymin ymax],'k')
%    plot([xmax xmax],[ymin ymax],'k')
%    plot([xmin xmax],[ymin ymin],'k')
%    plot([xmin xmax],[ymax ymax],'k')
%    hold off;

    figure(2); clf
    surf(pressX,pressY,press(:,:,k-1))
    shading interp
    view(2)
    axis equal
%    axis(ax);
    colorbar
%    caxis([cmin cmax])
    if confined
      hold on
      vec1 = [wallx;wallx(1,:)];
      vec2 = [wally;wally(1,:)];
      plot3(vec1,vec2,10*ones(size(vec1)),'k','linewidth',3);
      hold off
    end
    set(gca,'visible','off')
  end


  if options.savefig
    filename = ['./frames/image', sprintf('%04d',count),'.png'];
    count = count+1;
    figure(1);
    print(gcf,'-dpng','-r300',filename);
  end
  pause(0.01)
end

if (options.angle || options.jacobian)
  IA = zeros(ntime,1);
  jacobian = zeros(n,ntime);
  % compute inclination angle on an upsampled grid
  N = 1024;
  modes = [(0:N/2-1)';0;(-N/2+1:-1)'];

  for k = 1:ntime
    x = interpft(posx(:,:,k),N); y = interpft(posy(:,:,k),N);
    Dx = real(ifft(1i*modes.*fft(x)));
    Dy = real(ifft(1i*modes.*fft(y)));
    jac = sqrt(Dx.^2 + Dy.^2);
    tx = Dx./jac; ty = Dy./jac;
    nx = ty; ny = -tx;
    rdotn = x.*nx + y.*ny;
    rho2 = x.^2 + y.^2;

    J11 = 0.25*sum(rdotn.*(rho2 - x.*x).*jac)*2*pi/N;
    J12 = 0.25*sum(rdotn.*(-x.*y).*jac)*2*pi/N;
    J21 = 0.25*sum(rdotn.*(-y.*x).*jac)*2*pi/N;
    J22 = 0.25*sum(rdotn.*(rho2 - y.*y).*jac)*2*pi/N;

    J = [J11 J12; J21 J22];
    [V,D] = eig(J);

    [~,ind] = min(abs(diag(D)));
    IA(k) = atan(V(2,ind)/V(1,ind));
    jacobian(:,k) = interpft(jac,n);
  end
  IA(1) = pi/2;
end



