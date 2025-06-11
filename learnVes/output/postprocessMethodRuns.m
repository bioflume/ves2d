clear;
addpath ../../src/

prenames{1} = './N32_implicit_hedgehog.bin';
prenames{2} = './N32_implicit_rbf.bin';
prenames{3} = './N32_implicit_noNear.bin';
prenames{4} = './N32_explicit_hedgehog.bin';
prenames{5} = './N32_explicit_rbf.bin';
prenames{6} = './N32_explicit_noNear.bin';

names{1} = './res_N32_implicit_hedgehog.bin';
names{2} = './res_N32_implicit_rbf.bin';
names{3} = './res_N32_implicit_noNear.bin';
names{4} = './res_N32_explicit_hedgehog.bin';
names{5} = './res_N32_explicit_rbf.bin';
names{6} = './res_N32_explicit_noNear.bin';

op = poten(32);
% nCollidVes = zeros(50,6);
% for irun = 1 : 6
%   [vesx, vesy, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(names{irun});
%   X = [vesx(:,:,end);vesy(:,:,end)];
%   save(['finalConfigRun' num2str(irun) '.mat'],'X')
% end
for irun = 1 : 6
  [vesx1, vesy1, ten, time1, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(prenames{irun});

  [vesx2, vesy2, ten, time2, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(names{irun});
  vesx = zeros(NN,nv,102);
  vesy = zeros(NN,nv,102);
  for it = 1 : 52
  vesx(:,:,it) = vesx1(:,:,it);
  vesy(:,:,it) = vesy1(:,:,it);
  end
  for it = 1 : 50
  vesx(:,:,52+it) = vesx2(:,:,it+2);
  vesy(:,:,52+it) = vesy2(:,:,it+2);
  end
  
  time = [time1; time1(end) + time2(3:end)];

  if 0
  for istep = 1 : numel(time)
    disp(['Collision Check step: ' num2str(istep)])
    vesicle = capsules([vesx(:,:,istep);vesy(:,:,istep)],[],[],1,1,0);
    NearV2V = vesicle.getZone(vesicle,1);
    [icollisionVes,collidingIdcs,safeIdcs] = vesicle.collisionWVesOutput(NearV2V,0,op);
    nCollidVes(istep,irun) = numel(collidingIdcs);
  end
  end

  if 1
  numberOfFrames = numel(time);
  hFigure = figure;
  allTheFrames = cell(numberOfFrames,1);
  vidHeight = 576;
  vidWidth = 1024;
  allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
  allTheColorMaps = cell(numberOfFrames,1);
  allTheColorMaps(:) = {zeros(256, 3)};
  set(gcf, 'renderer', 'zbuffer');

  Vsize = 5;
  [xx,yy] = meshgrid(linspace(-0.5,5.5,50)',linspace(-0.5,5.5,50)');
  uu = sin(xx/Vsize*pi).*cos(yy/Vsize*pi); 
  vv = -cos(xx/Vsize*pi).*sin(yy/Vsize*pi);

  frameCount = 1;
  for k = 1 : 1 : numel(time)
    cla reset;
    kT = k;
    xvecT = [vesx(:,:,kT);vesx(1,:,kT)] ;
    yvecT = [vesy(:,:,kT);vesy(1,:,kT)];
  

    figure(1); clf; 
    h = plot(xvecT, yvecT, 'Color',[26/255 150/255 65/255 1],'linewidth',2);
    hold on
    hFill = fill(xvecT, yvecT,[26/255 150/255 65/255]);
    set(hFill,'EdgeColor',[26/255 150/255 65/255])
    if irun == 6 && k > 37
    h = plot(xvecT(:,[62;70]), yvecT(:,[62;70]), 'r','linewidth',2);
    hold on
    hFill = fill(xvecT(:,[62;70]), yvecT(:,[62;70]),'r');
    set(hFill,'EdgeColor','r')
    end
    
    hold on
    l = streamslice(xx,yy,uu,vv);
    set(l,'Color',[12/255,44/255,132/255, 0.75])
    set(l,'linewidth',0.5)
 
    plot(-0.5, -0.5, 'k.','markersize',0.001)
    plot(5, 5, 'k.','markersize',0.001)

    axis equal

    xlim([0 5])
    ylim([0 5])

    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'ztick',[]);

    title(kT)
    set(gca,'xcolor','w');
    set(gca,'ycolor','w');
    set(gca,'zcolor','w');
    box on
    set(gca,'visible','off')

 
 
    drawnow;
    myMovie(frameCount) = getframe(gca);
    frameCount = frameCount + 1;
  end

  fullFileName = ['~/Desktop/' names{irun} '.avi'];
  profile = 'Uncompressed AVI';
  writerObj = VideoWriter(fullFileName, profile);
  open(writerObj);
  % Write out all the frames.
  numberOfFrames = length(myMovie);
  for frameNumber = 1 : numberOfFrames 
     writeVideo(writerObj, myMovie(frameNumber));
  end
  close(writerObj);
  end
end