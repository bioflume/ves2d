

% Example implementation of memory mapping an NPY file using readNPYheader

filename = '/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/near_vel_fft_models/normParams/in_param_mode113_128.npy';

[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

figure;

if fortranOrder
    f = memmapfile(filename, 'Format', {dataType, arrayShape, 'd'}, 'Offset', totalHeaderLength);
    image(f.Data.d)
    
else
    % Note! In this case, the dimensions of the array will be transposed,
    % e.g. an AxBxCxD array becomes DxCxBxA. 
    f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
    
    tmp = f.Data.d;
    img = permute(tmp, length(arrayShape):-1:1); % note here you have to reverse the dimensions. 
    image(img./255)
end


%%
clear;

addpath ../../learnVes/dataPrepCodes/
inFile = '2024Oct_selften_input_downsample32.npy';
outFile = '2024Oct_selften_out_downsample32.npy';

% filename = '/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/near_vel_allModes_normParams/out_param_allmode.npy';
filename = inFile;
% filename = '~/Desktop/near_trained/out_param_downsample32_allmode.npy';
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
input = f.Data.d;

filename = outFile;
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
output = f.Data.d;

load ~/Documents/GitHub/ves2d/learnVes/dataPrepCodes/selfTensionDataSet.mat
myInput = [interpft(Xinput(1:end/2,:),32);interpft(Xinput(end/2+1:end,:),32)];
myOutput = interpft(Toutput,32);

for it = 1 : 156225
  figure(1);clf;
  subplot(1,2,1)
  plot(myInput(1:end/2,it),myInput(end/2+1:end,it),'k')
  hold on
  plot(input(:,1,it),input(:,2,it),'r')
  axis equal
  title(it)

  subplot(1,2,2)
  plot(linspace(0,1,32)',myOutput(:,it),'k')
  hold on
  plot(linspace(0,1,32)',output(:,it),'r')
  axis square

  
  pause


end


%%
% clear;
% in_param = zeros(32,4);
in_param = zeros(31,4);

% filename = '/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/near_vel_allModes_normParams/in_param_allmode.npy';
filename = '/Users/gokberk/Desktop/adv_trained/2024Oct_advfft_in_para_downsample_all_mode.npy';
% filename = '~/Desktop/near_trained/in_param_downsample32_allmode.npy';
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
tmp = f.Data.d;

for k = 1 : 31

in_param(k,:) = tmp(:,k)';

end


%%
clear;
out_param = zeros(32,2,12);

filename = '/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/near_vel_32modesfft_models/out_param_downsample32_allmode.npy';
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
tmp = f.Data.d;

for k = 1 : 32

out_param(k,:,:) = tmp(:,:,k)';

end




%%
% clear;
in_param = zeros(32,4);

filename = '/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/near_vel_32modesfft_models/in_param_downsample32_allmode.npy';
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
tmp = f.Data.d;

for k = 1 : 32

in_param(k,:) = tmp(:,k)';

end

%%
% clear;


filename = './TG_from_shanBIEM_-50.npy';
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);

tmp = f.Data.d;

vesx = zeros(128,128);
vesy = zeros(128,128);
X = zeros(256,128);
for k = 1 : 128
X(:,k) = tmp(k,:);
% X = X';
% vesx(:,:,k) = X(1:end/2,:);
% vesy(:,:,k) = X(end/2+1:end,:);
end

%% 
cmap = colormap('hsv');
cx = mean(vesx,1); cx = reshape(cx,128,100); 
cy = mean(vesy,1); cy = reshape(cy,128,100);
dcx = cx(:,2:end) - cx(:,1:end-1); dcx = [zeros(128,1) dcx];
dcy = cy(:,2:end) - cy(:,1:end-1); dcy = [zeros(128,1) dcy];

dr = sqrt(dcx.^2 + dcy.^2); 
for k = 1 : 100
[pdf_dr,xi] = ksdensity(dr(:,k));
figure(1);clf;
histogram(dr(:,k))
hold on
% plot(xi,pdf_dr,'k','linewidth',2)
xlim([0 0.05])
ylim([0 60])
% ylim([0 0.1])
axis square
title(k)
pause(0.1)
ax = gca;
exportgraphics(ax,['~/Desktop/figs/crashF' num2str(k) '.png'],'Resolution',300)
end
% for k = 30 : 100
% figure(1);clf; hold on;
% for ives = 1 : 128
%   plot(vesx(:,ives,k),vesy(:,ives,k),'Color',cmap(2*ives,:),'linewidth',2)
%   quiver(cx(ives,k),cy(ives,k),50*dcx(ives,k),50*dcy(ives,k),'Color',cmap(2*ives,:),'AutoScale','off')
% end
% axis equal
% xlim([-0.5 2.5])
% ylim([0 2])
% title(k)
% 
% ax = gca;
% exportgraphics(ax,['~/Desktop/figs/crashF' num2str(k) '.png'],'Resolution',300)
% pause(0.1)
% end
