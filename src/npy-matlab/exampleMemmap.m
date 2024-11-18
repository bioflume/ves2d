

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
% out_param = zeros(32,2,12);
out_param = zeros(127,4);

% filename = '/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/near_vel_allModes_normParams/out_param_allmode.npy';
filename = './2024Oct_advten_out_para_allmodes.npy';
% filename = '~/Desktop/near_trained/out_param_downsample32_allmode.npy';
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
tmp = f.Data.d;

for k = 1 : 127

out_param(k,:) = tmp(:,k)';

end




%%
% clear;
% in_param = zeros(32,4);
in_param = zeros(127,4);

% filename = '/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/near_vel_allModes_normParams/in_param_allmode.npy';
filename = './2024Oct_advten_in_para_allmodes.npy';
% filename = '~/Desktop/near_trained/in_param_downsample32_allmode.npy';
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
tmp = f.Data.d;

for k = 1 : 127

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


filename = './save_intm_var_Nov.npy';
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
