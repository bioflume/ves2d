

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

out_param = zeros(128,2,12);

filename = '/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/near_vel_fft_models/normParams/out_param_mode1_16.npy';
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
tmp = f.Data.d;

for k = 1 : 16

out_param(k,:,:) = tmp(:,:,k)';

end


filename = '/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/near_vel_fft_models/normParams/out_param_mode17.npy';
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
tmp = f.Data.d;

out_param(17,:,:) = tmp(:,:,1)';


filename = '/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/near_vel_fft_models/normParams/out_param_mode113_128.npy';
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);

f = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
tmp = f.Data.d;

for k = 1 : 16
out_param(112+k,:,:) = tmp(:,:,k)';
end