radScale = [2;4;8;16;32;64];
Ns = [16;32;64;128];

[nn,rr] = meshgrid(Ns,radScale);
errors = zeros(size(nn));


% for k = 1 : numel(nn(:))
% errors(k) = translationalTestAlpert(nn(k),rr(k));
% end

for k = 1 : numel(nn(:))
errors(k) = rotationalTestAlpert(nn(k),rr(k));
end


