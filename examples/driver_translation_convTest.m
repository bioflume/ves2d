%%
radScale = [2;4;8;16;32;64];
Ns = [16;32;64;128];
[nn,rr] = meshgrid(Ns,radScale);
for k = 1 : numel(nn(:))
    errors(k) = translationalTestSymmAlpert(nn(k),rr(k));
    pause
end


%%
radScale = [2;4;8;16;32;64];
Ns = [16;32;64;128];
[nn,rr] = meshgrid(Ns,radScale);
for k = 1 : numel(nn(:))
    errors(k) = rotationalTestSymmAlpert(nn(k),rr(k));
end


%% 
Ns = [16;32;64;128];
errors = zeros(4,3);
for k = 1 : numel(Ns)
    errors(k,:) = singleLayerMatrixTest(Ns(k));
end


%%
radScale = [2;4;8;16;32];
Ns = [16;32;64];

[nn,rr] = meshgrid(Ns,radScale);
errors = zeros(size(nn));
for k = 1 : numel(nn)
  errors(k,:) = ellipsoid_translationalTestSymmAlpert(nn(k),rr(k));
  pause
end

