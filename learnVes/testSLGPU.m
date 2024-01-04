clear; clc;

load ./ICs/VF35_81VesIC.mat; % 81 vesicles loaded
% X loaded, N = 64 
nv = size(X,2);
Ntests = [32;64;128;256;512];
sys_size = Ntests*81;

addpath ../src/
oc = curve;

% CPU and FMM (how many nodes?) tests are done on 

gpuTimesNew = zeros(numel(Ntests),1);
gpuTimesOld = zeros(numel(Ntests),1);



for it = 1 : numel(Ntests)
  disp(it)
  % Choose discretization
  N = Ntests(it);
  % upsample the vesicle resolution
  Xv = [interpft(X(1:end/2,:),N); interpft(X(end/2+1:end,:),N)];
  % build vesicle object
  vesicle = capsules(Xv, [], [], 1, 1, 0);
  % calculate bending force
  fBend = vesicle.tracJump(Xv,zeros(N,nv));
 
  
  % GPU Calculation
%   indexMat = {};
%   for k = 1 : vesicle.nv
%     indexMat{k} = ones(vesicle.N)*nan;
%   end
%   blkIndMat = gpuArray((single(blkdiag(indexMat{:}))));
%   saGPU = gpuArray(single(vesicle.sa));
%   fxGPU = gpuArray(single(fBend(1:end/2,:)));
%   fyGPU = gpuArray(single(fBend(end/2+1:end,:)));
%   xGPU = gpuArray(single(vesicle.X(1:end/2,:)));
%   yGPU = gpuArray(single(vesicle.X(end/2+1:end,:)));
%   stokesXGPUarr = gpuArray(single(zeros(vesicle.N*vesicle.nv,1)));
%   stokesYGPUarr = gpuArray(single(zeros(vesicle.N*vesicle.nv,1)));
%   stokesGPU = @() exactStokesSLGPUFastVect(xGPU, yGPU, saGPU, fxGPU, fyGPU, stokesXGPUarr, stokesYGPUarr, blkIndMat);
%   gpuTimesNew(it) = gputimeit(stokesGPU);
  

  rows = zeros(vesicle.N^2,vesicle.nv);
  cols = zeros(vesicle.N^2,vesicle.nv);
  for k = 1 : vesicle.nv
    [rr,cc] = (meshgrid((k-1)*N+1:k*N, (k-1)*N+1:k*N));
    rows(:,k) = rr(:);
    cols(:,k) = cc(:);
  end
  
  nrows = vesicle.N*vesicle.nv;
  indices = gpuArray(single((cols(:)-1)*nrows + rows(:))); 
  saGPU = gpuArray(single(vesicle.sa));
  fxGPU = gpuArray(single(fBend(1:end/2,:)));
  fyGPU = gpuArray(single(fBend(end/2+1:end,:)));
  xGPU = gpuArray(single(vesicle.X(1:end/2,:)));
  yGPU = gpuArray(single(vesicle.X(end/2+1:end,:)));
  stokesXGPUarr = gpuArray(single(zeros(vesicle.N*vesicle.nv,1)));
  stokesYGPUarr = gpuArray(single(zeros(vesicle.N*vesicle.nv,1)));
  stokesGPU = @() exactStokesSLGPUFastVectMemFriend(xGPU, yGPU, saGPU, fxGPU, fyGPU, stokesXGPUarr, stokesYGPUarr, indices);
  gpuTimesNew(it) = gputimeit(stokesGPU);
  %[stokesXGPUarr, stokesYGPUarr] = exactStokesSLGPUFastVectMemFriend(xGPU, yGPU, saGPU, fxGPU, fyGPU, stokesXGPUarr, stokesYGPUarr, indices);
%   stokesXGPUarr = reshape(stokesXGPUarr,[vesicle.N vesicle.nv]);
%   stokesYGPUarr = reshape(stokesYGPUarr,[vesicle.N vesicle.nv]);
  

  %stokesSLP = zeros(2*vesicle.N,vesicle.nv);
  %stokesSLP = exactStokesSLGPUVect(vesicle.X, vesicle.sa, fBend, stokesSLP);
  
  saGPU = gpuArray(single(vesicle.sa));
  fGPU = gpuArray(single(fBend));
  XGPU = gpuArray(single(vesicle.X));
  stokesGPUarr = gpuArray(single(zeros(2*vesicle.N,vesicle.nv)));
  stokesGPU = @() exactStokesSLGPUVect(XGPU, saGPU, fGPU, stokesGPUarr);
  gpuTimesOld(it) = gputimeit(stokesGPU);

end
save('comparison.mat','gpuTimesNew')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesXSLP, stokesYSLP] = exactStokesSLGPUFastVect(x, y, sa, fx, fy, stokesXSLP, stokesYSLP,blkIndMat)
N = size(x,1);
nv = size(x,2);
Nall = N*nv;

% multiply by arclength term
% do not even need this for loop
% all to all calculation, then apply arrayfun to remove itself
denx = fx.*sa*2*pi/N;
deny = fy.*sa*2*pi/N;
denx = denx(:); deny = deny(:);
ddenx = denx(:,ones(Nall,1))';
ddeny = deny(:,ones(Nall,1))';

x = x(:); y = y(:);
xsou = x(:,ones(Nall,1))';
ysou = y(:,ones(Nall,1))';
diffx = xsou'-xsou; %xtar = xsou'
diffy = ysou'-ysou;

dis2 = diffx.^2 + diffy.^2;

coeff = 0.5*log(dis2);
val = coeff.*ddenx;
stokesXSLP = stokesXSLP - sum(val+blkIndMat,2,"omitnan");
val = coeff.*ddeny;
stokesYSLP = stokesYSLP - sum(val+blkIndMat,2,"omitnan");

coeff = (diffx.*ddenx + diffy.*ddeny)./dis2;
val = coeff.*diffx;
stokesXSLP = stokesXSLP + sum(val+blkIndMat,2,"omitnan");
val = coeff.*diffy;
stokesYSLP = stokesYSLP + sum(val+blkIndMat,2,"omitnan"); 


end % exactStokesSLGPUFastVect

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesXSLP, stokesYSLP] = exactStokesSLGPUFastVectMemFriend(x, y, sa, fx, fy, stokesXSLP, stokesYSLP,indices)
N = size(x,1);
nv = size(x,2);
Nall = N*nv;

% multiply by arclength term
% do not even need this for loop
% all to all calculation, then apply arrayfun to remove itself
denx = fx.*sa*2*pi/N;
deny = fy.*sa*2*pi/N;
denx = denx(:); deny = deny(:);

x = x(:); y = y(:);
diffx = x(:,ones(Nall,1))-x(:,ones(Nall,1))'; %xtar = xsou'
diffy = y(:,ones(Nall,1))-y(:,ones(Nall,1))';

dis2 = diffx.^2 + diffy.^2;

coeff = 0.5*log(dis2);
coeff(indices) = 0;
stokesXSLP = stokesXSLP - coeff*denx;
stokesYSLP = stokesYSLP - coeff*deny;

coeff1 = (diffx.^2 ./ dis2);
coeff1(indices) = 0;
coeff2 = (diffx.*diffy)./dis2 ;
coeff2(indices) = 0;
stokesXSLP = stokesXSLP + coeff1*denx + coeff2*deny;

coeff1 = (diffx.*diffy)./dis2;
coeff1(indices) = 0;
coeff2 = (diffy.^2 ./ dis2);
coeff2(indices) = 0;
stokesYSLP = stokesYSLP + coeff1*denx + coeff2*deny; 


end % exactStokesSLGPUFastVectMemFriend

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stokesSLP = exactStokesSLGPUVect(X, sa, f, stokesSLP)
N = size(X,1)/2;
nv = size(X,2);
Ntar = N;
fx = f(1:end/2,:);
fy = f(end/2+1:end,:);

% multiply by arclength term

for k = 1:nv % vesicle of targets
  K = [(1:k-1) (k+1:nv)];
  % Loop over all vesicles except k
  allX = X(1:N,K); allX = allX(:);
  allY = X(N+1:2*N,K); allY = allY(:);
  allFx = fx(:,K);
  allFy = fy(:,K);
  allsa = sa(:,K);

  denx = allFx.*allsa*2*pi/N; denx = denx(:);
  deny = allFy.*allsa*2*pi/N; deny = deny(:);
  
  ddenx = denx(:,ones(Ntar,1))';
  ddeny = deny(:,ones(Ntar,1))';
  
  xsou = allX(:,ones(Ntar,1))';
  ysou = allY(:,ones(Ntar,1))';

  Nsou = numel(allX);
  xx = X(1:N,k);
  yy = X(N+1:2*N,k);
  xtar = xx(:,ones(Nsou,1));
  ytar = yy(:,ones(Nsou,1));
  
  diffx = xtar-xsou;
  diffy = ytar-ysou;
  dis2 = diffx.^2 + diffy.^2;
  
  coeff = 0.5*log(dis2);
  
  val = coeff.*ddenx;
  stokesSLP(1:N,k) = -sum(val,2);
  val = coeff.*ddeny;
  stokesSLP(N+1:2*N,k) = -sum(val,2);

  coeff = (diffx.*ddenx + diffy.*ddeny)./dis2;

  val = coeff.*diffx;
  stokesSLP(1:N,k) = stokesSLP(1:N,k) + sum(val,2);
  val = coeff.*diffy;
  stokesSLP(N+1:2*N,k) = stokesSLP(N+1:2*N,k) + sum(val,2); 
end % k

end % exactStokesSLGPUVect

