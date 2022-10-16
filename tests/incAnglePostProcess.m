addpath ../src/
% clear;

% ReducedArea = [.6 .7 .8 .9];
% ViscCont = [.25 .5 1 2 4 8 16];
% ShearRate = [.01 .1 1 10];
% np = 32;

% %- post process
% for ii = 1:length(ShearRate)
%   for jj = 1:length(ReducedArea)
%     for kk = 1:length(ViscCont)
%       sr = ShearRate(ii);
%       ra = ReducedArea(jj);
%       vc = ViscCont(kk);

%       fileName = ['InclinationAngle_' num2str(sr) '_' num2str(ra) '_' ...
%                   num2str(vc)];
      
%       if ( exist([fileName '.mat']) )
%         runFile = load([fileName '.mat']);
%         fileId = fopen([fileName '.bin'],'r');
%         Result = fread(fileId,'double');
%         fclose(fileId);

%         nv = 1;
%         Result = reshape(Result,5*nv*np+1,[]); 
%         Xv   = Result(1       :2*nv*np,:);
%         Time = Result(5*nv*np+1       ,:);
%         inclinationAngle = [];

%         for ll=1:size(Xv,2)
%           II = secondMoment(Xv(:,ll));
%           e = min(eig(II));
%           inclinationAngle(ll) = atan(II(1,2)/(II(2,2)-e));
%         end

%         data(ii,jj,kk).ShearRate = ShearRate(ii);
%         data(ii,jj,kk).ReducedArea = runFile.Ra;
%         data(ii,jj,kk).ViscCont = ViscCont(kk);
%         data(ii,jj,kk).Perimeter = runFile.L;
%         data(ii,jj,kk).Kappa = runFile.prams.kappa;
%         data(ii,jj,kk).IncAngle = inclinationAngle;
%         data(ii,jj,kk).Time = Time;
% %        data(ii,jj,kk).Xv = Xv;
        
%         disp([ii jj kk]);
% %         plot(Time,inclinationAngle);
% %         drawnow; 
% %         title(['Shear rate=' num2str(sr) ', reduced area=' num2str(ra) ...
% %                ', viscosity contrast=' num2str(vc)]);
% %         disp('Pausing');
% %         pause;
%       end
%       %[Ra A L] = reducedVolume(X);
%       %NonDimShear = 1 * (L / 2 / pi)^3 / prams.kappa;

%     end
%   end
% end

for ii=1:size(data,1)
  for jj=1:size(data,2)
    for kk=1:size(data,3)
      idxSet = {ii,jj,kk};
      disp('=============================');
      disp(data(idxSet{:}));
      disp('=============================');
      plot(data(idxSet{:}).Time,data(idxSet{: }).IncAngle);
      drawnow;
      pause;
    end 
  end
end

