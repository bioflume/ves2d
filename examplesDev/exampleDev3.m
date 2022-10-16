clear functions global;

% Get the initial marker point locations
n = 64;                                 % Number of discretization points
X = boundary(n,'curly');


% Setting modeling parameters
prams.T = 20;                            % Simulation time
prams.m = 800;                          % Number of time steps
prams.kappa = 1e-1;                     % Bending modulus
prams.order = 1;                        % order of the slover
prams.ts = prams.T/prams.m;             % step size (time)
prams.Incompressibility = 1;            % Incompressibility
prams.vInf = @(X) farFieldVel(X,'shear',0);
                                        % The far field velocity

% Setting Options
options.usePlot = 1;                    % To plot the intermidate states   
options.AxesHandle = gca;               % The axes for the  plot
options.axisOn = 0;                     % To show or hide plot axis
options.ProgressBar = 1;
options.saveData = 0;
options.dataStride = 17;

% Calling update function with all of the above parameters 
[Xfinal status] = Ves2D(X,prams,options);
% 
% fileId = fopen(status.fileName,'r');
% Result = fread(fileId,'double');
% fclose(fileId);
% 
% Result = reshape(Result,321,[]);       %The size of the saved data in
%                                        %each step (here 321) is included
%                                        %in the name of the file for
%                                        %future reference. The data is
%                                        %saved as a vector in the order of 
%                                        %[POAITION;TENSION;VELOCITY;TIME].
%                                        %The size of this vector is 5*n+1.
% 
% X = Result(1:2*n,:);
% sigma = Result(2*n+1:3*n,:);
% u = Result(3*n+1:5*n,:);
% 
% %Plotting the Results is one figure
% figure;
% for ii=1:12
%   subplot(3,4,ii); plot(Result([1:n 1],ii),Result([n+1:2*n n+1],ii), ...
%                         'Linewidth',2);
%      title(['t = ' num2str((ii-1)*17*prams.ts) 's']);
%      axis([-1.5 1.5 -2 2]);
%      if(mod(ii,4)~=1), set(gca,'ytick',[]);end
%      if(ii<9), set(gca,'xtick',[]);end
%  end    
% 
% % Plotting the bending enegry
% e(1) = Energy(X(:,1));
% for ii=1:size(Result,2)
%   e(ii+1) = Energy(X(:,ii));
% end
% figure; loglog(0:ii,e)
% title('Bending energy of the relaxing vesicle (log-log)');
% xlabel('Time step');ylabel('e')

% clear functions global
% 
% 
% qw = quadratureS(32,8,2*pi);
% for ii = 1:size(X,2)
% XX = reshape(X(:,ii),[],2);
% SS = sigma(:,end);
% uu = reshape(u(:,ii),[],2);
% 
% sa = sqrt(D1FT(XX(:,1)).^2 + D1FT(XX(:,2)).^2);
% x1 = DmFT(XX,1,[],1./sa);
% fk = -prams.kappa*DmFT(x1,3); X
% fs = DmFT(diag(sigmaN)*x1,1);
% f = fk+fs;
% 
% G = kernelS(XX(:),qw);
% uf = reshape(G*f(:),[],2);
% e = uu-uf;
% 
% disp(max(dot(e,e,2)));
% subplot(1,2,1); plot(XX(:,1),XX(:,2));
% hold on; quiver(XX(:,1),XX(:,2),uu(:,1),uu(:,2),0); hold off
% 
% subplot(1,2,2); plot(XX(:,1),XX(:,2));
% hold on; quiver(XX(:,1),XX(:,2),uf(:,1),uf(:,2),0); hold off
% 
% pause(.2);
% end