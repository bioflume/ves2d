
clear all;

% % fprintf('   Vel       Tracion     Pressure\n');
% % fprintf('--------------------------------------\n');
% % X = [];
% % 
% % flow = 'couette';
% % Ri = 5; Ro = 15;c = Ri*Ro/(Ro^2-Ri^2);
% % 
% % TRAC = @(x,y) -2*c*Ro/(x^2+y^2)*[0 1;1 0];
% % PRES = @(x,y) 0*x;
% % 
% % for a = 5:8
% %     prams.M = 2^a*[1 1]; np = sum(prams.M);
% %     prams.bd = @(ind,m) sampleBd(ind,m,1,'couette');
% % 
% %     prams.bc = @(x) forcing(x,flow);
% %     domain = fixedBound(prams.M,prams.bd,1);
% % 
% %     while(size(X,1)~=128)
% %         xrand = 30*rand(1000,2)-15;
% %         IN = inDomain(domain,xrand(:,1),xrand(:,2),.1);
% %         X = [X;xrand(IN,:)];
% %         X = X(1:min(size(X,1),128),:);
% %     end
% % 
% %     XX = []; n = [];
% %     for ind = 1:length(domain)
% %         XX = [XX;domain(ind).X];
% %         n = [n;domain(ind).n];
% %     end
% %     R = XX./[myNorm(XX) myNorm(XX)];
% % 
% %     [U trash mu] = farFieldVel(X,[],'confined',prams);
% %     U = reshape(U,[],2);
% %     Uref = forcing(X,flow);
% % 
% %     [st p] = evalTraction(prams.M,domain,mu);
% %     stPol = []; stPol(:,1) = dot(st,R,2);
% %     stPol(:,2) = dot(st,[-R(:,2) R(:,1)],2);
% %     %avgS2 = avgStress(prams.M,domain,mu,[],[]);
% %     for ind = 1:size(XX,1)
% %         stEx(ind,:) = sign(dot(n(ind,:),R(ind,:)))*TRAC(XX(ind,1),XX(ind,2))*[1;0];
% %         pEx(ind,1)  = PRES(XX(ind,1),XX(ind,2));
% %     end
% % 
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %viewer(domain);
% %     %hold on; quiver(X(:,1),X(:,2),U(:,1),U(:,2));
% %     %hold on; quiver(X(:,1),X(:,2),Uref(:,1),Uref(:,2),'r');
% %     Uerr = U-Uref; Uerr = max(myNorm(Uerr))/max(myNorm(Uref));
% %     Terr = stPol-stEx; Terr = max(myNorm(Terr))/max(myNorm(stEx));
% %     Perr = p - pEx; Perr = max(abs(Perr));
% % 
% %     fprintf('& %3.0f & %4.2e & %4.2e & - \\\\ \n',prams.M(1),Uerr,Terr);
% %     %figure;
% %     %subplot(1,2,1); plot(1:np,st(:,1),'-r',1:np,stEx(:,1),'-k');
% %     %subplot(1,2,2); plot(1:np,st(:,2),'-r',1:np,stEx(:,2),'-k');
% %     %pause
% %     %close all
% % end

%% General domain

% % Flow = {'poiseulle','cubic'};
% % fprintf('   Vel       Tracion     Pressure\n');
% % fprintf('--------------------------------------\n');
% % X = [];
% %
% % for dd = 1:size(Flow,2)
% %     flow = Flow{dd};
% %     switch flow
% %         case 'poiseulle'
% %             TRAC = @(x,y) [2*x 1-2*y;1-2*y 2*x];
% %             PRES = @(x,y) -2*x;
% %         case 'cubic'
% %             TRAC = @(x,y) [-6*x*y 3*(x^2+y^2);3*(x^2+y^2) -6*x*y];
% %             PRES = @(x,y) 6*x*y;
% %     end
% %
% %     for a = 5:8
% %         prams.M = 2^a*[1 1 1 1 1 1 1 1 1 1]; np = sum(prams.M);
% %         prams.bd = @(ind,m) sampleBd(ind,m,1,'stokes');
% %
% %         prams.bc = @(x) forcing(x,flow);
% %         domain = fixedBound(prams.M,prams.bd,1);
% %
% %         while(size(X,1)~=128)
% %             xrand = 10*rand(1000,2)-5;
% %             IN = inDomain(domain,xrand(:,1),xrand(:,2),.15);
% %             X = [X;xrand(IN,:)];
% %             X = X(1:min(size(X,1),128),:);
% %         end
% %
% %         XX = []; n = [];
% %         for ind = 1:length(domain)
% %             XX = [XX;domain(ind).X];
% %             n = [n;domain(ind).n];
% %         end
% %
% %         [U trash mu] = farFieldVel(X,[],'confined',prams);
% %         U = reshape(U,[],2);
% %         Uref = forcing(X,flow);
% %
% %         [st p] = evalTraction(prams.M,domain,mu);
% %         %avgS2 = avgStress(prams.M,domain,mu,[],[]);
% %         for ind = 1:size(XX,1)
% %             stEx(ind,:) = TRAC(XX(ind,1),XX(ind,2))*n(ind,:)';
% %             pEx(ind,1)  = PRES(XX(ind,1),XX(ind,2));
% %         end
% %
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         %viewer(domain);
% %         %hold on; quiver(X(:,1),X(:,2),U(:,1),U(:,2));
% %         %hold on; quiver(X(:,1),X(:,2),Uref(:,1),Uref(:,2),'r');
% %         Uerr = U-Uref; Uerr = max(myNorm(Uerr))/max(myNorm(Uref));
% %         Terr = st-stEx; Terr = max(myNorm(Terr))/max(myNorm(stEx));
% %         Perr = p - pEx; Perr = max(abs(Perr))/max(abs(pEx));
% %
% %         fprintf('& %4.2e & %4.2e & %4.2e \\\\ \n',Uerr,Terr,Perr);
% %         %figure;
% %         %subplot(1,2,1); plot(1:np,st(:,1),'-r',1:np,stEx(:,1),'-k');
% %         %subplot(1,2,2); plot(1:np,st(:,2),'-r',1:np,stEx(:,2),'-k');
% %         %pause
% %         %close all
% %     end
% % end
