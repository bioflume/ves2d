clear all;clc;
addpath ../src/

nv = 1;
Ri = 10;            
Ro = 20;            
omegaIn = 1;        
omegaOut = 0;       
reducedArea = .50;   
vesSize = 1;        

%%-- Simulation parameters and options
prams.T = 200;        
prams.n = 128;              
prams.ts= 1e-3;
prams.kappa = 1e-1;                           
prams.flowType = 'confined';             
prams.m = prams.T/prams.ts;                               
prams.order = 2;
prams.M = [256 128];

options.usePlot = 0;                          
options.progressBar = 0;                      
options.useGPU = 0;                           
options.showError = true;                     
options.saveData = true;
options.dataStride = floor(prams.m/100);
options.fileName = './results/A_couette_ra50_2.bin';


prams.bd = @(ind,m) sampleBd(ind,m,1,'couette','Ri',Ri,'Ro',Ro);
prams.bc = @(x) forcing(x,'couette','Ri',Ri,'Ro',Ro,'omegaIn',omegaIn, ...
                        'omegaOut',omegaOut);

prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,...
                                 prams,'direct',options.useGPU);
[t prams options] = monitor(prams,options);
domain = fixedBound(prams.M,prams.bd,1);
X = boundary(prams.n,'couette',domain,'uniform','nv',nv,'scale',vesSize, ...
             'reducedArea',reducedArea,'angle',pi/2,'center',[0;17.5]);


tic
[Xf status] = Ves2D(X,prams,options,@monitor);
toc 
exit;

%%Post process

lineSpec = {'k','c','r','b','m','g'};
RA = [50 60 70 80 90 95];

for idx=1:length(RA)
  ra = RA(idx);
  fileName = ['./results/A_couette_ra' num2str(ra) '.bin'];
  fileId = fopen(fileName,'r');
  Result = fread(fileId,'double');
  fclose(fileId);
  n = prams.n;
  Result = reshape(Result,5*nv*n+1+2*sum(prams.M)+3,[]); 
  
  Xv   = Result(1       :2*nv*n,:);
  Time = Result(5*nv*n+1       ,:);
  mu   = Result(5*nv*n+2:end   ,:);
  
  domain = fixedBound(prams.M,prams.bd,1);
  XX = domain(2).X; l = sqrt(dot(XX,XX,2));
  XX = XX./[l l];
  ds = domain(2).h*domain(2).jacob;
  mu = mu(2*prams.M(1)+1:2*sum(prams.M),:);
  for ii=1:size(mu,2)
    den = reshape(mu(:,ii),2,[])';
    den = den(:,1).*XX(:,2) - den(:,2).*XX(:,1);
    torque(ii) = 2*sum(den.*ds);
  end
  
  volFrac = nv*pi*vesSize^2/domain(1).area;
  effVisc = (1-(Ri/Ro)^2)/(4*pi*Ri*abs(omegaIn-omegaOut))*abs(torque);
  intVisc = (effVisc-1)/volFrac;
  
  Cx = mean(Xv(1  :  n,:)); Cx = reshape(Cx,nv,[]);
  Cy = mean(Xv(n+1:2*n,:)); Cy = reshape(Cy,nv,[]);
  [TH R] = cart2pol(Cx,Cy);
  
  %keyboard
  subplot(2,1,1); hold on;
  plot(Time,R,lineSpec{idx});
  axis([Time(1) Time(end) 12 15]);
  title('Centeroid of the vesicles vs. time');

  subplot(2,1,2); hold on;
  plot(Time,intVisc,lineSpec{idx});
  axis([Time(1) Time(end) 1 3.4]);
  title('Intrinsic viscosity vs. time');
  Leg{idx} = ['Reduced area = ' num2str(ra)];
end

subplot(2,1,1); legend(Leg{:}); grid on;
subplot(2,1,2); legend(Leg{:}); grid on;

%%Giovanni data
% load -ascii results/G_radial_distance_R10_20_O1.dat
% Rg = G_radial_distance_R10_20_O1;
% Tg = Rg(:,1);
% Tg = Tg(5:5:end);
% Rg = Rg(5:5:end,2);
% 
% subplot(1,2,1);
% plot(Time,R,'b',Tg,Rg,'g');
% title('The radial position');
% xlabel('Time'), ylabel('r'); 
% axis([0 11 14.8 15.05]);
% tt = sprintf('|R_A - R_G|/R_G= %2.2e',max(abs(R(1:length(Rg))'-Rg))/max(abs(Rg)));
% text(2,15,tt);
% 
% load -ascii results/G_effective_viscosity_R10_20_O1.dat
% evg = G_effective_viscosity_R10_20_O1;
% evg = evg(5:5:end,2);
% subplot(1,2,2);
% plot(Time,intVisc,'b',Tg,evg,'g');
% title('Intrinsic viscosity');
% xlabel('Time'), ylabel('\eta'); 
% axis([0 11 .2 1.8]);
% tt = sprintf('|\\eta_A - \\eta_G|/\\eta_G= %2.2e',max(abs(intVisc(1:length(evg))'-evg))/max(abs(evg)));
% text(2,1.2,tt);
% legend('Atlanta','Grenoble');

