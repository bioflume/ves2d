
initVes2D            %Initializing Matlab for the simulation

n = 64;              %Number of discretization points on each vesicle                    
nv = 2;              %Number of vesicles
Ri = 10;             %Couette flow inner circle radius
Ro = 20;             %Couette flow outer circle radius
omegaIn = 1;         %Inner circle angular velocity
omegaOut = 0;        %Outer cylinder angular velocity
reducedArea = .8;    %Vesicles' reduced area
volFrac = .01;       %The volume fraction (overrides nv if vesDist is set to volFrac)
vesDist = 'volFrac'; %The distribution of vesicles, 'uniform' or 'random' or 'volFrac'
vesSize = 1;         %Nondimensional size of the vesicle \sqrt(A/pi)
M = [256 128];       %Number of discretization point on the boundary [outer inner]
actionList = {'checkField','simulation','postProcess'};
                     %The actions to be taken = {'checkField','simulation','postProcess'}
                     %'simulation' is the default when actionList is not
                     %present or is empty.
                     
%%-- Simulation parameters and options
prams.T = 10;                                 %Simulation time horizon
prams.ts = .04;                               %Time step
prams.m = prams.T/prams.ts;                   %Number of time steps       
prams.kappa = 1e-1;                           %Bending modulus of the vesicles
prams.order = 1;                              %Time-stepping order

options.usePlot = 1;                          %Whether to plot real-time
options.progressBar = 0;                      %Show progress
options.saveData = 1;                         %Save data as a binary file
options.dataStride = ceil(prams.m/200);       %Frequency of saving data
options.fileName = './couetteFlow.bin';       %Data file name
options.useGPU = 0;                             
options.showError = true;                     %Showing the error in area.

%%-- Setting up the boundary
prams.flowType = 'confined';             
prams.M = M;                     
prams.bd = @(ind,m) sampleBd(ind,m,1,'couette','Ri',Ri,'Ro',Ro);
prams.bc = @(x) forcing(x,'couette','Ri',Ri,'Ro',Ro,'omegaIn',omegaIn,'omegaOut',omegaOut);;
prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,...
                                 prams,'direct',options.useGPU);

[t prams options] = monitor(prams,options);

%%-- Actions
if(~any(strcmp('actionList',who)) || isempty(actionList))
  actionList = {'simulation'};
end

for ii=1:length(actionList)
  
  switch actionList{ii}
    case 'checkField'
     %%-- Checking the velocity
     fprintf('\n\nChecking the velocity field in the absence of vesicles:\n');
     dr = Ro-Ri;
     domain = fixedBound(prams.M,prams.bd,1);
     [theta R] = meshgrid(linspace(0,pi/2,16),linspace(Ri+.1*dr,Ro-.1*dr,16));
     [xx yy] = pol2cart(theta(:),R(:));
 
     ur = prams.bc([xx yy]);
     [u  trash mu] = prams.vInf([xx;yy],[]); u = reshape(u,[],2);
     e = u-ur;
     fprintf('\n   Local inf error is %2.4e.\n',sqrt(max(dot(e,e,2)./dot(ur,ur,2))));
 
     viewer(domain);
     hold on; quiver(xx,yy,u(:,1),u(:,2));
     axis([0 Ro 0 Ro]);
     
    mu = reshape(mu(2*M(1)+1:2*sum(M)),2,[])';
     mu = mu(:,1).*domain(2).X(:,2) - mu(:,2).*domain(2).X(:,1);
     mu = mu./sqrt(dot(domain(2).X,domain(2).X,2));
     torque = 2*sum(domain(2).h*mu.*domain(2).jacob);
     
     effVisc = (1-(Ri/Ro)^2)/(4*pi*Ri*abs(omegaIn-omegaOut))*abs(torque);
     fprintf('   The viscosity of the fluid is %2.2f.\n\n',effVisc);
     disp('Press any key to continue ...'); pause;
     
   case 'simulation'
    %%-- Generating the distribution of vesicles
    domain = fixedBound(prams.M,prams.bd,1);
    X = boundary(n,'couette',domain,vesDist,volFrac,'nv',nv,'scale',vesSize,'reducedArea',reducedArea);
    viewer(X,[],prams,options);
    nv = size(X,2);
    
    %%-- Time simulation
    [XF status] = Ves2D(X,prams,options,@monitor);
    
   case 'postProcess'
    %%-- Reading data 
    if(any(strcmp(who,'status')))
      fileName = status.fileName;
    else
      fileName = input('Enter the name of the data file: ','s');
    end
      
    fileId = fopen(fileName,'r');
    Result = fread(fileId,'double');
    fclose(fileId);
    
    try
      Result = reshape(Result,5*nv*n+1+2*sum(M)+3,[]); 
    catch
      error(['The number of vesicles ''nv'' and number of discretization points ' ...
             '''n''\n are not compatible with the number of entries read ' ...
              'from file.']);
    end
    
    Xv   = Result(1       :2*nv*n,:);
    Time = Result(5*nv*n+1       ,:);
    mu   = Result(5*nv*n+2:end   ,:);
    
    domain = fixedBound(prams.M,prams.bd,1);
    XX = domain(2).X; l = sqrt(dot(XX,XX,2));
    XX = XX./[l l];
    ds = domain(2).h*domain(2).jacob;
    mu = mu(2*M(1)+1:2*sum(M),:);
    for ii=1:size(mu,2)
      den = reshape(mu(:,ii),2,[])';
      den = den(:,1).*XX(:,2) - den(:,2).*XX(:,1);
      torque(ii) = 2*sum(den.*ds);
    end

    volFrac = nv*pi*vesSize^2/domain(1).area;
    effVisc = (1-(Ri/Ro)^2)/(4*pi*Ri*abs(omegaIn-omegaOut))*abs(torque);
    intVisc = (effVisc-1)/volFrac;
    
    Xv = reshape(Xv,2*n,[]);
    subplot(3,1,1); 
    viewer(Xv,[],prams,options);
    title('Vesicle configuration');

    subplot(3,1,2); 
    %plot(Time,effVisc); axis([Time(1) Time(end) 1 1+1.2*(max(effVisc)-1)]);
    plot(Time,intVisc); axis([Time(1) Time(end) 1 1+1.2*(max(intVisc)-1)]);
    title('Effective viscosity vs. time');
    
    Cx = mean(Xv(1  :  n,:)); Cx = reshape(Cx,nv,[]);
    Cy = mean(Xv(n+1:2*n,:)); Cy = reshape(Cy,nv,[]);
    [TH R] = cart2pol(Cx,Cy);
    
    subplot(3,1,3);
    plot(Time,R);
    axis([Time(1) Time(end) Ri Ro]);
    title('Centeroid of the vesicles vs. time');
  end
end
