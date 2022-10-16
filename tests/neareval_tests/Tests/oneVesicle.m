%August 8, 2012
%Walid's code to do convergence study for 2dnear paper
%Cleaned up and udapted a whole bunch

clear all 
addpath ../../src
%This has vesicle-vesicle interactions
addpath ../../src/neareval
%This has getBins and getZone
addpath Fordriver
%This has some small routines of Walids
ifmm = 0;
%Whether we use the fmm or not
%example = 'density';
%example = 'extensional';
example = 'stokeslet';


% Convergence study where the exact solution is computed using a big number of points
% One vesicles are considered

N1 = 5;
A = invMatLagrange(N1);

%  exact solution 
%Nexact = 2^11;
Nves = 1;
Nexact = 2^10;

% Domain
shape = 1;
r1 = 1;
r2 = 1.5;
%r2 = 1;
%R = 1;
%a = 4;
%c3 = .2;
shapeparam = [r1 r2];
[x y deriv] = getDomain(Nexact,shape,shapeparam);
domain = [x ; y];
h = 1/Nexact*(sum(deriv)*2*pi/Nexact);
%arclength term


% Vself
qw = quadratureS(Nexact/2,8,2*pi);
G = kernelS([x' ; y'],qw);
%Need single-layer potential to compute vself
%and to solve for the density given a boundary condition

% Targets % 
if (strcmp(example,'stokeslet'))
  %dr = 0.0005;
  dr = 0.01;
else
  dr = -0.01;
end
%%dr = .0005;
%dr = -0.01;
[targetsx,targetsy,~] = getDomain(32,shape,shapeparam);
targetsx = (1+dr)*targetsx;
targetsy = (1+dr)*targetsy;
targets = [targetsx ; targetsy];

if strcmp(example,'density')
  [f1 f2] = getdensity(Nexact);
  %density is chosen a priori by user
  vself = G*[f1' ; f2'];
  vself = [vself(1:Nexact)' ; vself(Nexact+1:2*Nexact)'];

  f1 = f1.*deriv;
  f2 = f2.*deriv;
  f = [f1 ; f2];
  %Multiply density function by arclength

  % Near singular evaluation  
  [permute istart_box iend_box neighbors binof] = ...
    getBins([domain targets],Nexact*Nves,h);

  [zone,dist,closest,indcp,~,~,theta] = getZone(...
      [domain targets],zeros(2,Nexact*Nves),Nexact,1,h,...
    permute,istart_box,iend_box,neighbors,binof,A,targets);

  if(any(zone(Nexact*Nves+1:end)));
    disp('Using near-singular for near region of exact solution')
  end

  [uu,vv,uux,uuy] = eval_velocity_fmm(domain,f*2*pi/Nexact,...
    vself,Nexact,h,0,ifmm,targets,zone,dist,closest,indcp,theta);

  uexact = uux;
  vexact = uuy;
  %'Exact' solution using an over-refined grid
end


if strcmp(example,'extensional');
  uexact = -targetsx;
  vexact = targetsy;
end

if strcmp(example,'stokeslet');
  rho2 = targetsx.^2+targetsy.^2;
  uexact = -log(rho2)/2+(targetsx+targetsy)./rho2 .* targetsx;
  vexact = -log(rho2)/2+(targetsx+targetsy)./rho2 .* targetsy;
end

errx = [];
erry = [];
%Nset = [16 32 64 128 256 512];
Nset = [16 32 64 128];

for N = Nset
  disp(N)
  [x y deriv] = getDomain(N,shape,shapeparam);
  
  domain = [x ; y];
  qw = quadratureS(N/2,8,2*pi);
  G = kernelS([x' ; y'],qw);
  h = 1/N*(sum(deriv)*2*pi/N);
  %arclength term

  if (strcmp(example,'density'))
    [f1 f2] = getdensity(N);
    %density is chosen a priori by user
    vself = G*[f1' ; f2'];
    vself = [vself(1:N)' ; vself(N+1:2*N)'];
  end

  if (strcmp(example,'extensional'))
    bcx = -x; bcy = y;
    f = G\[bcx';bcy'];
    f1 = f(1:N)';
    f2 = f(N+1:end)';
    
    vself = [bcx;bcy];
  end

  if (strcmp(example,'stokeslet'))
    rho2 = x.^2 + y.^2;
    bcx = -log(rho2)/2+(x+y)./rho2 .* x;
    bcy = -log(rho2)/2+(x+y)./rho2 .* y;
    f = G\[bcx';bcy'];
    f1 = f(1:N)';
    f2 = f(N+1:end)';
    
    vself = [bcx;bcy];
  end

  f1 = f1.*deriv;
  f2 = f2.*deriv;
  f = [f1 ; f2];
    
  [permute istart_box iend_box neighbors binof] = ...
    getBins([domain targets],N*Nves,h);

  [zone,dist,closest,indcp,~,~,theta] = getZone([domain targets],...
    zeros(2,N*Nves),N,1,h,...
    permute,istart_box,iend_box,neighbors,binof,A,targets);

  if(any(zone(Nexact*Nves+1:end)-1));
    disp('Not using near-singular for small values of N')
  end

  [uu,vv,uux,uuy] = eval_velocity_fmm(domain,f*2*pi/N,...
    vself,N,h,0,ifmm,targets,zone,dist,closest,indcp,theta);

  u = uux;
  v = uuy;
  errx = [errx norm(u - uexact,inf)];
  erry = [erry norm(v - vexact,inf)];
end

err = max(errx,erry);
%Take the maximum of the error in the x and y direction

%loglog(Nset,err)
%grid
%title('Convergence of the nearly singular integration in log-log scale')
%xlabel('N')
%ylabel('Relative error') 
%hold on
%I=floor(log2(Nset));
%I = find(Nset==2.^I);
%loglog(Nset(I),err(I),'.')
%loglog(Nset(I),err(I),'o')
%set(gca,'xtick',[32 64 128 256 512])
%set(gca,'xlim',[26 512+16])
%set(gca,'fontsize',16)
%set(get(gca,'xlabel'),'fontsize',16)
%set(get(gca,'ylabel'),'fontsize',16)
%set(get(gca,'title'),'fontsize',16)


% Nset1 = Nset(end-10:end);
% err1 = err(end-10:end);
% P=polyfit(log(Nset1),log(err1),1);
% display(['The order of convergence is ' num2str(-P(1))]);






