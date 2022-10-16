function vf = velocity1(tk,pt,time,options,ntrac,nv,...
    posx,posy,fx,fy,den1,den2,const,n)

%disp(tk)
message = ['Time is ' num2str(tk,'%6.3e')];
disp(message)
fid = fopen('time.log','a');
fprintf(fid,'%s\n',message);
fclose(fid);

dt = time(2)-time(1);
ntime = numel(time);
nk = floor((tk-time(1))/dt);
%p = 6;
p = 2;
i = max(1,nk-p/2+1);
%This is not going to be an integer unless p is even.
%What is going on here???
if i+p > ntime
  ntime-p;
end
A = invMatLagrange(5);
C = invMatLagrange(p);
%t = (tk-time(i+1))/((p-1)*dt);
t = (tk-time(i))/((p-1)*dt);
%Value that the Lagrange interpolant will have to be evaluated at

ptx = pt(1:ntrac)';
pty = pt(ntrac+1:2*ntrac)';
%x and y coordinate of the tracers
Vx = [];
Vy = [];
%These will store the velocity at the interpolating time steps


qw = quadratureS(n/2,8,2*pi);
%Quadrature weights and locations 

for k = 0:p-1
  if (strcmp(options.file,'example1_data.bin'))
      vx = pty(end,:);
      vy = zeros(1,ntrac);  
  elseif (strcmp(options.file,'example2_data.bin'))
      vx = zeros(1,ntrac);
      vy = zeros(1,ntrac);  
  elseif (strcmp(options.file,'example3_data.bin'))
      vx = pty(end,:);
      vy = zeros(1,ntrac);  
  elseif (strcmp(options.file,'example4_data.bin'))
      vx = -ptx(end,:);
      vy = pty(end,:);
  elseif (strcmp(options.file,'example6_data.bin'))
      vx = sin(ptx(end,:)).*cos(pty(end,:));
      vy = -cos(ptx(end,:)).*sin(pty(end,:));
  end
  %background velocity at tracer points
    
  vself = [];
  for j=1:nv
    xx = posx(:,j,i+k)';
    yy = posy(:,j,i+k)';
    sa = sqrt(D1FT(xx').^2 + D1FT(yy').^2);
    sa = sa';
    G = kernelS([xx';yy'],qw);
    %SLP due to vesicle j.  Needed for vself
    vselfTemp = G*[fx(:,j,i+k);fy(:,j,i+k)];
    vself = [vself  [vselfTemp(1:n)'; vselfTemp(n+1:2*n)']];
    %Store vself for vesicle j at time i+k
    fx(:,j,i+k) = fx(:,j,i+k)'.*sa*2*pi/n;
    fy(:,j,i+k) = fy(:,j,i+k)'.*sa*2*pi/n;
  end
  h = 1/n*sum(sa)*2*pi/n;

  x = posx(:,:,i+k);
  y = posy(:,:,i+k);
  domain = [x(:)' ; y(:)'];
  ffx = fx(:,:,i+k);
  ffy = fy(:,:,i+k);
  f = [ffx(:)' ; ffy(:)'];
  Xeval = [ptx ; pty];
  %Get positions and charges of all vesicles at time i+k
  col_format = 0;
  fmm = 1;

  [permute istart iend neighbors binof] = ...
      getBins([domain Xeval],n*nv,h);
  [zonev,distv,closestv,indcpv,~,~,theta] = ...
        getZone([domain Xeval],zeros(2,n*nv),n,nv,h,permute,...
        istart,iend,neighbors,binof,A,Xeval);
  %Need zones, bins, etc to do near-singular integration
  %Don't need to pass vself as this is computed on the fly in 
  %eval_velocity_fmm

  [~,~,uux,uuy] = eval_velocity_fmm(domain,f,vself,n,h,...
      col_format,fmm,Xeval,zonev,distv,closestv,indcpv,theta); 
  %evaluate the velocity field due to the vesicles on the tracers
    
  if (options.confined)

    Dbd = [];

    Ri = 10;
    Ro = 20;
    opts.bd = @(ind,m) sampleBd(ind,m,1,'couette','Ri',Ri,'Ro',Ro);
    bound = @(n) fixedBound(n,opts.bd,1);
    %Get information about solid walls

    mu = [den1(:,i+k); den2(:,i+k); const(:,i+k)];
    mu2 = zeros(2*sum(options.M)+3*options.nholes,1);
    %arrange mu in format that evalFieldDirect needs
    istart = 1;
    iend = sum(options.M);
    mu2(1:2:2*iend) = mu(istart:iend);
    mu2(2:2:2*iend) = mu(iend+1:2*iend);
    mu2(2*iend+1:end) = mu(2*iend+1:end);
    nearsing = 1;

    bb = bound(options.M);
    for kk=1:length(options.M)
      v = bb(kk).X;
      Dbd{kk} = kernelD(v(:));
    end
    %Build double-layer potential for all components of the solid walls

    istart1 = 1;
    nd = numel(options.M);
    for k=1:nd
      iend1 = istart1 + options.M(k) - 1;
      sabd = D1FT(bb(k).X(istart1:iend1,1)).^2 + ...
        D1FT(bb(k).X(istart1:iend1,2)).^2;
      arc(k) = sum(sqrt(sabd))*2*pi/options.M(k)/options.M(k);
      %arclength term of solid wall k

      [permute,istart_box,iend_box,neighbors,binof] = ...
          getBins([bb(k).X' Xeval],options.M(k),arc(k));

      [zonebv{k},distbv{k},closestbv{k},indcpbv{k},~,~,thetabv{k}] = ...
          getZone([bb(k).X' Xeval],zeros(2,options.M(k)),options.M(k),1,arc(k),...
          permute,istart_box,iend_box,neighbors,binof,A,Xeval);
      %Get zones, bins, etc for solid walls.  vself is not needed here as
      %it is computed on the fly inside evalFieldDirect
    end


    vself = [];
    istart = 1;
    for kk = 1:numel(options.M)
      iend = istart + options.M(kk) - 1;
      muTemp = [mu(istart:iend);mu((istart:iend)+sum(options.M))];
      vselfTemp = (-1/2*eye(2*options.M(kk)) + Dbd{kk})*muTemp;
      vself = [vself; [vselfTemp(1:end/2);vselfTemp(end/2+1:end)]];
      istart = iend + 1;
    end
    %Calculate vself which is the velocity due to the solid walls
    %evaluated on the solid walls
    Bmu = evalFieldDirect(options.M,mu2,vself,bound,ptx,pty,[],nearsing,...
        zonebv,distbv,closestbv,indcpbv,thetabv);
    %Compute the velocity due to the solid walls

    vx = Bmu(1:ntrac)';
    vy = Bmu(ntrac+1:2*ntrac)';
  end


  vx = vx + uux;
  vy = vy + uuy;
  %Update the velocity due to the vesicles with the velocity due
  %to the solid walls
  %vx = vx;
  %vy = vy;

  Vx = [Vx ; vx];
  Vy = [Vy ; vy];
  %Storing the velocity at all of the interpolation time steps
        
end

for i=1:ntrac
  vx = Vx(:,i);
  vy = Vy(:,i);
  Px = C*vx;
  Py = C*vy;
  vf(:,i) = [polyval(Px,t) polyval(Py,t)]';
end
%Interpolate using 1D Lagrange interpolation

vf = vf(:);
vf = [vf(1:2:end);vf(2:2:end)];
%rearrange velocity in the required manner

%x = pt(1:ntrac);
%y = pt(ntrac+1:2*ntrac);
%u = y.*(1/3-400/3./(x.^2+y.^2));
%v = x.*(-1/3+400/3./(x.^2+y.^2));
%An exact solution for a Couette flow with inner solid wall having
%constant velocity while outer wall fixed.  No vesicles in simulation; solid
%walls are the only thing that drive the flow.





