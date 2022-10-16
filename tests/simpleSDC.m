clf
addpath ../src

%benOp = 0;
<<<<<<< HEAD
benOp = 2;
tbenOp = 4;
=======
%benOp = 2;
benOp = 4;
>>>>>>> 4251dbc133f68e0e32a4cd43e83130e8a5a18c1e
% number of arclength derivatives in bending operator


T = 100;
<<<<<<< HEAD
M = 10*2.^(0:8);
nsdc = 20;
=======
%M = 10*2.^(0:8);
M = 10*2.^(9:12);
nsdc = 5;
>>>>>>> 4251dbc133f68e0e32a4cd43e83130e8a5a18c1e

oc = curve;
N = 192;
theta = (0:N-1)'*2*pi/N;
%theta = oc.redistributeArcLength(cos(theta),3*sin(theta));
D1 = real(fft1.fourierDiff(N));


Xprov = zeros(2*N,1,3);
eX = zeros(2*N,1,3);
vel = zeros(2*N,1,3);
residual = zeros(2*N,1,3);
A = zeros(2*N,2*N,3);

time = [0]; 
ea = [0];
op = poten(N);
fprintf('***************************************\n')
fprintf('Number of points on interface is %i\n',N);
fprintf('Time horizon is                  %i\n',T);
fprintf('Number of SDC corrections is     %i\n',nsdc);
fprintf('Order of spatial derivative is   %i\n',benOp);
fprintf('***************************************\n')
for kk = 1:numel(M);
  X = [1*cos(theta);3*sin(theta)];
  [~,a0,~] = oc.geomProp(X);
  m = M(kk);
  dt = T/m;
  for ntime = 1:m
    Xprov(:,:,1) = X;

    % compute provisional solution
    for j = 1:2
      vesicle = capsules(Xprov(:,:,j),[],[],1,1,0);
      G = op.stokesSLmatrix(vesicle); 
      if benOp == 0
        Ben = -eye(2*N);
      elseif benOp == 2
        arcDeriv = vesicle.isa(:,ones(N,1)).*D1;
        Ben = kron(eye(2),arcDeriv^2);
      elseif benOp == 4
        arcDeriv = vesicle.isa(:,ones(N,1)).*D1;
        Ben = -kron(eye(2),arcDeriv^4);
      end
      A(:,:,j) = G*Ben;
      Xprov(:,:,j+1) = (eye(2*N) - dt/2*A(:,:,j))\Xprov(:,:,j);
    end

    % build vesicle at final GL point
    vesicle = capsules(Xprov(:,:,3),[],[],1,1,0);
    G = op.stokesSLmatrix(vesicle); 
    if benOp == 0
      Ben = -eye(2*N);
    elseif benOp == 2
      arcDeriv = vesicle.isa(:,ones(N,1)).*D1;
      Ben = kron(eye(2),arcDeriv^2);
    elseif benOp == 4
      arcDeriv = vesicle.isa(:,ones(N,1)).*D1;
      Ben = -kron(eye(2),arcDeriv^4);
    end
    A(:,:,3) = G*Ben;

    % compute velocity
    for j = 1:3
      vel(:,:,j) = A(:,:,j)*Xprov(:,:,j);
    end
    % integrate velocity
    velInt = lobattoInt(vel);

    % compute residual
    for n = 1:3
      residual(:,:,n) = Xprov(:,:,1) - Xprov(:,:,n) + ...
        dt/2 * velInt(:,:,n);
    end


    % Do SDC corrections
    for k = 1:nsdc
      for j = 2:3
        eX(:,:,j) = (eye(2*N) - dt/2*A(:,:,j))\...
            (residual(:,:,j) - residual(:,:,j-1) + eX(:,:,j-1));

      end
      Xprov = Xprov + eX;

      % update operator G*Bending
      for j = 1:3
        vesicle = capsules(Xprov(:,:,j),[],[],1,1,0);
        G = op.stokesSLmatrix(vesicle); 
        if benOp == 0
          Ben = -eye(2*N);
        elseif benOp == 2
          arcDeriv = vesicle.isa(:,ones(N,1)).*D1;
          Ben = kron(eye(2),arcDeriv^2);
        elseif benOp == 4
          arcDeriv = vesicle.isa(:,ones(N,1)).*D1;
          Ben = -kron(eye(2),arcDeriv^4);
        end
        A(:,:,j) = G*Ben;
      end

      % compute velocity
      for j = 1:3
        vel(:,:,j) = A(:,:,j)*Xprov(:,:,j);
      end
      % integrate velocity
      velInt = lobattoInt(vel);

      % compute residual
      for j = 1:3
        residual(:,:,j) = Xprov(:,:,1) - Xprov(:,:,j) + ...
          dt/2 * velInt(:,:,j);
      end
    end

    X = Xprov(:,:,end);
    [~,a,~] = oc.geomProp(X);
    ea(ntime+1) = abs(a - a0)/abs(a0);
%    time(ntime+1) = time(ntime) + dt;
%    plot(X(1:N),X(N+1:end),'r')
%    axis equal
%    axis([-2 2 -3 3])
%    title(num2str(time(end),'%4.2e'));
%    pause(0.01)
%    pause
  end
  fprintf('Number of time step is %i\n',m);
  fprintf('Error in area is       %4.2e\n\n',ea(end));
end











