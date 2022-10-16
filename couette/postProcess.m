function [Xv Time C effVisc options prams] = postProcess(fileSpec, varargin)
% arguments : 'TimePlot', 'Trajectory', 'ErrorPlot', 'EffectiveViscosity'

  fname = Stringify(fileSpec);

  %- Return values
  Xv = []; Time = []; C = []; effVisc = []; options = []; params = [];
 
  %- Reading from file
  if ( ~exist([fname '.mat']) )
    return;
  else 
    load([fname '.mat']);
  end

  if ( exist(options.fileName) )
    fileId = fopen(options.fileName,'r');
  elseif ( exist (['./results/' options.fileName]) )
    fileId = fopen(['./results/' options.fileName],'r');
  else
    error(['The file "' options.fileName '" does not exist.']);
  end    

  Result = fread(fileId,'double');
  fclose(fileId);

  if (~isempty(prams.n) ), n = prams.n;else, prams.n = n; end
  if (~isempty(prams.nv) ), nv = prams.nv;else, prams.nv = nv; end

  
  Result = reshape(Result,5*nv*n+1+2*sum(prams.M)+3,[]); 
  Xv   = Result(1       :2*nv*n,:);
  Time = Result(5*nv*n+1       ,:);
  mu   = Result(5*nv*n+2:end   ,:);  

  %- Fixing the time
  Dt = diff([0 Time]);
  ind = find(Dt < 0);
  Dt(ind) = Time(ind + 1);
  Time = cumsum(Dt);
  
  %- Analysis
  nstep = size(Xv,2);
  
  Energy = zeros(nstep,1);
  C      = zeros(nstep, 2*nv);
  Area   = zeros(nstep,nv);
  Length = zeros(nstep,nv);
  
  effVisc = zeros(nstep,1);

  for ii=1:size(Xv,2)
    Xs = reshape(Xv(:,ii),[],nv);
    [trash Area(ii,:) Length(ii,:) c] = reducedVolume(Xs);
    C(ii,:) = [c(1,:) c(2,:)];
    
    xHat = fft(Xs(1:n,:))/2/pi; 
    yHat = fft(Xs(n+1:2*n,:))/2/pi;
    xyHat = [xHat(:) yHat(:)]; 
    Energy(ii) = sum(dot(xyHat,xyHat,2));
  end

  if ( any(strcmp(varargin,'TimePlot')) )

    ii = find(strcmp(varargin,'TimePlot'));
    if ( length(varargin)>ii && any(ishandle(varargin{ii+1})) )
      axes(varargin{ii+1});
    else
      figure;
    end
    nameformat = ceil(log10(nstep));
    nameformat = ['%0' num2str(nameformat) 'u'];
    for ii=1:size(Xv,2)
      Xs = reshape(Xv(:,ii),[],nv);
      cla;
      viewer(Xs,[],prams,options);
      title(num2str(Time(ii),'%02.2f'));
      drawnow;
      
      jj = find(strcmp(varargin,'SaveFig'));
      if(~isempty(jj))
        saveas(gcf,[varargin{jj+1} num2str(ii,nameformat)],'jpg');
      end
      drawnow;
      pause(.04);
    end
  end
  
  if ( any(strcmp(varargin,'Trajectory')) )

    ii = find(strcmp(varargin,'Trajectory'));
    if ( length(varargin)>ii && any(ishandle(varargin{ii+1})) )
      axes(varargin{ii+1});
    else
      figure;
    end
    
    [T R] = cart2pol(C(:,1:nv),C(:,nv+1:end));
    plot(gca,Time, R);
    drawnow;
  end

  if ( any(strcmp(varargin,'ErrorPlot')) )

    ii = find(strcmp(varargin,'ErrorPlot'));
    if ( length(varargin)>ii && any(ishandle(varargin{ii+1})) )
      axes(varargin{ii+1});
    else;
      figure;
    end
    
    subplot(1,3,1); plot(Energy/Energy(1));
    title('l^2 energy');

    Area   = abs(1-Area./repmat(Area(1,:),nstep,1));
    Length = abs(1-Length./repmat(Length(1,:),nstep,1));
    subplot(1,3,2); semilogy(max(Area,[],2)); title('Error in Area');
    subplot(1,3,3); semilogy(max(Length,[],2)); title('Error in length');
    drawnow;
  end

  if ( any(strcmp(varargin,'EffectiveViscosity')) )

    M = prams.M;
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

    %volFrac = nv*pi*vesSize^2/domain(1).area;
    effVisc = (1-(Ri/Ro)^2)/(4*pi*Ri*abs(omegaIn-omegaOut))*abs(torque);
    
    ii = find(strcmp(varargin,'EffectiveViscosity')); 
    if ( length(varargin)>ii && any(ishandle(varargin{ii+1})) )
      axes(varargin{ii+1});
    else
      figure;
    end
    
    plot(Time,effVisc); axis([Time(1) Time(end) 1 1+1.2*(max(effVisc)-1)]);
    title('Effective viscosity vs. time');
    drawnow;
  end