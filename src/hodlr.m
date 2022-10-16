classdef hodlr

properties

    lev_max;
    lmax;
    etol;
    memsize;
    profile;
    
    L;
    SL;
    LS;
    invS1;
    invS2;
    Q;
    R;
    D;
    l;
    partition;
    W;
    Invz;
    Invd;
    QR;
    H;

    newIDs
    oldIDs
    
    memmapL;
    memmapQ;
    memmapR;
    memmapD;
    memmapW;
    memmapInvz;
    memmapInvd;
    memmapM21;
    memmapM22;
    idWall2gids;
   
    isAdaptive;
    isOutofCore;
    logFile;
    verbose;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = hodlr(lev_max,etol,memsize,isOutofCore,isAdaptive,logFile)
   o.lev_max = lev_max;
   o.etol = etol;
   o.l = 1000*ones(o.lev_max,1);
   o.memsize = memsize;
   o.isOutofCore = isOutofCore;
   o.isAdaptive = isAdaptive;
   o.logFile = logFile;
   o.profile = false;
   o.verbose = false;
   fid = fopen(o.logFile,'w');
   fprintf(fid,'***********************************************\n');
   fprintf(fid,'HODLR is called.\n');
   fprintf(fid,'Maximum level is %d\n',lev_max);
   fprintf(fid,'Error tolerance is %2.3e\n',etol);
   fclose(fid);
end %hodlr: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [o]=hodlr_no_FMM(o,matFiles)     % matrix, low-rank, partitioning, no. of levels

    if o.profile
        thodlr0 = tic; 
    end
    
    logFileID = fopen(o.logFile,'a');  
    
    if o.isOutofCore
        n = size(o.memmapL.Data.L,1);
    else
        n = size(o.L,1);
    end
    ltest = 10;
    lev = o.lev_max+1;
    p_ind(1) = 0;  p_ind(2^(lev-1)+1) = n;
    p_ind(2:2^(lev-1)) = sort(o.partition(1:2^(lev-1)-1));
    m = max(p_ind(2:end)-p_ind(1:end-1));

    if(o.isOutofCore)
        filename = [matFiles 'O3.bin'];
        fileID = fopen(filename,'w');
        fclose(fileID);
    else
        o.Q = cell(2^(o.lev_max+1)-2,1);
        o.R = cell(2^(o.lev_max+1)-2,1);
        o.D = cell(2^o.lev_max,1);
        O3 = zeros(n,m);
    end
    
    Qtmp = zeros(n,o.lmax);
    O1 = zeros(n,o.lmax);
    O2 = zeros(n,o.lmax); 

    normA = zeros(2^o.lev_max,1);
    p_ind = zeros(2^o.lev_max+1,1);

    for lev = 1:o.lev_max
        fprintf(logFileID,'-----------------------------\n');
        fprintf(logFileID,'At level %d of HODLR\n',lev);
        if o.verbose
        fprintf('At level %d of HODLR\n',lev);
        end
        
        p_ind(1) = 0;  p_ind(2^(lev)+1) = n;
        p_ind(2:2^(lev))=sort(o.partition(1:2^(lev)-1));
        minblocksize=min(p_ind(2:2^(lev)+1)-p_ind(1:2^(lev)));
        lguess=min(o.l(lev)/2,minblocksize);
        lguessnew=0;
        eguess=realmax;
        p_ind(1)=0;  p_ind(2^(lev-1)+1)=n;
        p_ind(2:2^(lev-1))=sort(o.partition(1:2^(lev-1)-1));
        while(eguess>o.etol)  
            if lguess==lguessnew
                if(lguess>minblocksize)
                    fprintf(logFileID,'Cannot meet accuracy requirement at level %d\n',lev);
                    if o.verbose
                    fprintf('Cannot meet accuracy requirement at level %d\n',lev);
                    end
                    break;
                end;
                lguess= 2*lguess;
                Qtmp = zeros(n,lguess);
                O1 = zeros(n,lguess);
                O2 = zeros(n,lguess);
            end;
            lguess=min(2*lguess,minblocksize);
            if o.verbose
            fprintf('lguess = %d\n',lguess);
            end
            fprintf(logFileID,'-----------------------------\n');
            fprintf(logFileID,'lguess = %d\n',lguess);
            
            if o.profile
                tAtimesG0 = tic; 
            end
            
            [Y1,Y2]=o.AtimesG(O1,O2,p_ind,lguess,lev,n);

            if o.profile
                t = toc(tAtimesG0);
                fprintf('Multiplying matrix with random gaussians, size %d, takes %2.2e seconds\n',lguess,t);
                fprintf(logFileID,'Multiplying matrix with random gaussians, size %d, takes %2.2e seconds\n',lguess,t);
            end

            if o.profile
                tQR0 = tic; 
            end
            
            for i=1:(2^(lev-1))
                ll=p_ind(i)+1;  rr=p_ind(i+1);  xx=o.partition(2^(lev-1)+i-1);
                [Qtmp(ll:xx,1:lguess),Rtmp]=qr(Y1(ll:xx,1:lguess),0);
                sig=sort(abs(diag(Rtmp)),'descend');
                normA(2*i-1)=sig(1);
                [Qtmp(xx+1:rr,1:lguess),Rtmp]=qr(Y2(xx+1:rr,1:lguess),0);
                sig=sort(abs(diag(Rtmp)),'descend');
                normA(2*i)=sig(1);
            end;
            clear Rtmp;
            if o.profile
                t = toc(tQR0);
                fprintf('QR factorization to get orthonormal basis for column space of the matrix takes %2.2e seconds\n',t);
                fprintf(logFileID,'QR factorization to get orthonormal basis for column space of the matrix takes %2.2e seconds\n',t);
            end
            
            if o.profile
                tAtimesG1 = tic; 
            end
            
            [Y1,Y2]=o.AtimesG(O1,O2,p_ind,ltest,lev,n);

            if o.profile
                t = toc(tAtimesG1);
                fprintf('Multiplying matrix with random gaussians, size %d, takes %2.2e seconds\n',ltest,t);
                fprintf(logFileID,'Multiplying matrix with random gaussians, size %d, takes %2.2e seconds\n',ltest,t);
            end

            if o.profile
                teguess = tic; 
            end
            
            eguess=o.check_error(Qtmp,Y1,Y2,n,lev,lguess,ltest,normA);

            if o.profile
                t = toc(teguess);
                fprintf('Error in approximating column space = %2.2e for low-rank = %d takes %2.2e seconds\n',eguess,lguess,t);
                fprintf(logFileID,'Error in approximating column space = %2.2e for low-rank = %d takes %2.2e seconds\n',eguess,lguess,t);
            end
            
            lguessnew=min(2*lguess,minblocksize);
        end;
        lleft=0;  lright=lguess;
        while(lright-lleft>10)
            lmid=floor((lleft+lright)/2);
            
            if o.profile
                teguess = tic; 
            end
            
            e = o.check_error(Qtmp,Y1,Y2,n,lev,lmid,ltest,normA);

            if o.profile
                t = toc(teguess);
                fprintf('Error in approximating column space = %2.2e for low-rank = %d takes %2.2e seconds\n',e,lmid,t);
                fprintf(logFileID,'-----------------------------\n');
                fprintf(logFileID,'Error in approximating column space = %2.2e for low-rank = %d takes %2.2e seconds\n',e,lmid,t);
            end
            
            if(e>o.etol)
                lleft=lmid;
            else
                lright=lmid;
            end;
        end;
        o.l(lev)=lright;
        jstart=1+sum(o.l(1:lev-1));
        jend=jstart+o.l(lev)-1;

        for i=1:(2^(lev-1))
            ll=p_ind(i)+1;  rr=p_ind(i+1);  xx=o.partition(2^(lev-1)+i-1);
            O1(ll:xx,1:o.l(lev)) = Qtmp(ll:xx,1:o.l(lev));
            O1(xx+1:rr,1:o.l(lev)) = zeros(rr-xx,o.l(lev));
            O2(ll:xx,1:o.l(lev)) = zeros(xx-ll+1,o.l(lev));
            O2(xx+1:rr,1:o.l(lev)) = Qtmp(xx+1:rr,1:o.l(lev));
        end;
            
        if o.profile
            tAtransposeG0 = tic; 
        end
            
        [Y1,Y2] = o.AtransposetimesG(O1,O2,p_ind,o.l(lev),lev,n);
        
        if o.profile
            t = toc(tAtransposeG0);
            fprintf('Multiplying transpose of matrix with orthogonal basis, size %d, takes %2.2e seconds\n',o.l(lev),t);
            fprintf(logFileID,'Multiplying transpose of matrix with orthogonal basis, size %d, takes %2.2e seconds\n',o.l(lev),t);
        end

        if(o.isOutofCore)
            Rtmp = zeros(o.l(lev),n);
            filename = [matFiles 'Q.bin'];
            fileIDQ = fopen(filename,'a');
            filename = [matFiles 'R.bin'];
            fileIDR = fopen(filename,'a');
            
            if o.profile
                tQR0 = tic; 
            end
            
            for i=1:(2^(lev-1))
                ll=p_ind(i)+1;  rr=p_ind(i+1);  xx=o.partition(2^(lev-1)+i-1);
                [Qt1,Rtmp(1:o.l(lev),ll:rr-xx+ll-1)]=qr(Y2(xx+1:rr,1:o.l(lev))',0);
                [Qt2,Rtmp(1:o.l(lev),rr-xx+ll:rr)]=qr(Y1(ll:xx,1:o.l(lev))',0);
                Qtmp(ll:xx,1:o.l(lev))=Qtmp(ll:xx,1:o.l(lev))*Qt1;
                Qtmp(xx+1:rr,1:o.l(lev))=Qtmp(xx+1:rr,1:o.l(lev))*Qt2;
            end;
            

            if o.profile
                t = toc(tQR0);
                fprintf('2nd QR factorization of the matrix takes %2.2e seconds\n',t);
                fprintf(logFileID,'2nd QR factorization of the matrix takes %2.2e seconds\n',t);
            end
            
            fwrite(fileIDQ,Qtmp(:,1:o.l(lev)),'double');
            fwrite(fileIDR,Rtmp','double');
            fclose(fileIDQ);
            fclose(fileIDR);
            filename = [matFiles 'Q.bin'];        
            o.memmapQ=memmapfile(filename, ...
              'Format', {'double' [n sum(o.l(1:lev))] 'M'},     ...
              'Writable', true);
            filename = [matFiles 'R.bin'];
            o.memmapR=memmapfile(filename, ...
              'Format', {'double' [n sum(o.l(1:lev))] 'M'},     ...
              'Writable', true);
            clear Qt1 Qt2 Y1 Y2 Rtmp;
        else
            for i=1:(2^(lev-1))
                ll=p_ind(i)+1;  rr=p_ind(i+1);  xx=o.partition(2^(lev-1)+i-1);
                [Qt1,o.R{2^lev - 2 + 2 * i - 1}]=qr(Y2(xx+1:rr,1:o.l(lev))',0);
                [Qt2,o.R{2^lev - 2 + 2 * i - 0}]=qr(Y1(ll:xx,1:o.l(lev))',0);
                o.Q{2^lev - 2 + 2 * i - 1}=Qtmp(ll:xx,1:o.l(lev))*Qt1;
                o.Q{2^lev - 2 + 2 * i - 0}=Qtmp(xx+1:rr,1:o.l(lev))*Qt2;
            end;            
            clear Qt1 Qt2;
        end
    end;
    clear O1 O2 Qtmp Rtmp Y1 Y2;
    
    lev=o.lev_max+1;
    p_ind(1)=0;  p_ind(2^(lev-1)+1)=n;
    p_ind(2:2^(lev-1))=sort(o.partition(1:2^(lev-1)-1));
    
    if(o.isOutofCore)
        maxbsize = floor(1e9*o.memsize/(8*n));
        bsize = min(m,maxbsize);
        
        if o.profile
            tD0 = tic; 
        end

        filename = [matFiles 'O3.bin'];
        fileID = fopen(filename,'a');
        for j = 1:floor(m/bsize)
            O3tmp = zeros(n,bsize);
            jstart = 1+(j-1)*bsize;
            jend = j*bsize;     
            for i=1:(2^(lev-1))
                ll = p_ind(i)+1;  rr = p_ind(i+1);
                if(rr-ll+1<jend)
                    O3tmp(ll+jstart-1:rr,1:rr-ll-jstart+2) = eye(rr-ll-jstart+2,rr-ll-jstart+2);
                else
                    O3tmp(ll+jstart-1:rr,:) = eye(rr-ll-jstart+2,jend-jstart+1);
                end;
            end;
            fwrite(fileID,O3tmp,'double');
        end;
        jstart = 1+j*bsize;
        O3tmp = zeros(n,m-bsize*floor(m/bsize));
        for i=1:(2^(lev-1))
            ll = p_ind(i)+1;  rr = p_ind(i+1);
            if(rr-ll+1<m)
                O3tmp(ll+jstart-1:rr,1:rr-ll-jstart+2) = eye(rr-ll-jstart+2,rr-ll-jstart+2);
            else
                O3tmp(ll+jstart-1:rr,:) = eye(rr-ll-jstart+2,m-jstart+1);
            end;
        end;        
        fwrite(fileID,O3tmp,'double');
        fclose(fileID);

        if o.profile
            t = toc(tD0);
            fprintf('Creating indentity matrices takes %2.2e seconds\n',t);
            fprintf(logFileID,'-----------------------------\n');
            fprintf(logFileID,'Creating indentity matrices takes %2.2e seconds\n',t);
        end
            
        memmapO3=memmapfile(filename, ...
          'Format', {'double' [n m] 'M'},     ...
          'Writable', true);

        if o.profile
            tD0 = tic; 
        end
        
        filename = [matFiles 'D.bin'];
        fileID = fopen(filename,'a');
        Dtmp = zeros(n,bsize);
        for j = 1:floor(m/bsize)
            jstart = 1+(j-1)*bsize;
            jend = j*bsize;            
            for i = 1:floor(n/bsize)
                istart = 1+(i-1)*bsize;
                iend = i*bsize;
                Dtmp(istart:iend,:) = o.memmapL.Data.L(istart:iend,:)*memmapO3.Data.M(:,jstart:jend);
            end;
            istart = 1+i*bsize;
            Dtmp(istart:end,:) = o.memmapL.Data.L(istart:end,:)*memmapO3.Data.M(:,jstart:jend);
            fwrite(fileID,Dtmp,'double');
        end;
        jstart = 1+j*bsize;
        Dtmp = zeros(n,m-bsize*floor(m/bsize));
        for i = 1:floor(n/bsize)
            istart = 1+(i-1)*bsize;
            iend = i*bsize;
            Dtmp(istart:iend,:) = o.memmapL.Data.L(istart:iend,:)*memmapO3.Data.M(:,jstart:end);
        end;
        istart = 1+i*bsize;
        Dtmp(istart:end,:) = o.memmapL.Data.L(istart:end,:)*memmapO3.Data.M(:,jstart:end); 
        fwrite(fileID,Dtmp,'double');
        fclose(fileID);  
        
        o.memmapD=memmapfile(filename, ...
          'Format', {'double' [n m] 'M'},     ...
          'Writable', true);
      
        for i = 1:floor(m/bsize)
            istart = 1+(i-1)*bsize;
            iend = i*bsize;
            o.memmapD.Data.M(:,istart:iend) = o.memmapD.Data.M(:,istart:iend) - o.mult(memmapO3.Data.M(:,istart:iend),lev-1,n,0);
        end;
        istart = 1+i*bsize;
        o.memmapD.Data.M(:,istart:end) = o.memmapD.Data.M(:,istart:end) - o.mult(memmapO3.Data.M(:,istart:end),lev-1,n,0);

        if o.profile
            t = toc(tD0);
            fprintf('Extracting diagonal matrices takes %2.2e seconds\n',t);
            fprintf(logFileID,'-----------------------------\n');
            fprintf(logFileID,'Extracting diagonal matrices takes %2.2e seconds\n',t);
        end
        
    else
        for i=1:(2^(lev-1))
            ll = p_ind(i)+1;  rr = p_ind(i+1);
            O3(ll:rr,1:rr-ll+1) = eye(rr-ll+1,rr-ll+1);
        end;
        
        Dtmp = o.L*O3;
        Dtmp = Dtmp-o.mult(O3,lev-1,n,0);
        for i=1:(2^(lev-1))
            ll = p_ind(i)+1;  rr = p_ind(i+1);
            o.D{i} = Dtmp(ll:rr,1:rr-ll+1);
        end;
    end;
    
    clear O3 Dtmp;
    
    o.lmax = max(o.l);
    
    if o.verbose
    fprintf('The maximum block size is chosen as %2d\n',o.lmax);
    end
    fprintf(logFileID,'-----------------------------\n');
    fprintf(logFileID,'The maximum block size is chosen as %2d\n',o.lmax);

    if o.profile
        t = toc(thodlr0);
        fprintf('Hodlr takes %2.2e seconds\n',t);
        fprintf(logFileID,'Hodlr takes %2.2e seconds\n',t);
    end
    fclose(logFileID);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o=factorize(o,n,matFiles)
    logFileID = fopen(o.logFile,'a');
    
    if o.profile
        tfact0 = tic; 
    end
    
    lsum=sum(o.l);
	cols_lev=zeros(o.lev_max,1);
    for i=1:o.lev_max
        cols_lev(i)=o.l(i)*2^(i-1);
    end
    
	Wp=zeros(n,2*o.lmax);
    Vp=zeros(2*o.lmax,n);
    
    p_ind=zeros(2^o.lev_max+1,1);
    lev=o.lev_max+1;
    p_ind(1)=0;  p_ind(2^(lev-1)+1)=n;
    p_ind(2:2^(lev-1))=sort(o.partition(1:2^(lev-1)-1));
    m=max(p_ind(2:end)-p_ind(1:end-1));
    if(o.isOutofCore)
        filename = [matFiles 'Invd.bin'];
        fileIDInvd = fopen(filename,'a');
    else
        o.W=cell(2^(o.lev_max + 1) - 2 , 1);
        o.Invz=cell(2^o.lev_max - 1 , 1);
        o.Invd=cell(2^o.lev_max , 1);
    end;

    

    if o.profile
        tLU0 = tic; 
    end
        
    for i=1:2^o.lev_max
        ll=p_ind(i)+1;  rr=p_ind(i+1);
        if(o.isOutofCore)
            Invdtmp=zeros(m,rr-ll+1);
            Invdtmp(1:rr-ll+1,1:rr-ll+1)=inv(o.memmapD.Data.M(ll:rr,1:rr-ll+1));
            fwrite(fileIDInvd,Invdtmp,'double');
        else
            o.Invd{i}=inv(o.D{i});
        end     
    end;
    
    if o.profile
        t = toc(tLU0);
        fprintf('Direct inverse at last level takes %2.2e seconds\n',t);
        fprintf(logFileID,'-----------------------------\n');
        fprintf(logFileID,'Direct inverse at last level takes %2.2e seconds\n',t);
    end
    
    if(o.isOutofCore)
        clear Invdtmp;
        fclose(fileIDInvd);
        filename = [matFiles 'Invd.bin'];
        o.memmapInvd=memmapfile(filename, ...
          'Format', {'double' [m n] 'M'},     ...
          'Writable', true);

        filenameInvz = [matFiles 'Invz.bin'];
        fileIDInvz = fopen(filenameInvz,'a');
        filenameW = [matFiles 'W.bin'];
        fileIDW = fopen(filenameW,'a');
    end
    
    for lev=o.lev_max:-1:1
        fprintf(logFileID,'At level %d of factorize\n',lev);
        if o.verbose
        fprintf('At level %d of factorize\n',lev);
        end
        
        p_ind(1)=0;  p_ind(2^(lev-1)+1)=n;
        p_ind(2:2^(lev-1))=sort(o.partition(1:2^(lev-1)-1));
        jstart=1+sum(o.l(1:lev-1));
        jend=jstart+o.l(lev)-1;
        
        if(o.isOutofCore)        
            Wtmp=zeros(n,o.l(lev));
            Invztmp=zeros(2*o.lmax,2*o.l(lev));
        end

        if o.profile
            tfact1 = tic; 
        end
        
        for i=1:2^(lev-1)
            ll=p_ind(i)+1;  rr=p_ind(i+1);  xx=o.partition(2^(lev-1)+i-1);
            istart=1+2*(sum(cols_lev(1:lev-1))+o.l(lev)*(i-1));
            iend=istart+2*o.l(lev)-1;
            if(o.isOutofCore)
                
                if o.profile
                    tload0 = tic; 
                end
        
                Q1=o.memmapQ.Data.M(ll:xx,jstart:jend);
                Q2=o.memmapQ.Data.M(xx+1:rr,jstart:jend);
                R1=o.memmapR.Data.M(ll:rr-xx+ll-1,jstart:jend);
                R2=o.memmapR.Data.M(rr-xx+ll:rr,jstart:jend);
                if o.profile
                    t = toc(tload0);
                    fprintf('Retrieving Q and R from memory takes %2.2e seconds\n',t);
                    fprintf(logFileID,'-----------------------------\n');
                    fprintf(logFileID,'Retrieving Q and R from memory takes %2.2e seconds\n',t);
                end
                
                if o.profile
                    tsolve0 = tic; 
                end
                
                Wtmp(ll:xx,:)=o.solve(Q1,n,lev,2*i-1);
                Wtmp(xx+1:rr,:)=o.solve(Q2,n,lev,2*i);
                
                if o.profile
                    t = toc(tsolve0);
                    fprintf('Solves take %2.2e seconds\n',t);
                    fprintf(logFileID,'-----------------------------\n');
                    fprintf(logFileID,'Solves take %2.2e seconds\n',t);
                end
                
                if o.profile
                    tassign0 = tic; 
                end
                
                Wp(1:rr-ll+1,1:2*o.l(lev))=[zeros(xx-ll+1,o.l(lev)) Wtmp(ll:xx,:);...
                                    Wtmp(xx+1:rr,:) zeros(rr-xx,o.l(lev))];
                Vp(1:2*o.l(lev),1:rr-ll+1)=[R2' zeros(o.l(lev),rr-xx);...
                                     zeros(o.l(lev),xx-ll+1) R1'];
                if o.profile
                    t = toc(tassign0);
                    fprintf('Assignment in SMWF takes %2.2e seconds\n',t);
                    fprintf(logFileID,'-----------------------------\n');
                    fprintf(logFileID,'Assignment in SMW takes %2.2e seconds\n',t);
                end
                
                if o.profile
                    tsmw0 = tic; 
                end
                
                Invztmp(1:2*o.l(lev),:)=inv(eye(2*o.l(lev))+Vp(1:2*o.l(lev),1:rr-ll+1)*Wp(1:rr-ll+1,1:2*o.l(lev))); 

                if o.profile
                    t = toc(tsmw0);
                    fprintf('LU of SMWF takes %2.2e seconds\n',t);
                    fprintf(logFileID,'-----------------------------\n');
                    fprintf(logFileID,'LU of SMWF takes %2.2e seconds\n',t);
                end
                
                if o.profile
                    twrite0 = tic; 
                end
                
                fwrite(fileIDInvz,Invztmp,'double');
                
                if o.profile
                    t = toc(twrite0);
                    fprintf('Writing factors to file takes %2.2e seconds\n',t);
                    fprintf(logFileID,'-----------------------------\n');
                    fprintf(logFileID,'Writing factors to file takes %2.2e seconds\n',t);
                end
            else
                o.W{2^lev - 2 + 2 * i - 1}=o.solve(o.Q{2^lev - 2 + 2 * i - 1},n,lev,2*i-1);
                o.W{2^lev - 2 + 2 * i - 0}=o.solve(o.Q{2^lev - 2 + 2 * i - 0},n,lev,2*i);
                Wp(1:rr-ll+1,1:2*o.l(lev))=[zeros(xx-ll+1,o.l(lev)) o.W{2^lev - 2 + 2 * i - 1};...
                                    o.W{2^lev - 2 + 2 * i - 0} zeros(rr-xx,o.l(lev))];
                Vp(1:2*o.l(lev),1:rr-ll+1)=[o.R{2^lev - 2 + 2 * i - 0} zeros(o.l(lev),rr-xx);...
                                     zeros(o.l(lev),xx-ll+1) o.R{2^lev - 2 + 2 * i - 1}];
                o.Invz{2^(lev - 1) - 1 + i}=inv(eye(2*o.l(lev))+Vp(1:2*o.l(lev),1:rr-ll+1)*Wp(1:rr-ll+1,1:2*o.l(lev)));         
            end;
        end; 
        if(o.isOutofCore)
            fwrite(fileIDW,Wtmp,'double');

            o.memmapW=memmapfile(filenameW, ...
              'Format', {'double' [n sum(o.l(o.lev_max:-1:lev))] 'M'},     ...
              'Writable', true);
          
            o.memmapInvz=memmapfile(filenameInvz, ...
              'Format', {'double' [2*o.lmax 2*sum(cols_lev(o.lev_max:-1:lev))] 'M'},     ...
              'Writable', true);
        end
        
        if o.profile
            t = toc(tfact1);
            fprintf('Factorize at level %d takes %2.2e seconds\n',lev,t);
            fprintf(logFileID,'-----------------------------\n');
            fprintf(logFileID,'Factorize at level %d takes %2.2e seconds\n',lev,t);
        end
    end;
    if(o.isOutofCore)
        fclose(fileIDW);        
        fclose(fileIDInvz);
    end
    if o.profile
        t = toc(tfact0);
        fprintf('Factorize takes %2.2e seconds\n',t);
        fprintf(logFileID,'-----------------------------\n');
        fprintf(logFileID,'Factorize takes %2.2e seconds\n',t);
    end
    fclose(logFileID);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y1,Y2] = AtimesG(o,O1,O2,p_ind,lcurrent,lev,n)
    if o.verbose
    fprintf('Inside AtimesG\n');  
    logFileID = fopen(o.logFile,'a'); 
    fprintf(logFileID,'Inside AtimesG\n');
    end
    
    if o.profile
        tinit0 = tic; 
    end
    
    for i = 1:(2^(lev-1))
        ll = p_ind(i)+1;  rr = p_ind(i+1);  xx = o.partition(2^(lev-1)+i-1);
        O1(ll:xx,1:lcurrent) = normrnd(0,1,xx-ll+1,lcurrent);
        O1(xx+1:rr,1:lcurrent) = zeros(rr-xx,lcurrent);
        O2(ll:xx,1:lcurrent) = zeros(xx-ll+1,lcurrent);
        O2(xx+1:rr,1:lcurrent) = normrnd(0,1,rr-xx,lcurrent);
    end;
    
    Y1 = zeros(n,lcurrent);
    Y2 = zeros(n,lcurrent);
    
    if o.profile
        t = toc(tinit0);
        fprintf('Generating random matrices takes %2.2e seconds\n',t);
        fprintf(logFileID,'-----------------------------\n');
        fprintf(logFileID,'Generating random matrices takes %2.2e seconds\n',t);
    end

    if(o.isOutofCore)
        maxbsize = floor(1e9*o.memsize/(8*n));
        bsize = min(n,maxbsize);
        if o.profile
            fprintf('Block size = %d\n',bsize);
            fprintf(logFileID,'Block size = %d\n',bsize);
        end
        for i = 1:floor(n/bsize)

            istart=1+(i-1)*bsize;
            iend=i*bsize;
            if o.profile
                tbl0 = tic; 
            end
            
            row = o.memmapL.Data.L(istart:iend,:);
            
            if o.profile
                t = toc(tbl0);
                fprintf('Retriving block %d takes %2.2e seconds\n',i,t);
                fprintf(logFileID,'-----------------------------\n');
                fprintf(logFileID,'Retriving block %d takes %2.2e seconds\n',i,t);
            end
            
            if o.profile
                tmult0 = tic; 
            end
            Y1(istart:iend,:) = row*O2(:,1:lcurrent);
            Y2(istart:iend,:) = row*O1(:,1:lcurrent);
            
            if o.profile
                t = toc(tmult0);
                fprintf('Multiplication of block %d takes %2.2e seconds\n',i,t);
                fprintf(logFileID,'-----------------------------\n');
                fprintf(logFileID,'Multiplication of block %d takes %2.2e seconds\n',i,t);
            end
        end;
        
        if o.profile
            tmult0 = tic; 
        end

        istart=1+i*bsize;
        row = o.memmapL.Data.L(istart:end,:);
        Y1(istart:end,:) = row*O2(:,1:lcurrent);
        Y2(istart:end,:) = row*O1(:,1:lcurrent);
        
        if o.profile
            t = toc(tmult0);
            fprintf('Multiplication of block %d takes %2.2e seconds\n',i+1,t);
            fprintf(logFileID,'-----------------------------\n');
            fprintf(logFileID,'Multiplication of block %d takes %2.2e seconds\n',i+1,t);
        end
        clear row; 
    else
        Y1 = o.L*O2(:,1:lcurrent);
        Y2 = o.L*O1(:,1:lcurrent);    
    end;
    
    if o.profile
        tmult0 = tic; 
    end

    Y1 = Y1-o.mult(O2(:,1:lcurrent),lev-1,n,0);     % Y1 = A12 * G2
    Y2 = Y2-o.mult(O1(:,1:lcurrent),lev-1,n,0);     % Y2 = A21 * G1   

    if o.profile
        t = toc(tmult0);
        fprintf('Subtracting previous levels of approximation from mult takes %2.2e seconds\n',t);
        fprintf(logFileID,'-----------------------------\n');
        fprintf(logFileID,'Subtracting previous levels of approximation from mult takes %2.2e seconds\n',t);
    end
    if o.verbose || o.profile
    fclose(logFileID);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y1,Y2] = AtransposetimesG(o,O1,O2,p_ind,lcurrent,lev,n)
    Y1 = zeros(n,lcurrent);
    Y2 = zeros(n,lcurrent);
    if(o.isOutofCore)
        maxbsize = floor(1e9*o.memsize/(8*n));
        bsize = min(n,maxbsize);
        for i = 1:floor(n/bsize)
            istart = 1+(i-1)*bsize;
            iend = i*bsize;
            row = o.memmapL.Data.L(:,istart:iend)';
            Y1(istart:iend,:) = row*O2(:,1:lcurrent);
            Y2(istart:iend,:) = row*O1(:,1:lcurrent);
        end;
        istart = 1+i*bsize;
        row = o.memmapL.Data.L(:,istart:end)';
        Y1(istart:end,:) = row*O2(:,1:lcurrent);
        Y2(istart:end,:) = row*O1(:,1:lcurrent);
        clear row;
    else
        Y1 = o.L'*O2(:,1:lcurrent);
        Y2 = o.L'*O1(:,1:lcurrent);
    end;
    Y1 = Y1-o.mult_transpose(O2(:,1:lcurrent),lev-1,n,0);     % Y1 = A12 * G2
    Y2 = Y2-o.mult_transpose(O1(:,1:lcurrent),lev-1,n,0);     % Y2 = A21 * G1   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eguess = check_error(o,Qtmp,Y1,Y2,n,lev,lcurrent,ltest,normA) 
    eguess = 0;
    p_ind = zeros(2^(lev-1)+1,1);
    p_ind(1) = 0;  p_ind(2^(lev-1)+1) = n;
    p_ind(2:2^(lev-1)) = sort(o.partition(1:2^(lev-1)-1));    
    for i = 1:(2^(lev-1))
        ll = p_ind(i)+1;  rr = p_ind(i+1);  xx = o.partition(2^(lev-1)+i-1);
        for j = 1:ltest
            eguess = max(norm(Qtmp(ll:xx,1:lcurrent)*(Qtmp(ll:xx,1:lcurrent)'*Y1(ll:xx,j))-Y1(ll:xx,j))/(norm(Y1(ll:xx,j))*normA(2*i-1)),eguess);
            eguess = max(norm(Qtmp(xx+1:rr,1:lcurrent)*(Qtmp(xx+1:rr,1:lcurrent)'*Y2(xx+1:rr,j))-Y2(xx+1:rr,j))/(norm(Y2(xx+1:rr,j))*normA(2*i)),eguess);
        end;
    end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=mult(o,x,lev_max,n,dflag)

    y=zeros(n,size(x,2));
    p_ind=zeros(2^lev_max+1,1);
   
%     maxbsize = floor(1e9*o.memsize/(8*n));
%     bsize = min(n,maxbsize);
%     Qtmp = o.memmapQ.Data.M(:,1:maxbsize);

    for lev=1:lev_max
        p_ind(1)=0;  p_ind(2^(lev-1)+1)=n;
        p_ind(2:2^(lev-1))=sort(o.partition(1:2^(lev-1)-1));
        jstart=1+sum(o.l(1:lev-1));
        jend=jstart+o.l(lev)-1;        
        for i=1:2^(lev-1)
            ll=p_ind(i)+1;  rr=p_ind(i+1);  xx=o.partition(2^(lev-1)+i-1);
            if(o.isOutofCore)
                y(ll:xx,:)=y(ll:xx,:)+o.memmapQ.Data.M(ll:xx,jstart:jend)*(o.memmapR.Data.M(ll:rr-xx+ll-1,jstart:jend)'*x(xx+1:rr,:)); 
            else
                y(ll:xx,:)=y(ll:xx,:)+o.Q{2^lev - 2 + 2 * i - 1}*(o.R{2^lev - 2 + 2 * i - 1}*x(xx+1:rr,:)); 
            end
        end
        for i=1:2^(lev-1)
            ll=p_ind(i)+1;  rr=p_ind(i+1);  xx=o.partition(2^(lev-1)+i-1);
            if(o.isOutofCore)
                y(xx+1:rr,:)=y(xx+1:rr,:)+o.memmapQ.Data.M(xx+1:rr,jstart:jend)*(o.memmapR.Data.M(rr-xx+ll:rr,jstart:jend)'*x(ll:xx,:)); 
            else
                y(xx+1:rr,:)=y(xx+1:rr,:)+o.Q{2^lev - 2 + 2 * i - 0}*(o.R{2^lev - 2 + 2 * i - 0}*x(ll:xx,:)); 
            end
        end
    end
    lev=lev_max+1;
    p_ind(1)=0;  p_ind(2^(lev-1)+1)=n;
    p_ind(2:2^(lev-1))=sort(o.partition(1:2^(lev-1)-1));
    if(dflag)
        for i=1:2^lev_max
            ll=p_ind(i)+1;  rr=p_ind(i+1);
            if(o.isOutofCore)
                y(ll:rr,:)=y(ll:rr,:)+o.memmapD.Data.M(ll:rr,1:rr-ll+1)*x(ll:rr,:);
            else
                y(ll:rr,:)=y(ll:rr,:)+o.D{i}*x(ll:rr,:);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=mult_transpose(o,x,lev_max,n,dflag)
   y=zeros(n,size(x,2));
   p_ind=zeros(2^lev_max+1,1);
   for lev=1:lev_max
        p_ind(1)=0;  p_ind(2^(lev-1)+1)=n;
        p_ind(2:2^(lev-1))=sort(o.partition(1:2^(lev-1)-1));
        jstart=1+sum(o.l(1:lev-1));
        jend=jstart+o.l(lev)-1;      
        for i=1:2^(lev-1)
            ll=p_ind(i)+1;  rr=p_ind(i+1);  xx=o.partition(2^(lev-1)+i-1);
            if(o.isOutofCore)
                y(ll:xx,:)=y(ll:xx,:)+o.memmapR.Data.M(rr-xx+ll:rr,jstart:jend)*(o.memmapQ.Data.M(xx+1:rr,jstart:jend)'*x(xx+1:rr,:)); 
            else
                y(ll:xx,:)=y(ll:xx,:)+o.R{2^lev - 2 + 2 * i - 0}'*(o.Q{2^lev - 2 + 2 * i - 0}'*x(xx+1:rr,:)); 
            end
        end;
        for i=1:2^(lev-1)
            ll=p_ind(i)+1;  rr=p_ind(i+1);  xx=o.partition(2^(lev-1)+i-1);
            if(o.isOutofCore)
                y(xx+1:rr,:)=y(xx+1:rr,:)+o.memmapR.Data.M(ll:rr-xx+ll-1,jstart:jend)*(o.memmapQ.Data.M(ll:xx,jstart:jend)'*x(ll:xx,:)); 
            else
                y(xx+1:rr,:)=y(xx+1:rr,:)+o.R{2^lev - 2 + 2 * i - 1}'*(o.Q{2^lev - 2 + 2 * i - 1}'*x(ll:xx,:)); 
            end
        end;
    end;
    lev=lev_max+1;
    p_ind(1)=0;  p_ind(2^(lev-1)+1)=n;
    p_ind(2:2^(lev-1))=sort(o.partition(1:2^(lev-1)-1));
    if(dflag)
        for i=1:2^lev_max
            ll=p_ind(i)+1;  rr=p_ind(i+1);
            if(o.isOutofCore)
                y(ll:rr,:)=y(ll:rr,:)+o.memmapD(ll:rr,1:rr-ll+1)'*x(ll:rr,:);
            else
                y(ll:rr,:)=y(ll:rr,:)+o.D{i}'*x(ll:rr,:);
            end
        end;
    end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w=solve(o,u,n,lev,k)
    if o.verbose
    logFileID = fopen(o.logFile,'a');
    fprintf('Inside solve\n');  
    fprintf(logFileID,'Inside solve\n');
    end

    if o.profile
        tsolve0 = tic; 
    end
    
%     invLU = @(L,U,w) (U\(L\w));
	cols_lev = zeros(o.lev_max,1);
    for i = 1:o.lev_max
        cols_lev(i) = o.l(i)*2^(i-1);
    end
    j=o.lev_max;
    w = zeros(size(u,1),size(u,2));
%     Wp=zeros(n,2*l);
%     Vp=zeros(2*l,n);
%    invZVw = zeros(2*o.lmax,size(u,2));
    p_ind = zeros(2^o.lev_max+1,1);
    p_ind(1) = 0;  p_ind(2^j+1) = n;
    p_ind(2:2^j) = sort(o.partition(1:2^j-1));
    ss = p_ind((k-1)*2^(j-lev)+1);
    for i=1:2^(j-lev)
        ll = p_ind(i+(k-1)*2^(j-lev))+1;  rr = p_ind(i+1+(k-1)*2^(j-lev));
        if(o.isOutofCore)

            if o.profile
                tLUd0 = tic; 
            end

            Invd1 = o.memmapInvd.Data.M(1:rr-ll+1,ll:rr);

            if o.profile
                t = toc(tLUd0);
                fprintf('Retreiving Invd takes %2.2e seconds\n',t);
                fprintf(logFileID,'-----------------------------\n');
                fprintf(logFileID,'Retreiving Invd takes %2.2e seconds\n',t);
            end

            if o.profile
                tInvdapply0 = tic; 
            end
            w(ll-ss:rr-ss,:) = Invd1*u(ll-ss:rr-ss,:);        % inv(D)*u
            if o.profile
                t = toc(tInvdapply0);
                fprintf('Applying Invd takes %2.2e seconds\n',t);
                fprintf(logFileID,'-----------------------------\n');
                fprintf(logFileID,'Applying Invd takes %2.2e seconds\n',t);
            end            
        else

            if o.profile
                tInvdapply0 = tic; 
            end
            w(ll-ss:rr-ss,:) = o.Invd{(k - 1) * 2^(j - lev) + i}*u(ll-ss:rr-ss,:);   % inv(D)*u     
            if o.profile
                sz = size(o.Invd{(k - 1) * 2^(j - lev) + i},1);
                t = toc(tInvdapply0);
                fprintf('Applying Invd of size %d takes %2.2e seconds\n',sz,t);
                fprintf(logFileID,'-----------------------------\n');
                fprintf(logFileID,'Applying Invd of size %d takes %2.2e seconds\n',sz,t);
            end      
        end
    end;  
    if(o.isOutofCore)
        clear Ld1 Ud1;
    end
    for j=o.lev_max-1:-1:lev
        if o.verbose
        fprintf(logFileID,'At level %d of solve\n',j);
        fprintf('At level %d of solve\n',j);
        end
        p_ind(1)=0;  p_ind(2^j+1)=n;
        p_ind(2:2^j)=sort(o.partition(1:2^j-1));
        ss=p_ind((k-1)*2^(j-lev)+1);
        jstart=1 + sum(o.l(1:j));
        jend=jstart + o.l(j+1)-1;
        if(o.isOutofCore)
            jstartW = 1 + sum(o.l(o.lev_max:-1:j+2));
            jendW = jstartW + o.l(j+1)-1;
        end
        for i=1:2^(j-lev)
            ll=p_ind(i+(k-1)*2^(j-lev))+1;  rr=p_ind(i+1+(k-1)*2^(j-lev));  xx=o.partition(2^j+i+(k-1)*2^(j-lev)-1);
            if(o.isOutofCore)
                istart = 1+2*(sum(cols_lev(o.lev_max:-1:j+2))+o.l(j+1)*(k-1)*2^(j-lev)+o.l(j+1)*(i-1));
                iend = istart+2*o.l(j+1)-1;
                
                if o.profile
                    tRLUzW0 = tic; 
                end   
                R1 = o.memmapR.Data.M(ll:rr-xx+ll-1,jstart:jend);
                R2 = o.memmapR.Data.M(rr-xx+ll:rr,jstart:jend);
                Invz1 = o.memmapInvz.Data.M(1:2*o.l(j+1),istart:iend);
                W1 = o.memmapW.data.M(ll:xx,jstartW:jendW);
                W2 = o.memmapW.Data.M(xx+1:rr,jstartW:jendW);
                if o.profile
                    t = toc(tRLUzW0);
                    fprintf('Retreiving Invz, R, W takes %2.2e seconds\n',t);
                    fprintf(logFileID,'-----------------------------\n');
                    fprintf(logFileID,'Retreiving Invz, R, W takes %2.2e seconds\n',t);
                end   
                
                if o.profile
                    tRLUzWapply0 = tic; 
                end     
                
                invZVw(1:2*o.l(j+1),:) = Invz1*[R2'*w(ll-ss:xx-ss,:);...     % (I - W*inv(Z)*V)*w     W=inv(D)*U 
                                R1'*w(xx-ss+1:rr-ss,:)];
                w(ll-ss:rr-ss,:) = w(ll-ss:rr-ss,:) - [W1*invZVw(o.l(j+1)+1:2*o.l(j+1),:);W2*invZVw(1:o.l(j+1),:)];
                
                if o.profile
                    t = toc(tRLUzWapply0);
                    fprintf('Using Invz, R, W in SMWF takes %2.2e seconds\n',t);
                    fprintf(logFileID,'-----------------------------\n');
                    fprintf(logFileID,'Using Invz, R, W in SMWF takes %2.2e seconds\n',t);
                end   
            else
                istart = 1 + 2*(sum(cols_lev(1:j)) + o.l(j+1)*(k-1)*2^(j-lev)+o.l(j+1)*(i-1));
                iend = istart + 2*o.l(j+1) - 1;
                
                if o.profile
                    tInvzVwapply0 = tic; 
                end 
                invZVw = o.Invz{2^j - 1 + (k - 1) * 2^(j - lev) + i}*[o.R{2^(j + 1) - 2 + 2 * (k - 1) * 2^(j - lev) + 2 * i - 0}*w(ll-ss:xx-ss,:);...     % (I - W*inv(Z)*V)*w     W=inv(D)*U 
                                o.R{2^(j + 1) - 2 + 2 * (k - 1) * 2^(j-lev) + 2 * i - 1}*w(xx-ss+1:rr-ss,:)];
                            
                if o.profile
                    t = toc(tInvzVwapply0);
                    fprintf('Using Invz, R in SMWF takes %2.2e seconds\n',t);
                    fprintf(logFileID,'-----------------------------\n');
                    fprintf(logFileID,'Using Invz, R in SMWF takes %2.2e seconds\n',t);
                end 
                
                if o.profile
                    tWInvzVwapply0 = tic; 
                end

                w(ll-ss:rr-ss,:) = w(ll-ss:rr-ss,:) - [o.W{2^(j + 1) - 2 + 2 * (k - 1) * 2^(j - lev) + 2 * i - 1}*invZVw(o.l(j+1)+1:2*o.l(j+1),:);
                                                       o.W{2^(j + 1) - 2 + 2 * (k - 1) * 2^(j - lev) + 2 * i - 0}*invZVw(1:o.l(j+1),:)];                
                if o.profile
                    t = toc(tWInvzVwapply0);
                    fprintf('Using W in SMWF takes %2.2e seconds\n',t);
                    fprintf(logFileID,'-----------------------------\n');
                    fprintf(logFileID,'Using W in SMWF takes %2.2e seconds\n',t);
                end 
            end
        end;
    end;
            
    if o.profile
        t = toc(tsolve0);
        fprintf('solve takes takes %2.2e seconds\n',t);
        fprintf(logFileID,'-----------------------------\n');
        fprintf(logFileID,'solve takes takes %2.2e seconds\n',t);
    end
    
    if o.profile || o.verbose
    fclose(logFileID);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = solveReqs(o,memmapM11,memmapM12,n,k,matFiles,idebug)
    % Write the timings to screen and file
    logFileID = fopen(o.logFile,'a');

    
    nl = n-k;
    maxbsize = floor(1e9*o.memsize/(8*n));
    
    
    if o.profile
      tQR0 = tic; 
    end
    
    [o.QR,~] = qr(memmapM12.Data.M,0);
    
    if o.profile
      t = toc(tQR0);
      fprintf('QR of M12 takes %2.2e seconds\n',t);
      fprintf(logFileID,'-----------------------------\n');
      fprintf(logFileID,'QR of M12 takes %2.2e seconds\n',t);
    end
    
    oh = householder;
    if o.profile
      tHHVects0 = tic;
    end
    
    o.H = oh.nrqr(o.QR);
    
    if o.profile
      t = toc(tHHVects0);
      fprintf('Getting Householder vectors takes %2.2e seconds\n',t);
      fprintf(logFileID,'-----------------------------\n');
      fprintf(logFileID,'Getting Householder vectors takes %2.2e seconds\n',t);
    end

    if(o.isOutofCore)
        filename = [matFiles 'tmpL.bin'];
        fileID = fopen(filename,'w');
        bsize = min(n,maxbsize);
        for i = 1:floor(n/bsize)
            istart=1+(i-1)*bsize;
            iend=i*bsize;
            col = memmapM11.Data.M(:,istart:iend);
            col = oh.applypreQNt(col,o.H);
            fwrite(fileID,col,'double');
        end;
        istart = 1+i*bsize;
        col = memmapM11.Data.M(:,istart:end);
        col = oh.applypreQNt(col,o.H);
        fwrite(fileID,col,'double');
        fclose(fileID);
        clear col;
        memmaptmpL = memmapfile(filename, ...
          'Format', {'double' [n n] 'tmpL'},     ...
          'Writable', true);

        for i = 1:floor(n/bsize)
            istart = 1+(i-1)*bsize;
            iend = i*bsize;
            row = memmaptmpL.Data.tmpL(istart:iend,:);
            row = oh.applypostQN(row,o.H);
            memmaptmpL.Data.tmpL(istart:iend,:) = row;
        end;
        istart=1+i*bsize;
        row = memmaptmpL.Data.tmpL(istart:end,:);
        row = oh.applypostQN(row,o.H);
        memmaptmpL.Data.tmpL(istart:end,:) = row;
        clear row;
        

        filename = [matFiles 'L.bin'];
        fileID = fopen(filename,'w');
        bsize = min(nl,maxbsize);
        for i = 1:floor(nl/bsize)
            istart = 1+k+(i-1)*bsize;
            iend = k+i*bsize;
            col = memmaptmpL.Data.tmpL(k+1:n,istart:iend);
            fwrite(fileID,col,'double');
        end;
        istart = 1+k+i*bsize;
        col = memmaptmpL.Data.tmpL(k+1:n,istart:end);
        fwrite(fileID,col,'double');    
        fclose(fileID);
        
        o.memmapL = memmapfile(filename, ...
          'Format', {'double' [nl nl] 'L'},     ...
          'Writable', true);
      
        tmp = zeros(n,k);
        bsize = min(n,maxbsize);
        for i = 1:floor(n/bsize)
            istart = 1+(i-1)*bsize;
            iend = i*bsize; 
            row = memmapM11.Data.M(istart:iend,:);
            tmp(istart:iend,:)=row*o.QR;
        end;
        istart = 1+i*bsize;
        row = memmapM11.Data.M(istart:end,:);
        tmp(istart:end,:)=row*o.QR;   
        S = o.QR'*tmp;

        o.LS = oh.applypreQNt(tmp,o.H);
        o.LS = o.LS(k+1:n,:);
        tmp = zeros(k,n);
        for i = 1:floor(n/bsize)
            istart = 1+(i-1)*bsize;
            iend = i*bsize; 
            col = memmapM11.Data.M(:,istart:iend);
            tmp(:,istart:iend) = o.QR'*col;
        end;
        istart = 1+i*bsize;
        col = memmapM11.Data.M(:,istart:end);
        tmp(:,istart:end) = o.QR'*col;

        o.SL = oh.applypostQN(tmp,o.H);
        o.SL = o.SL(:,k+1:n);
        clear tmp;
        clear col;
    else
        M11 = memmapM11.Data.M;  

        o.L = oh.applypreQNt(M11,o.H);
        o.L = oh.applypostQN(o.L,o.H);
        o.L = o.L(k+1:end,k+1:end);

        tmp = M11*o.QR;
        S = o.QR'*tmp;
        o.LS = oh.applypreQNt(tmp,o.H);
        o.LS = o.LS(k+1:n,:);
    
        tmp = o.QR'*M11;
        o.SL = oh.applypostQN(tmp,o.H);
        o.SL = o.SL(:,k+1:n);
        clear M11;
    end
    

    o = o.simple_partition(nl);
    
    if(o.isOutofCore)
        filename = [matFiles 'Q.bin'];
        fileID = fopen(filename,'w');
        fclose(fileID);
        filename = [matFiles 'R.bin'];
        fileID = fopen(filename,'w');
        fclose(fileID);
        filename = [matFiles 'D.bin'];
        fileID = fopen(filename,'w');
        fclose(fileID);
        filename = [matFiles 'Invd.bin'];
        fileID = fopen(filename,'w');
        fclose(fileID);
        filename = [matFiles 'Invz.bin'];
        fileID = fopen(filename,'w');
        fclose(fileID);
        filename = [matFiles 'W.bin'];
        fileID = fopen(filename,'w');
        fclose(fileID);
    end 
    
    if o.profile
        tHODLR0 = tic;
    end

    [o] = o.hodlr_no_FMM(matFiles);
    if o.profile
        t = toc(tHODLR0);
        fprintf('HODLR takes %2.2e seconds\n',t);
        fprintf(logFileID,'-----------------------------\n');
        fprintf(logFileID,'HODLR takes %2.2e seconds\n',t);
        tFact0 = tic;
    end
    o = o.factorize(nl,matFiles);
    if o.profile
        t = toc(tFact0);
        fprintf('Factorization takes %2.2e seconds\n',t);
        fprintf(logFileID,'-----------------------------\n');
        fprintf(logFileID,'Factorization takes %2.2e seconds\n',t);
        tic;
    end

    LLS = o.solve(o.LS,nl,0,1);
    tmp = [zeros(k,size(LLS,2));LLS];
%     o.invS1 = inv(o.QR'*memmapM12.Data.M);
    o.invS1 = (o.QR'*memmapM12.Data.M)\eye(k);
    
    M21 = o.memmapM21.Data.M;
    M22 = o.memmapM22.Data.M;
%     o.invS2 = inv(M21*(o.QR-oh.applypreQN(tmp,o.H))+M22*(o.invS1*(o.SL*LLS-S)));
    o.invS2 = (M21*(o.QR-oh.applypreQN(tmp,o.H))+M22*(o.invS1*(o.SL*LLS-S)))\eye(k);
    
    clear tmp;
    
    if o.profile
        t = toc;
        fprintf('Computing invS1 and invS2 takes %2.2e seconds\n',t);
        fprintf(logFileID,'-----------------------------\n');
        fprintf(logFileID,'Computing invS1 and invS2 takes %2.2e seconds\n',t);
        fclose(logFileID);
    end
    
    if idebug
        bsize = min(nl,maxbsize);
        errhodlr = 0;
        for i = 1:10
          randvec = rand(nl,1);
          y1 = o.mult(randvec,o.lev_max,nl,1);
          y2 = zeros(nl,1);
          if(o.isOutofCore)
              for j = 1:floor(nl/bsize)
                  jstart = 1+(j-1)*bsize;
                  jend = j*bsize; 
                  row = o.memmapL.Data.L(jstart:jend,:);
                  y2(jstart:jend,1) = row*randvec;
              end;
              jstart = 1+j*bsize; 
              row = o.memmapL.Data.L(jstart:end,:);
              y2(jstart:end,1) = row*randvec;
          else
              y2 = o.L*randvec;
          end;
          if errhodlr<norm(y1-y2)/norm(y2)
              errhodlr = norm(y1-y2)/norm(y2);
          end
        end
    
        errfact = 0;
        for i = 1:10
          randvec = rand(nl,1);
          y1 = zeros(nl,1);
          if(o.isOutofCore)
              for j=1:floor(nl/bsize)
                  jstart = 1+(j-1)*bsize;
                  jend = j*bsize; 
                  row = o.memmapL.Data.L(jstart:jend,:);
                  y1(jstart:jend,1) = row*randvec;
              end;
              jstart = 1+j*bsize;
              row = o.memmapL.Data.L(jstart:end,:);
              y1(jstart:end,1) = row*randvec;
          else
              y1 = o.L*randvec;
          end;
          y2 = o.solve(y1,nl,0,1);
          if errfact<norm(y2-randvec)/norm(randvec)
              errfact = norm(y2-randvec)/norm(randvec);
          end
        end
        
        if ~isempty(o.logFile)
          fid = fopen(o.logFile,'a');
          fprintf(fid,'***********************************************\n');
          fprintf(fid,'HODLR is used, errors and timings are reported:\n');
          fprintf(fid,'Error in hodlr:                %4.2e \n',errhodlr);
          fprintf(fid,'Error in factorization:        %4.2e \n',errfact);
          fprintf(fid,'***********************************************\n');

          fprintf('***********************************************\n');
          fprintf('HODLR is used, errors and timings are reported:\n');
          fprintf('Error in hodlr:                %4.2e \n',errhodlr);
          fprintf('Error in factorization:        %4.2e \n',errfact);
          fprintf('***********************************************\n');
          fclose(fid);
        else
          fprintf('***********************************************\n');
          fprintf('HODLR is used, errors and timings are reported:\n');
          fprintf('Error in hodlr:                %4.2e \n',errhodlr);
          fprintf('Error in factorization:        %4.2e \n',errfact);
          fprintf('***********************************************\n');
        end
    end
    
    % Save the matrices needed for the solve to the disk
    filename = [matFiles 'H.bin'];
    fileID = fopen(filename,'w');
    fwrite(fileID,size(o.H),'double');
    fwrite(fileID,o.H,'double');
    fclose(fileID);
    
    filename = [matFiles 'QR.bin'];
    fileID = fopen(filename,'w');
    fwrite(fileID,size(o.QR),'double');
    fwrite(fileID,o.QR,'double');
    fclose(fileID);
    
    filename = [matFiles 'invS2.bin'];
    fileID = fopen(filename,'w');
    fwrite(fileID,size(o.invS2),'double');
    fwrite(fileID,o.invS2,'double');
    fclose(fileID);
    
    filename = [matFiles 'invS1.bin'];
    fileID = fopen(filename,'w');
    fwrite(fileID,size(o.invS1),'double');
    fwrite(fileID,o.invS1,'double');
    fclose(fileID);
    
    filename = [matFiles 'SL.bin'];
    fileID = fopen(filename,'w');
    fwrite(fileID,size(o.SL),'double');
    fwrite(fileID,o.SL,'double');
    fclose(fileID);
    
    filename = [matFiles 'LS.bin'];
    fileID = fopen(filename,'w');
    fwrite(fileID,size(o.LS),'double');
    fwrite(fileID,o.LS,'double');
    fclose(fileID);
            
    filename = [matFiles 'l2.bin'];
    fileID = fopen(filename,'w');
    fwrite(fileID,size(o.l),'double');
    fwrite(fileID,o.l,'double');
    fclose(fileID);
    
    filename = [matFiles 'partition.bin'];
    fileID = fopen(filename,'w');
    fwrite(fileID,size(o.partition),'double');
    fwrite(fileID,o.partition,'double');
    fclose(fileID);
    
    if(o.isOutofCore)
        filename = [matFiles 'R.bin'];
        fileID = fopen(filename,'w');
        fwrite(fileID,size(o.R),'double');
        fwrite(fileID,o.R,'double');
        fclose(fileID);

        filename = [matFiles 'W.bin'];
        fileID = fopen(filename,'w');
        fwrite(fileID,size(o.W),'double');
        fwrite(fileID,o.W,'double');
        fclose(fileID);

        filename = [matFiles 'Invz.bin'];
        fileID = fopen(filename,'w');
        fwrite(fileID,size(o.Invz),'double');
        fwrite(fileID,o.Invz,'double');
        fclose(fileID);

        filename = [matFiles 'Invd.bin'];
        fileID = fopen(filename,'w');
        fwrite(fileID,size(o.Invd),'double');
        fwrite(fileID,o.Invd,'double');
        fclose(fileID);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = solveReqs_wMats(o,matFiles)
    
    % Load matrices from the disk to the memory
    filename = [matFiles 'H.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.H = memmap.Data.M(:,:);
    
    filename = [matFiles 'QR.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.QR = memmap.Data.M(:,:);
    
    filename = [matFiles 'invS2.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.invS2 = memmap.Data.M(:,:);
    
    filename = [matFiles 'invS1.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.invS1 = memmap.Data.M(:,:);
    
    filename = [matFiles 'SL.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.SL = memmap.Data.M(:,:);
    
    filename = [matFiles 'LS.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.LS = memmap.Data.M(:,:);
    
    filename = [matFiles 'R.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.R = memmap.Data.M(:,:);
    
    filename = [matFiles 'l.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.l = memmap.Data.M(:,:);
    
    filename = [matFiles 'partition.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.partition = memmap.Data.M(:,:);
    
    filename = [matFiles 'W.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.W = memmap.Data.M(:,:);
    
    filename = [matFiles 'Lz.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.Lz = memmap.Data.M(:,:);
    
    filename = [matFiles 'Uz.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.Uz = memmap.Data.M(:,:);
    
    filename = [matFiles 'Ld.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.Ld = memmap.Data.M(:,:);
    
    filename = [matFiles 'Ud.bin'];
    fileID = fopen(filename,'r');
    matSize = fread(fileID,2,'double');
    fclose(fileID);
    memmap = memmapfile(filename,'Offset',16,'Format',{'double',...
        [matSize(1) matSize(2)],'M'});
    o.Ud = memmap.Data.M(:,:);

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lhs=schur_complement(o,rhs,n,k)
    
    oh=householder;
    
    C=o.memmapM21.Data.M;
    F=o.memmapM22.Data.M;

    f1=rhs(1:n);
    f2=rhs(n+1:end);
    clear rhs;

    nl=n-k;

    r1=oh.applypreQNt(f1,o.H);
    r1=r1(k+1:n,:);
    r2=o.QR'*f1;
    Lr1=o.solve(r1,nl,0,1);
    tmp=[zeros(k,size(Lr1,2));Lr1];

    uR=o.invS2*(f2-(C*oh.applypreQN(tmp,o.H))+F*o.invS1*(o.SL*Lr1-r2));
    uN=r1-o.LS*uR;
    uN=o.solve(uN,nl,0,1);
    uN=[zeros(k,size(uN,2));uN];
    lambda=F\(f2-C*(o.QR*uR)-C*(oh.applypreQN(uN,o.H)));

    u=o.QR*uR+oh.applypreQN(uN,o.H);
    lhs=[u;lambda];    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = simple_partition(o,n)
    o.partition = zeros(2^o.lev_max-1,1);
    p_ind = zeros(2^o.lev_max+1,1);
    for j = 1:o.lev_max
        p_ind(1)=0;  p_ind(2^(j-1)+1)=n;
        p_ind(2:2^(j-1))=sort(o.partition(1:2^(j-1)-1));
        o.partition(2^(j-1):2^j-1,1)=floor((p_ind(2:2^(j-1)+1)+p_ind(1:2^(j-1)))/2);
    end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idWall2gids = pnts2qtree(~,XwallsInt,XwallsExt,idebug)
% idWall2gids = pnts2qtree(Xwalls,idebug) takes the points on the walls and
% reorder them based on Z-ordering. It constructs a quadtree and orders the
% points. The output idWall2gids(order of the points,2) returns what wall
% and what point on that wall the ordered point belongs to. This ordering
% is necessary for Fast Direct Solver and used only if the fast direct
% solver is requested by the user.
%
% This function requires the following functions:
% [] qtree.m
% [] morton_id.m

NbdInt = size(XwallsInt,1)/2; % number of points per wall
NbdExt = size(XwallsExt,1)/2;
nvbd   = size(XwallsInt,2)+1; % number of walls
n = NbdInt*(nvbd-1) + NbdExt; % number of points in total

% Shift Xwalls to (0,1)^2 domain b/c qtree code requires that
minExtx = min(XwallsExt(1:NbdExt));
minExty = min(XwallsExt(NbdExt+1:end));
minExt = min(minExtx,minExty);

XExtShift = zeros(size(XwallsExt));
XExtShift(1:NbdExt) = XwallsExt(1:NbdExt) -...
  (minExt-0.1);
XExtShift(NbdExt+1:end) = XwallsExt(NbdExt+1:end) - ...
  (minExt-0.1);

maxExtx = max(XExtShift(1:NbdExt));
maxExty = max(XExtShift(NbdExt+1:end));
maxExt = max(maxExtx,maxExty);

XExtNorm = zeros(size(XwallsExt));
XExtNorm(1:NbdExt) = XExtShift(1:NbdExt)/...
  (maxExt+0.1);
XExtNorm(NbdExt+1:end) = XExtShift(NbdExt+1:end)/...
  (maxExt+0.1);

XIntShift = zeros(size(XwallsInt));
XIntShift(1:NbdInt,:) = XwallsInt(1:NbdInt,:) -...
  (minExt-0.1);
XIntShift(NbdInt+1:end,:) = XwallsInt(NbdInt+1:end,:) - ...
  (minExt-0.1);

XIntNorm = zeros(size(XwallsInt));
XIntNorm(1:NbdInt,:) = XIntShift(1:NbdInt,:)/...
  (maxExt+0.1);
XIntNorm(NbdInt+1:end,:) = XIntShift(NbdInt+1:end,:)/...
  (maxExt+0.1);

xnormInt = XIntNorm(1:NbdInt,:);
ynormInt = XIntNorm(NbdInt+1:end,:);
xnorm = [XExtNorm(1:NbdExt);xnormInt(:)]; 
ynorm = [XExtNorm(NbdExt+1:end);ynormInt(:)];

% Order the points for qtree code
pnts = [xnorm(:)';ynorm(:)']; % in the form of (2,n)
maxPntsPerNode = 1;      % points per box
maxLevel       = 40;     % maximum tree depth

% Construct the class qtree
oq = qtree;
oq.insert_points(1:n,pnts,maxPntsPerNode,maxLevel);

% Get the leaves
leaves = oq.leaves();

% Assign the global ids to the orders
idx = 1;
gids_order = zeros(n,1);
for l = 1 : length(leaves)
  if ~isempty(leaves{l})
    ngids = numel(leaves{l}.gids);  
    gids_order(idx:idx+ngids-1) = leaves{l}.gids;
    idx = idx+ngids;
  end
end

% Find the ordered global ids corresponding wall and point on that wall
idWall2gids = zeros(n,2); % (n,1): which point,(on) (n,2): which wall
for id = 1 : n
  pointID = gids_order(id);
  if pointID <= NbdExt
    idWall2gids(id,1) = gids_order(id);
    idWall2gids(id,2) = 1;
  else
    pointID = pointID-NbdExt;
    pointNum = mod(pointID,NbdInt);
    if pointNum == 0
      idWall2gids(id,1) = pointNum+NbdInt;
    else
      idWall2gids(id,1) = pointNum;
    end
    idWall2gids(id,2) = ceil(pointID/NbdInt)+1;
  end 
end

% Plot the quadtree, walls, points, and show the ordering
if idebug

    figure(5); clf;
    oq.plottree;
    hold on
    plot(XExtNorm(1:NbdExt),XExtNorm(NbdExt+1:end),'b','linewidth',3)
    plot(XIntNorm(1:NbdInt,:),XIntNorm(NbdInt+1:end,:),'b','linewidth',3)
    axis equal
    plot(xnorm(:),ynorm(:),'g.','markersize',10)
    for id = 1 : n
        plot(xnorm(gids_order(id)),ynorm(gids_order(id)),'rs','markersize',15)
        pause(0.05)
    end
       
%     % Alternatively, on the actual grid
%     plot(Xwalls(1:end/2,:),Xwalls(end/2+1:end,:),'k','linewidth',2)
%     for l = 1 : n
%         plot(Xwalls(idWall2gids(l,1),idWall2gids(l,2)),...
%             Xwalls(idWall2gids(l,1)+Nbd,idWall2gids(l,2)),...
%             'rs','markersize',20) 
%         pause(0.1)
%     end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % methods
end % hodlr
    
