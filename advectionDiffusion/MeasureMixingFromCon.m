function mixing = MeasureMixingFromCon
fid = fopen('/scratch/gokberk/Case411/DiffIC/Case411_VesVel_Con_2Rings.bin','r');
val = fread(fid,6,'double');
N_radii = val(1)-1;
N_theta = val(2)-1;
ntime = val(3);
Th = val(4);
radius = [val(5) val(6)];

[~,~,r1D,theta,~,~,~] = generateGrid(N_theta,N_radii,radius);
rMat = repmat(r1D,1,N_theta);

fidWrite = fopen('/scratch/gokberk/Case411/DiffIC/TwoRingsVes_CG.bin','w');
fwrite(fidWrite,[N_radii+1;N_theta+1;ntime+1],'double');
fwrite(fidWrite,Th,'double');
fwrite(fidWrite,radius,'double');
fclose(fidWrite);
    
fidLog = fopen('/scratch/gokberk/Case411/DiffIC/TwoRingsVes_CGLog.log','w');
fclose(fidLog);
        
for tt = 1 : ntime
    C = fread(fid,val(1)*val(2),'double');
    C = reshape(C,val(1),val(2));
    
    if tt == 1
         mmean = 1/(pi*(radius(2)^2-radius(1)^2)) * trapz(theta,trapz(r1D,C(:,1:end-1).*rMat,1));
    end
   tic
   [CG] = MeasureMixingLaplacian(C(:,1:end-1)-mmean,N_radii,N_theta,radius,2);
   tim = toc;
   mixing(tt) = CG;
   
   writeData('/scratch/gokberk/Case411/DiffIC/TwoRingsVes_CG.bin',CG);
   keepDiary('/scratch/gokberk/Case411/DiffIC/TwoRingsVes_CGLog.log',ntime,tt,tim);

  
   % ---------------------------------------------------------------------
end
end

function writeData(fileName,dataArr)
fid = fopen(fileName,'a');
fwrite(fid,dataArr,'double');
fclose(fid);
end

function keepDiary(fileName,ntime,tst,tim)
fid = fopen(fileName,'a');
fprintf(fid,'%s\n',['Elapsed Time at step = ' num2str(tst) ' of ' num2str(ntime) ' = ' num2str(tim) ' sec']);

fclose(fid);
end