set(0,'DefaultAxesFontSize',24)
set(0,'DefaultAxesFontName', 'Helvetica')

% CASE = '/scratch/gokberk/PeNormTest_Con.bin';
CASE = '/scratch/gokberk/Case211/NBC/LayerIC/Case211_VesVel_Con_Pe4.bin';
fileNameV{1} = [CASE];
% fileNameV{1} = [CASE '_VesVel_Con_PeInf.bin'];
fileNameV{2} = [CASE '_VesVel_Con_Pe2.bin'];
fileNameV{3} = [CASE '_VesVel_Con_Pe3.bin'];
fileNameV{4} = [CASE '_VesVel_Con_Pe4.bin'];
fileNameV{5} = [CASE '_VesVel_Con_Pe5.bin'];
fileNameV{6} = [CASE '_VesVel_Con_Pe6.bin'];
fileNameV{7} = [CASE '_VesVel_Con_Pe7.bin'];
fileNameV{8} = [CASE '_VesVel_Con_Pe8.bin'];

fileNameA{1} = [CASE '_VesVel_Con_Pe1.bin'];
fileNameA{2} = [CASE '_AntVel_Con_Pe2.bin'];
fileNameA{3} = [CASE '_AntVel_Con_Pe3.bin'];
fileNameA{4} = [CASE '_AntVel_Con_Pe4.bin'];
fileNameA{5} = [CASE '_AntVel_Con_Pe5.bin'];
fileNameA{6} = [CASE '_AntVel_Con_Pe6.bin'];
fileNameA{7} = [CASE '_AntVel_Con_Pe7.bin'];
fileNameA{8} = [CASE '_AntVel_Con_Pe8.bin'];

% Read File
fid = fopen(fileNameV{1},'r');
val = fread(fid,6,'double');
fclose(fid);
N_radii = val(1)-1;
N_theta = val(2)-1;
ntime   = val(3);
Th      = val(4);
radius  = [val(5) val(6)];

% Generate Space and Time Domains
[x,y,r1D,theta,~,~,~] = generateGrid(N_theta,N_radii,radius);
rMat = repmat(r1D,1,N_theta);
time = linspace(0,Th,ntime)';

% PDF data
pdfPoints = 401;
datCon = linspace(0,1,pdfPoints);
pdfCon_Ves = zeros(pdfPoints,ntime);
pdfCon_Ant = zeros(pdfPoints,ntime);

mmean_Ves = zeros(1,ntime);
CM_Ves    = zeros(3,ntime);

mmean_Ant = zeros(1,ntime);
CM_Ant    = zeros(3,ntime);
% CM(1) = 2nd moment
% CM(2) = 3rd moment
% CM(3) = 4th moment

CMDIFF = zeros(3,ntime,8);
% Plot PDF?
plotPDF = 0;

for i = 1 : 1
    fidVes = fopen(fileNameV{i},'r');
    val1 = fread(fidVes,6,'double');
    
%     fidAnt = fopen(fileNameA{i},'r');
%     val2 = fread(fidAnt,6,'double');

    tic
    for tt = 1 : ntime
        % Get the concentration
        C1 = fread(fidVes,val1(1)*val1(2),'double');
        C1_mat = reshape(C1,val1(1),val1(2));
        C1_mat = C1_mat(:,1:end-1);
        
%         C2 = fread(fidAnt,val2(1)*val2(2),'double');
%         C2_mat = reshape(C2,val2(1),val2(2));
%         C2_mat = C2_mat(:,1:end-1);
        
        % Compute the mean
        mmean_Ves(tt) = 1/(pi*(radius(2)^2-radius(1)^2)) * trapz(theta,trapz(r1D,C1_mat.*rMat,1));
%         mmean_Ant(tt) = 1/(pi*(radius(2)^2-radius(1)^2)) * trapz(theta,trapz(r1D,C2_mat.*rMat,1));

        % Compute moments
        for mom = 1 : 3
            CM_Ves(mom,tt) = 1/(pi*(radius(2)^2-radius(1)^2)) * trapz(theta,trapz(r1D,(C1_mat).^(mom+1).*rMat,1));
%             CM_Ant(mom,tt) = 1/(pi*(radius(2)^2-radius(1)^2)) * trapz(theta,trapz(r1D,(C2_mat-mmean_Ant(tt)).^(mom+1).*rMat,1));
            
%             CMDIFF(mom,:,i) = CM_Ant(mom,:)./CM_Ves(mom,:);
        end
        
        
        % Get the pdf (do this for each Pe number separately)
        if plotPDF
            pdfCon_Ves(:,tt) = ksdensity(C1_mat(:),datCon)/pdfPoints;
            pdfCon_Ant(:,tt) = ksdensity(C2_mat(:),datCon)/pdfPoints;
            
            plot(datCon,pdfCon_Ves(:,tt),'b','linewidth',4);
            hold on
            plot(datCon,pdfCon_Ant(:,tt),':r','linewidth',4);
            legend('w/ vesicles', 'w/o vesicles');
            pause(0.1)
            hold off
        end

    end
    toc
end
