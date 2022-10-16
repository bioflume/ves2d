load stability.mat;
close all;
figure;
schemeDate = s1t2; %Choices s1t1 s1t2 s2t1

order = {'1','2','3','4'};
col = 'rgbcmyk';
mark = 'o*.xsd^v<>ph';
shearRate = {'.5','1','1.5','2','3','4','5'};
for i = 1:4
    ind = (schemeDate(:,1) == i);
    rang = 1+floor(7*rand-eps);
    marker = 1+floor(12*rand-eps);
    hold on; plot(schemeDate(ind,2),schemeDate(ind,5),...
        ['-' col(rang) mark(marker)]);
end
title('scheme 1 -- ts vs shear rate');
xlabel('shear rate');ylabel('ts');
legend(order{:});

figure;
for i = [.5 1 1.5 2 3 4 5]
    ind = (schemeDate(:,2) == i);
    rang = 1+floor(7*rand-eps);
    marker = 1+floor(12*rand-eps);
    hold on; plot(schemeDate(ind,1),schemeDate(ind,5),...
        ['-' col(rang) mark(marker)]);
end
title('scheme 1 -- ts vs order');
xlabel('order');ylabel('ts');
legend(shearRate{:});
