clc;
color= 'rgbk'
for lambda = 4:4:16
    for lobe = [0 2:5]
        n = 128;
        [g X] = sampleBd(2,n,1,'spectrum',lobe);
        % hold on; plotVec(X); grid on; axis equal
        D = kernelD(X(:));
        D = eye(2*n)-2*(1-lambda)*D/(1+lambda);
        L = sort(abs(eig(D)));
        hold on; plot(L,color(lambda/4))
        cond(D)
    end
end
