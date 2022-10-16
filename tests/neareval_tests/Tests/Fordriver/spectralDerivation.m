function zft=spectralDerivation(zf)
N=length(zf);
zft=fftshift(1i*[-N/2:N/2-1].*fftshift(zf));