function [x y deriv ] = getDomain(N,shape,shapeparam)

theta = (0:N-1)'*2*pi/N;
zf=zeros(1,N); 

if shape == 1
    r1 = shapeparam(1);
    r2 = shapeparam(2);
    zf(2) = 1/2*(r1+r2);
    zf(N) = 1/2*(r1-r2);
    
   
elseif shape == 2
    
    c3 = shapeparam(1);    

    zf=zeros(1,N); 
    zf(2) = 1;
    zf(4) = c3;
    
elseif shape == 3
    
    R = shapeparam(1);
    a = shapeparam(2);
    
    zf=zeros(1,N); 
    zf(2) = R;
    zf(a+2)=R/4;
    zf(N-a+2)=R/4;

end

z = N*ifft(zf); % Points of the shape in complex format
zft = spectralDerivation(zf); % Spectral derivation of zf
zt = N*ifft(zft);
x = real(z);
y = imag(z);
Dx= real(zt);
Dy = imag(zt);
deriv = sqrt(Dx.^2+Dy.^2);
