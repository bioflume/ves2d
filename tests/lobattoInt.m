function fint = lobattoInt(f)

N = size(f,1)/2;
ntime = size(f,3);
fint = zeros(2*N,1,ntime);


if ntime == 2
  t = [-1 1];
  for n = 1:2
    fint(:,1,n) = fint(:,1,n) + ...
        t(n).^2*0.25*(-f(:,1,1)+f(:,1,2));

    fint(:,1,n) = fint(:,1,n) + ...
        t(n)*0.5*(f(:,1,1)+f(:,1,2));

    fint(:,1,n) = fint(:,1,n) + ...
        (0.75*f(:,1,1) + 0.25*f(:,1,2));
  end
  % order 1 Gauss-Lobatto or equispaced

elseif ntime == 3
  t = [-1 0 1];
  for n = 1:3
    fint(:,1,n) = fint(:,1,n) + ...
        t(n).^3*(0.1666666666666667*(f(:,1,1)+f(:,1,3)) - ...
        0.3333333333333333*f(:,1,2));

    fint(:,1,n) = fint(:,1,n) + ...
        t(n).^2*(-0.25*(f(:,1,1)-f(:,1,3)));

    fint(:,1,n) = fint(:,1,n) + ...
        t(n)*f(:,1,2);

    fint(:,1,n) = fint(:,1,n) + ...
        0.4166666666666667*f(:,1,1) + ...
        0.6666666666666667*f(:,1,2) - ...
        0.0833333333333333*f(:,1,3);
  end
  % order 2 Gauss-Lobatto or equispaced
end



