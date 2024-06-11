
function Xrot = rotateVes(X,center,theta)

 Xrot = zeros(size(X));
 x = X(1:end/2); y = X(end/2+1:end);
 Xrot(1:end/2) = (x-center(1))*cos(theta) - (y-center(2))*sin(theta) + center(1);
 Xrot(end/2+1:end) = (x-center(1))*sin(theta) + (y-center(2))*cos(theta) + center(2);

end