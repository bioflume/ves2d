classdef householder 
% simple recursive implementation of householder qr vactorization
% 
% NOTES
%   * the routine is slow
%   * the routine is not optimized for round off
%   * for robust factorizations, use MATLAB's qr() 
%    
%
% usage example:
%   
%  m = 4; n = 2;   
%  A = rand(m,n);
%  o = householder;
%  [Q,R]=o.qr(A);
%   
% see the documentation of the individual methods

methods
    

%/******************************************************/
  function u = reflect(this, w, v)
% function u = reflect(w, v)
% create householder reflection u so that w gets aligned to v
% the norm of v should be 1.
    
    wn = norm(w);
    if norm(w)<eps, u=0*w; return; end;
    
    u = w - wn*v;
    t = norm(u);
    % check for zero u, for example if v is parallel to w.
    if t<eps(wn), u=0*w; return; end
    u = u/t;
  end


%/******************************************************/
  function u = reflect_e1(this, w)
% function u = reflect_e1(w)
% create householder reflection u so that w gets aligned to e1
    
    wn = norm(w);
    if norm(w)<eps, u=0*w; return; end;
    
    u=w;
    sigma = sign(w(1))*wn;
    u(1) = w(1) + sigma;
    un = sqrt(2*sigma*u(1));
    if un<eps(wn), u=0*w; return; end
    u=u/un;
  end




%/******************************************************/
  function B=transform(this, A, u)
% function B=transform(A, u)      
% apply Householder reflection 1-2uu' on each column of the matrix A
    if norm(u) == 0, B=A; return; end;

    vt = u'*A;
    B = A - 2*(u*vt);
   end
%/******************************************************/
  function [R,U] = getUonly(this,A)
    [m,n] = size(A); assert(m>=n);

    [R,U] = this.rqr(A);
    R = R(1:n,1:n);
  end
%/******************************************************/
  function [Q,R,U]=qr(this,A)
% function [Q,R,U]=qr(A)      
% compute the QR factorization of a matrix A in C^(m x n), where
% m >= n. U holds all the Householder vectors.

    % initializations and checks
    [m,n] = size(A); assert(m>=n);
    
    % call the recursive version of the method
    [R,U]=this.rqr(A);

    % build the Q matrix. 
    Q=this.getQ(U);

    % kill the zero blocks (MATLAB's qr(A,0));
    R=R(1:n,1:n);
  end

%/******************************************************/
% function U=nrqr(A)      
% compute the QR factorization of a matrix A in C^(m x n), where
% m >= n. U holds all the Householder vectors. This is the non-recursive qr
function U=nrqr(this,A)
    [m,n] = size(A);
    U = A*0;
    for i=1:n
        u = this.reflect_e1(A(i:m,i));
        A(i:m,i:n) = this.transform(A(i:m,i:n),u);
        U(i:m,i) = u;
    end
end
  
%/******************************************************/
% function [R,U]=rqr(A)      
% compute the QR factorization of a matrix A in C^(m x n), where
% m >= n. U holds all the Householder vectors. This is the recurcive qr
function [R,U]=rqr(this,A)
    [m,n] = size(A);
    U = A*0; 
    u = this.reflect_e1(A(:,1));
    R = this.transform(A,u);
    U(:,1) = u;
    if n==1, return; end;  % base case
    [Rr, Ur] = this.rqr( R(2:m,2:n) );
    R(2:m,2:n) = Rr;
    U(2:m,2:n) = Ur;
  end

%/******************************************************/
function [v] = Q(this,U,w)
% function [v] = Q(U,w)
% Given the Householder vectors U, apply them to a vector w
% so that v = Q*w;
%
  [m,n] = size(U);
  v=w;
  for k=n:-1:1
      in = v(k:end);
      u  = U(k:end,k);
      out = this.transform(in, u);
      v(k:end)=out;
  end
end

%/******************************************************/
function Q=getQ(this, U)
%function Q=getQ(U)
% Given the Householder vectors U (obtain by the qr method), construct
% the Q matrix.

  [m,n] = size(U);
  Im = eye(m);
  Q=0*U;
  for k=1:n
     Q(:,k) = this.Q(U,Im(:,k));
  end
end

%/******************************************************/
function QN=getQN(this, U)
%function Q=getQ(U)
% Given the Householder vectors U (obtain by the qr method), construct
% the Q matrix.

  [m,n] = size(U);
  Im = eye(m);
  QN=0*U;
  idx = 1;
  for k=n+1:m
     QN(:,idx) = this.Q(U,Im(:,k));
     idx = idx+1;
  end
end

%/******************************************************/
function A=applypreQNt(this,A,U)
%function Q=getQ(U)
% Given the Householder vectors U (obtain by the qr method), construct
% the Q matrix.

    [m,n] = size(U);
    for i=1:n
        u=U(:,i);
        vt = u'*A;
        A = A - 2*(u*vt);
    end
end

%/******************************************************/
function A=applypreQN(this,A,U)
%function Q=getQ(U)
% Given the Householder vectors U (obtain by the qr method), construct
% the Q matrix.

    [m,n] = size(U);
    for i=n:-1:1
        u=U(:,i);
        vt = u'*A;
        A = A - 2*(u*vt);
    end
end

%/******************************************************/
function A=applypostQN(this,A,U)
%function Q=getQ(U)
% Given the Householder vectors U (obtain by the qr method), construct
% the Q matrix.

    [m,n] = size(U);
    for i=1:n
        u=U(:,i);
        vt = A*u;
        A = A - 2*(vt*u');
    end
end

end % methods

end % class
      
  
  

   
