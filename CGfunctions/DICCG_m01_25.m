
function[x,flag,relres,ii,resvec]=DICCG_m2(A,b,Z,tol,maxit,M1,M2,x0,varargin)
%DPCG  Deflated Preconditioned Conjugate Gradients Method.
%   X = DPCG(A,B,Z) attempts to solve the system of linear equations PA*X=PB for
%   X. The N-by-N coefficient matrix A must be symmetric and positive
%   definite, the right hand side column vector B must have length N, the 
%   matrix X must have N rows.
%
%   X = DPCG(AFUN,B) accepts a function handle AFUN instead of the matrix A.
%   AFUN(X) accepts a vector input X and returns the matrix-vector product
%   A*X. In all of the following syntaxes, you can replace A by AFUN.
%
%   X = DPCG(A,B,Z,TOL) specifies the tolerance of the method. If TOL is []
%   then PCG uses the default, 1e-6.
%
%   X = DPCG(A,B,Z,TOL,MAXIT) specifies the maximum number of iterations. If
%   MAXIT is [] then PCG uses the default, min(N,20).
%
%   X = DPCG(A,B,Z,TOL,MAXIT,M) and X = PCG(A,B,TOL,MAXIT,M1,M2) use symmetric
%   positive definite preconditioner M or M=M1*M2 and effectively solve the
%   system inv(M)*A*X = inv(M)*B for X. If M is [] then a preconditioner
%   is not applied. M may be a function handle MFUN returning M\X.
%
%   X = DPCG(A,B,Z,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0 is
%   [] then PCG uses the default, an all zero vector.
%
%   [X,FLAG] = DPCG(A,B,...) also returns a convergence FLAG:
%    0 DPCG converged to the desired tolerance TOL within MAXIT iterations
%    1 DPCG iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 DPCG stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during PCG became too
%      small or too large to continue computing.
%
%   [X,FLAG,RELRES] = DPCG(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = DPCG(A,B,...) also returns the iteration number
%   at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = DPCG(A,B,...) also returns a vector of the
%   estimated residual norms at each iteration including NORM(B-A*X0).
%
%   Example:
%      n1 = 21; A = gallery('moler',n1);  b1 = A*ones(n1,1);
%      tol = 1e-6;  maxit = 15;  M = diag([10:-1:1 1 1:10]);
%      [x1,flag1,rr1,iter1,rv1] = dpcg(A,b1,tol,maxit,M);
%   Or use this parameterized matrix-vector product function:
%      afun = @(x,n)gallery('moler',n)*x;
%      n2 = 21; b2 = afun(ones(n2,1),n2);
%      [x2,flag2,rr2,iter2,rv2] = dpcg(@(x)afun(x,n2),b2,tol,maxit,M);
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTAB, BICGSTABL, CGS, GMRES, LSQR, MINRES, QMR,
%   SYMMLQ, TFQMR, ICHOL, FUNCTION_HANDLE.

%   Copyright 1984-2013 The MathWorks, Inc.

if (nargin < 3)
    error(message('MATLAB:dpcg:NotEnoughInputs'));
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    if (m ~= n)
        error(message('MATLAB:dpcg:NonSquareMatrix'));
    end
    if ~isequal(size(b),[m,1])
        error(message('MATLAB:dpcg:RSHsizeMatchCoeffMatrix', m));
    end
%     if ~isequal(size(Z,1),m)
%         error(message('MATLAB:dpcg:ZsizeMatchCoeffMatrix', m));
%     end
else
    m = size(b,1);
    n = m;
    if ~iscolumn(b)
        error(message('MATLAB:dpcg:RSHnotColumn'));
    end
end

% Assign default values to unspecified parameters
if (nargin < 4) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol <= eps
    warning(message('MATLAB:dpcg:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:dpcg:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if (nargin < 5) || isempty(maxit)
    maxit = min(n,20);
end
size(A)
size(Z)
% if size(A) == size(Z)
%     Z1=Z;
% else
%     Z1=eye(size(A));
%     Z1(1:size(Z,1),1:size(Z,2))=Z;
% end
%  Z=Z1;
 if size(A) == size(Z)
     A1=A;
     b1=b;
     
 else
   %  A1=eye(size(Z));
     A1=A(1:size(Z,1),1:size(Z,1));
     b1=b(1:size(Z,1));
 end
  A=A1;


E = Z' * A * Z;
EI = sparse(inv(E));
% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                     % Norm of rhs vector, b
%[Pb]=dvect(Z,EI,A,b);
%n2Pb = norm(Pb);                     % Norm of rhs vector, b
if (n2b == 0)                      % if    rhs vector is all zeros
    x = zeros(n,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    ii = 0;                      % no iterations need be performed
    resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 3)
        itermsg('dpcg',tol,maxit,0,flag,ii,NaN);
    end
    return
end

if ((nargin >= 6) && ~isempty(M1))
    existM1 = 1;
    [m1type,m1fun,m1fcnstr] = iterchk(M1);
    if strcmp(m1type,'matrix')
        if ~isequal(size(M1),[m,m])
            error(message('MATLAB:dpcg:WrongPrecondSize', m));
        end
    end
else
    existM1 = 0;
    m1type = 'matrix';
end

if ((nargin >= 7) && ~isempty(M2))
    existM2 = 1;
    [m2type,m2fun,m2fcnstr] = iterchk(M2);
    if strcmp(m2type,'matrix')
        if ~isequal(size(M2),[m,m])
            error(message('MATLAB:dpcg:WrongPrecondSize', m));
        end
    end
else
    existM2 = 0;
    m2type = 'matrix';
end

if ((nargin >= 8) && ~isempty(x0))
    if ~isequal(size(x0),[n,1])
        error(message('MATLAB:dpcg:WrongInitGuessSize', n));
    else
        x = x0;
    end
else
    x = zeros(n,1);
end

if ((nargin > 8) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error(message('MATLAB:dpcg:TooManyInputs'));
end

[Pb]=dvect(Z,EI,A,b);
% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
n2b=norm(b);                % Relative tolerance
tolb = tol * n2b;                  % Relative tolerance
n2Pb = norm(Pb);
tolPb = tol * n2Pb;  
lb = iterapp('mldivide',m1fun,m1type,m1fcnstr,b,varargin{:});
plb = iterapp('mldivide',m1fun,m1type,m1fcnstr,Pb,varargin{:});
n2lb = norm(lb);
n2Plb = norm(plb);
r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
normr = norm(r);                   % Norm of residual
normr_act = normr;
x0=x(1:size(Z,1));
if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    ii = 0;
    resvec = normr;
    
    if (nargout < 2)
        itermsg('dpcg',tol,maxit,0,flag,ii,relres);
    end
    return
end
%Deflated residual
[r]=dvect(Z,EI,A,r);

normr = norm(r);                   % Norm of residual
normr_act = normr;


resvec = zeros(maxit+1,1);         % Preallocate vector for norm of residuals
resvec(1,:) = normr;               % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of minimum residual
rho = 1;
stag = 0;                          % stagnation of the method
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;



%lb=M1\b;
%plb=M1\Pb;

%r0=b-A*x0;
%[r0]=dvect(Z,EI,A,r0);
%r0=M1\r0;
%p0=M2\r0;
nor=abs(lb'*lb);


     if existM1
        r0 = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
        if ~all(isfinite(r0))
            flag = 2;
            return
        end
    else % no preconditioner
        r0 = r;
    end
%nor=1;
for ii=1:maxit
         if existM2
         p0 = iterapp('mldivide',m2fun,m2type,m2fcnstr,r0,varargin{:});
         if ~all(isfinite(p0))
             flag = 2;
            break
         end
     else % no preconditioner
         p0 = r0;
         end
         
    rho1 = rho;
    rho = r0' * r0;
    if ((rho == 0) || isinf(rho))
        flag = 4;
        break
    end
    if (ii == 1)
        p = p0;
    else
        beta = rho / rho1;
        if ((beta == 0) || isinf(beta))
            flag = 4;
            break
        end
        p = p0 + beta * p;
    end
     q = iterapp('mtimes',afun,atype,afcnstr,p,varargin{:});
         [ap]=dvect(Z,EI,A,q);
         [apt]=tdvect(Z,EI,A,q);
         pq = p' * q;
     if ((pq <= 0) || isinf(pq))
        flag = 4;
        break
    else
        %alpha = (r0'*r0) / pq;
        alpha = rho/ pq;
    end    
     if isinf(alpha)
        flag = 4;
        break
     end
    
     if existM1
         r1 = iterapp('mldivide',m1fun,m1type,m1fcnstr,ap,varargin{:});
         if ~all(isfinite(r1))
             flag = 2;
             break
         end
     else % no preconditioner
         r1 = ap;
     end
     x=x0+alpha*p;
     %rl = iterapp('mldivide',m1fun,m1type,m1fcnstr,ap,varargin{:});
     r=r0-alpha*(r1);
    
     
     
     %alpha=(r0'*r0)/((ap)'*p0);   
    % x=x0+alpha*p0;
     %r=r0-alpha*(M1\ap);
     %beta=(r'*r)/(r0'*r0);
     %p=M2\r+beta*p0;  
     p0=p;
     r0=r;

     color=[0.1 0.5 0.5];

 
     %problem is when computing norm(r) is not the same as abs(r'*r)
      normr = norm(r)
    normr_act = normr;
    resvec(ii+1,1) = normr;
     

      relres=abs(r'*r)/n2lb
      resvec(1,ii)=norm(r);
     figure(100)
     hl1=semilogy(ii,relres,'s','Color',color);
     hold on
       flag=0;

     if (relres>=tol)
         flag=1;
     end
     if flag==0
         break
     end     
     x0=x;
end
[x]=tdvect(Z,EI,A,x);
[Qb]=Qm(Z,EI,b);
x=Qb+x;

end
function[Px]=dvect(Z,EI,A,x)
[Qx]=Qm(Z,EI,x);
Px=x-A*Qx;
end
function[Px]=tdvect(Z,EI,A,x)
Px=A'*x;
Px=Qm(Z,EI,Px);
Px=x-Px;
end
function[Qx]=Qm(Z,EI,x)
Qx=Z'*x;
Qx=EI*Qx;
Qx=Z*Qx;
end

