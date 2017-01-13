close all
clear all
clc
%Create a symmetric positive definite matrix, A.

N = 100;
A = delsq(numgrid('S',N));
%Create an incomplete Cholesky factorization as a preconditioner for pcg. Use a constant vector as the right hand side. As a baseline, execute pcg without a preconditioner.

b = ones(size(A,1),1);
tol = 1e-8;
maxit = 100;
Z=diag(A);
[x0,fl0,rr0,it0,rv0] = cgs(A,b,tol,maxit);
%Note that fl0 = 1 indicating that pcg did not drive the relative residual to the requested tolerance in the maximum allowed iterations. Try the no-fill incomplete Cholesky factorization as a preconditioner.
L1 = ichol(A);
[x1,fl1,rr1,it1,rv1] = cgs(A,b,tol,maxit,L1,L1');
[x2,fl2,rr2,it2,rv2] = dcgs(A,b,Z,tol,maxit,L1,L1');
semilogy(rv0./norm(b),'r.');
hold on
semilogy(rv1./norm(b),'b.');
semilogy(rv2./norm(b),'k.');
legend('No Preconditioner','ICCG','DICCG');