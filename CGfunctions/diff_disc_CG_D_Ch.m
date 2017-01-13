%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function
clear all
close all
clc
x=5;
y=5;
l=7;
s0=100;
s=10;
n=x*y*l-2*x;
[a,b,z]=matrixf(x,y,l,s0,s);
 iteration=500;
 tol=10^-7;
%  figure
% 
% spy(z)

%break
 xi(1:size(b,1))=rand;
 fprintf('CG')
 conda=condest(a)
 [va,da]=eigs(a,n);
ab=abs(a);
 figure
 plot(diag(da),'o')
 title('CG')
r=sum(ab,1);
n1=max(r);
ab=abs(inv(a));
r=sum(ab,1);
n1i=max(r);
cond1=n1*n1i
condeff2=max(diag(real(da)))/min(diag(real(da)))
% m=va*inv(va'*a*va)*va';
% [va,da]=eigs(m,n);
%  figure
%  plot(diag(m),'o')
title('M')
ma=m*a;
[va,da]=eigs(ma,n);
 figure
 plot(diag(ma),'o')
min(diag(ma))
max(diag(ma))
 title('M-1a')
%  clear z
%  z(:,1)=va(:,1);
%  z(:,2)=va(:,2);
%  z(:,3)=va(:,3);
%  z(:,4)=va(:,4);
  
 figure
 plot(diag(da),'o')
 title('CG')
b=b';
l=ichol(a);
[x1,iter1,e1,hline1]=CGF(a,b,xi,iteration,tol);
[x2,iter2,e2,hline2]=CGCh(a,b,xi,iteration,tol,l);
[x3,iter3,e3,hline3]=DCGF(a,b,xi,iteration,tol,z);
[x4,iter4,e4,hline4]=DCGChF(a,b,xi,iteration,tol,z,l);
ylabel('log(Error)')
xlabel('Iteration')
legend([hline1,hline2,hline3,hline4],'CG','CGCH','DCG','DCGCh')
 fprintf('\n  Method      Iteration #    error   \n');
  fprintf('\n CG      %8d      %10.0d\n',iter1, max(e1));
  fprintf('\n CGCh %8d      %10.0d\n',iter2, max(e2));
  fprintf('\n DCG          %8d           %1.0d\n',iter3, max(e3));
  fprintf('\n DCGCh         %8d           %1.0d\n',iter4, max(e4));
x5=a\b';
l=7;
 h=0;
 T(1,1:y)=1;
T(l*y,1:y)=0;
 T1=T;
 T2=T;
 T3=T;
 T4=T;
 T5=T;

  % x5=a\b;
load sol5_7
a = sol.A;
b=sol.b;
p=sol.pressure;
h=0;
  for i=2:x*l-1
       for j=1:y
           h=h+1;
       T1(i,j)=x1(h);
        T2(i,j)=x2(h);
        T3(i,j)=x3(h);
        T4(i,j)=x4(h);
        T5(i,j)=x5(h);
      error(h)=(x1(h)-p(175-h))/p(175-h);
       end 
  end
  h=0;
for i=1:x*l
    for j=1:y
        h=h+1;
  T6(i,j)=p(175-h+1);
    end
end
   figure
   subplot(3,2,1)
  mesh(T1)
   title('CG')
   subplot(3,2,2)
   mesh(T2)
   title('CGCh')
   subplot(3,2,3)
   mesh(T3)
   title('DCG')
   subplot(3,2,4)
  mesh(T4)
    title('DCGCh')
   subplot(3,2,5)
   mesh(T5)
  title('a|b')
  
     
  
  figure
   subplot(2,1,1)
  mesh(T1)
   title('CG')
  subplot(2,1,2)
   mesh(T5)
  title('solver')
  figure
  mesh(abs((T1(:)-T6(:))/T6(:)))