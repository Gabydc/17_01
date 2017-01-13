%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function
clear all
close all
clc
x=10;
y=10;
l=7;
s0=100;
s=10;



[a,b,z]=matrixf(x,y,l,s0,s);
 iteration=500;
 tol=10^-7;
 figure
spy(a)
%break
 xi(1:size(b,1))=rand;
 h=0;
for j=1:l
    
 for i=1:x*y
h=h+1;
     z1(h,j)=1;

 end
end
 spy(z1)
b=b';
load sol5_7
a1 = sol.A;
b1=sol.b;
L=ichol(a);
Z=z;
xi1(1:size(b1,1),1)=rand;
[x1,iter1,e1,hline1]=CGF(a,b,xi,iteration,tol);
[x11,fl0,rr0,it0,rv0] = pcg(a,b',tol,iteration,[],[],xi');
[x2,iter2,e2,hline2]=CGCh(a,b,xi,iteration,tol,L);
[x21,fl1,rr1,it1,rv1] = pcg(a,b',tol,iteration,L,L',xi');
[x3,iter3,e3,hline3]=DCGChF(a,b,xi,iteration,tol,z,L);
[x31,fl2,rr2,it2,rv2] = dpcg(a,b',Z,tol,iteration,L,L',xi');
ylabel('log(Error)')
xlabel('Iteration')
%legend([hline1,hline2,hline3,hline4],'CG','CGCH','DCG','DCGCh')
 fprintf('\n  Method      Iteration #    error   \n');
  fprintf('\n CG      %8d      %10.0d\n',iter1, max(e1));
   fprintf('\n CGm      %8d      %10.0d\n',it0, max(rr0));
  fprintf('\n CGCh    %8d      %10.0d\n',iter2, max(e2));
  fprintf('\n CGChm      %8d      %10.0d\n',it1, max(rr1));
  fprintf('\n DCG     %8d           %1.0d\n',iter3, max(e3));
  fprintf('\n DCGm      %8d      %10.0d\n',it2, max(rr2));
 

x5=a\b';
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
a1 = sol.A;
b1=sol.b;
p=sol.pressure;
h=0;
  for i=2:x*l-1
       for j=1:y
           h=h+1;
       T1(i,j)=x1(h);
        T11(i,j)=x11(h);
        T2(i,j)=x2(h);
        T21(i,j)=x21(h);
        T3(i,j)=x3(h);
        T31(j,j)=x31(h);
        T5(i,j)=x5(h);
    
       end 
  end


   figure
   subplot(1,2,1)
  mesh(T1)
   title('CG')
   subplot(2,2,2)
   mesh(T11)
   title('CGm')
   figure
   subplot(1,2,1)
  mesh(T2)
   title('ICCG')
   subplot(2,2,2)
   mesh(T21)
   title('ICCGm')
      figure
   subplot(1,2,1)
  mesh(T3)
   title('DICCG')
   subplot(2,2,2)
   mesh(T31)
   title('DICCGm')
   figure
   mesh(T5)
  title('a|b')
  figure
  semilogy(rv0./norm(b),'r.');
hold on
semilogy(rv1./norm(b),'b.');
semilogy(rv2./norm(b),'k.');
legend('No Preconditioner','ICCG','DICCG');
 