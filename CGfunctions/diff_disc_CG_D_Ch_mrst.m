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

dir='';
def=1;
fc=15;

[a,b,z]=matrixf(x,y,l,s0,s);
 maxit=500;
 tol=10^-5;
 figure
spy(z)

%break
 xi(1:size(b,1))=rand;
b=b';
x5=a\b';
L=ichol(a);
xi1(1:size(b,1),1)=rand;
[x1,iter1,e1,hline1]=CGF(a,b,xi,maxit,tol);
[xc,fl0,rr0,it0,rv0] = pcg(a,b',tol,maxit,L,L',xi');
[x2,iter2,e2,hline2]=CGCh(a,b,xi,maxit,tol,L);
[x3,iter3,e3,hline3]=DICCG(a,b',xi',maxit,tol,z,L,0,'uno',l,dir,def,fc);
[xd,fl2,rr2,it2,rv2] = dpcg(a,b',z,tol,maxit,L,L',xi');
[x4,iter4,e4,hline4]=DCGChF(a,b,xi,maxit,tol,z,L);

ylabel('log(Error)')
xlabel('Iteration')
legend([hline1,hline2,hline3,hline4],'CG','CGCH','DCG','DCGCh')
 fprintf('\n  Method      Iteration #    error   \n');
  fprintf('\n CG      %8d      %10.0d\n',iter1, max(e1));
  fprintf('\n CGCh    %8d      %10.0d\n',iter2, max(e2));
  fprintf('\n DICCG     %8d           %1.0d\n',iter3, max(e3));
  fprintf('\n  DCGCh_F     %8d           %1.0d\n',iter4, max(e4));
  fprintf('\n pcg      %8d      %10.0d\n',it0, rr0);
  fprintf('\n  deflation1 %8d      %10.0d\n',it2, rr2);


 h=0;
 T(1,1:y)=1;
T(l*y,1:y)=0;
 T1=T;
 T2=T;
 T3=T;
 T4=T;
 T5=T;
   Tc=T;
  Td=T;


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
        T2(i,j)=x2(h);
        T3(i,j)=x3(h);
        T4(i,j)=x4(h);
        T5(i,j)=x5(h);
        Tc(i,j)=xc(h);
        Td(i,j)=xd(h);
      error(h)=(x1(h)-p(175-h))/p(175-h);
       end 
  end
  Te=T5-Td;
  Te1=T5-T3;
  Te2=Td-T3;
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
   title('DICCG')
   subplot(3,2,4)
  mesh(T4)
    title('DCGCh')
   subplot(3,2,5)
   mesh(Tc)
  title('CGCh_f')
   subplot(3,2,6)
   mesh(Td)
  title('DCGCh_f')
  
     
  
  figure
   subplot(3,1,1)
  mesh(Te)
   title('errordefmat')
  subplot(3,1,2)
   mesh(Te1)
  title('erordefown')
    subplot(3,1,3)
   mesh(Te2)
  title('erordefmatvsdefown')
%   figure
%   mesh(abs((T1(:)-T6(:))/T6(:)))