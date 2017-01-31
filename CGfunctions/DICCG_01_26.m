
function[x,flag,relres,iter,resvec]=DICCG_01_26(a,b,z,tol,iter)
 xi=rand(size(a,1),1);
% b=b1(1:size(z,1),1);
h2=0;
h21=0;
size(z)
%size(a1)


%a=a1(1:size(z,1),1:size(z,1));
%full(a1)
size(a)
n=size(a,2);
if size(a) == size(z)
    Z1=z;
else
    Z1=eye(size(a));
    Z1(1:size(z,1),1:size(z,2))=z;
end
z=Z1;
e=z'*a*z;
q=z/e*z';
pd=eye(size(a*q))-a*q;
clear e
%b=b';
n=size(a,2);
r1=b-a*xi;
r1=pd*r1;
l=ichol(a);
r0=(l\r1);

p0=(l'\r0);
r0=r0;
fprintf('Dach')
[vdach,da]=eigs(inv(l*l')*pd*a,n);
conddach=condest(inv(l*l')*pd*a)
ab=abs(inv(l*l')*pd*a);
ab=sum(ab,1);
n1=max(ab);
ab=abs(inv(inv(l*l')*pd*a));
ab=sum(ab,1);
n1i=max(ab);
cond1=n1*n1i
condeff2=max(diag(real(da)))/min(diag(real(da)))
figure
 plot(diag(real(da)),'*')
 title('ddach');
for iter=1:iter
    w=pd*a*p0;
     alpha=(r0'*r0)/(w'*p0);    
     x=xi+alpha*p0;
     r=r0-alpha*(l\w); 
     beta=(r'*r)/(r0'*r0);
     p=(l'\r)+beta*p0;  
     p0=p;
     r0=r;
    if x==0
       e=0;
     else
          e=abs((x-xi)./x)*100;
     relres=sqrt(e'*a*e);
     resvect=relres;
      color=[0.1 0.5 0.5];
     figure(10)
     hline=plot(iter,relres,'*','Color',color);
     hold on
     end
       flag=0;
     if (relres>=tol)
         flag=1;
     end
     if flag==0
         break
     end     
     
    
     xi=x;
end
 x=q*b+pd'*xi;
 

