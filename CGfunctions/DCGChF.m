
function[xf,iter,ee,hline]=DCGChF(a,b,xi,iteration,tol,z,l)
e=z'*a*z;
q=z/e*z';
pd=eye(size(a*q))-a*q;

b=b';
n=size(a,2);
r1=b-a*xi';
r1=pd*r1;

r0=(l\r1);

p0=(l'\r0)';
r0=r0';
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
for iter=1:iteration
    w=pd*a*p0';
     alpha=(r0*r0')/(w'*p0');    
     xf=xi+alpha*p0;
     r=r0-alpha*(l\w)'; 
     beta=(r*r')/(r0*r0');
     p=(l'\r')'+beta*p0;  
     p0=p;
     r0=r;
    if xf==0
       e=0;
     else
          e=abs((xf'-xi')./xf')*100;
     ee=sqrt(e'*a*e);
      color=[0.1 0.5 0.5];
     figure(10)
     hline=plot(iter,ee,'*','Color',color);
     hold on
     end
       flag=0;
     if (ee>=tol)
         flag=1;
     end
     if flag==0
         break
     end     
     
    
     xi=xf;
end
 xf=q*b+pd'*xi';
