
function[xf,iter,ee,hline]=DCGF1(a,b,xi,iteration,tol,z) 
e=z'*a*z;
q=z/e*z';
pd=eye(size(a*q))-a*q;
b=b';
xi=xi';
n=size(a,2);
r=-(a*xi-b);
r=pd*r;
p=r;
fprintf('Da')
[vda,da]=eigs(pd*a,n);
condad=condest(pd*a)
ab=abs(pd*a);
ab=sum(ab,1);
n1=max(ab);
ab=abs(inv(pd*a));
ab=sum(ab,1);
n1i=max(ab);
cond1=n1*n1i
condeff2=max(diag(real(da)))/min(diag(real(da)))
figure
 plot(diag(real(da)),'p')
 title('dda');
for iter=1:iteration

    alpha=r'*r/((pd*a*p)'*p);
    xf=xi+alpha*p;
    r1=r-alpha*pd*a*p;
    beta=(r1'*r1)/(r'*r);
    p1=r1+beta*p;
    
     if xf==0
         e=0;
     else
          e=abs((xi'-xf')./xf')*100;
     ee=sqrt(e*a*e'); 
     figure(10)
     color=[0.9 0.8 0.2];
     hline=plot(iter,log(ee),'s','Color',color);
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
         r=r1;
         p=p1;
         
    
end

xf=q*b+pd'*xf;


