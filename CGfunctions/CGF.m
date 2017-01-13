
function[xf,iter,ee,hline]=CGF(a,b,xi,iteration,tol) 
b=b';
xi=xi';
n=size(a,2);
r=-(a*xi-b);
p=r;


for iter=1:iteration
    alpha=r'*r/(p'*a*p);
    xf=xi+alpha*p;
    r1=r-alpha*a*p;
    beta=(r1'*r1)/(r'*r);
    p1=r1+beta*p;
    
     if xf==0
         e=0;
     else
          e=abs((xi'-xf')./xf')*100;
     ee=sqrt(e*a*e');
     color=[0.7 0.3 0.6];
         figure(10)
     hline=plot(iter,log(ee),'*','Color',color);
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
