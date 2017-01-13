
function[xf,iter,ee,hline]=CGCh(a,b,xi,iteration,tol,l) 
b=b';
n=size(a,2);
r1=b-a*xi';

r0=l\r1;
p0=l'\r0;
r0=r0';
p0=p0';
condch=condest(inv(l*l')*a)
[vach,dach]=eigs(inv(l*l')*a)

for iter=1:iteration
     alpha=(r0*r0')/((a*p0')'*p0') ;   
     xf=xi+alpha*p0;
     r=r0-alpha*((inv(l)*a)'*p0')';
     beta=(r*r')/(r0*r0');
     p=(l'\r')'+beta*p0 ;
     p0=p;
     r0=r;
    if xf==0
       e=0;
     else
          e=abs((xf'-xi')./xf')*100;
     ee=sqrt(e'*a*e);
     
     color=[0.2 0.8 0.6];
     figure(1)
     hline=plot(iter,log(ee),'o','Color',color);
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