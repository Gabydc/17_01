
function[x,flag,relres,iter,resvec]=DICCG_01_25_1(a,b,z,tol,iter)
 xi=zeros(size(a,1),1);
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
    Z1=z;
    for j=(size(z,1)+1):size(a,1)
        rt=j-size(z,1);
        
    for i=1:(size(a,1)-size(z,2))
        if rt==i
    Z1(j,i)=1;
        end
    end
    end
end
full(a)
z=eye(size(a));
full(Z1)
  z=Z1;
e=z'*a*z;
eigs(e);
ei=(inv(e));
eigs(ei)
l=ichol(a);
[pb]=deflatevect(z,ei,a,b);
lb=l\b;
plb=l\pb;
nor1=abs(lb'*lb);
r0=b-a*xi;
[r0]=deflatevect(z,ei,a,r0);
r0=l\r0;
p0=l'\r0;
nor=abs(lb'*lb);
%nor=1;
for iter=1:iter
         [ap]=deflatevect(z,ei,a,a*p0);
         [apt]=tdeflatevect(z,ei,a,p0);
     alpha=(r0'*r0)/((ap)'*p0);   
     x=xi+alpha*p0;
     r=r0-alpha*(l\ap);
     beta=(r'*r)/(r0'*r0);
     p=l'\r+beta*p0;  
     p0=p;
     r0=r;
     color=[0.1 0.5 0.5];

      relres=abs(r'*r)/nor;
      resvec(:,iter)=relres;
     figure(100)
     hl1=semilogy(iter,relres,'p','Color',color);
     hold on
       flag=0;
     
     if (relres>=tol)
         flag=1;
     end
     if flag==0
         break
     end     
     xi=x;
end

[x]=tdeflatevect(z,ei,a,x);

qb=z'*b;
qb=ei*qb;
qb=z*qb;
x=qb+x;
% b1(size(z,1)+1)
% a1(size(z,1)+1,size(z,1)+1)
% x(1,size(z,1)+1)=b1(size(z,1)+1)/a1(size(z,1)+1,size(z,1)+1);
% x(1,size(z,1)+2)=b1(size(z,1)+2)/a1(size(z,1)+2,size(z,1)+2);