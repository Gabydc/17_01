function[x]=Ux(l,y) 

n=size(l,2);
x(n)=y(n)/l(n,n);
 for i=n-1:-1:1
     s=0;
         for j=i+1:n
             s=s+l(j,i)*x(j);
         end
     x(i)=(y(i)-s)/l(i,i); 
 end
 end
 