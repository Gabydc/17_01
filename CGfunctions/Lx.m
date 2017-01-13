function[y]=Lx(l,b) 
y(1)=b(1)/l(1,1);
n=size(l,2);
 for i=2:n
     s=0;
         for j=1:i-1
             s=s+l(i,j)*y(j);
         end
     y(i)=(b(i)-s)/l(i,i); 
 end
 end
 