
% a= [  4 3 1
%       3 5 2
%       1 2 6];
%  b=[3125 3650 2800]; 

function[l]=icholeskyf(a) 
l(1,1)=sqrt(a(1,1));
n=size(a,2);
 for k=2:n
     for i=1:k-1
 if a(k,i)~=0
         s=0;
         if i~=1
         for j=1:i-1
             s=s+l(i,j)*l(k,j);
         end
         end
         l(k,i)=(a(k,i)-s)/l(i,i);
     end
     l(k,k)=a(k,k);
     s=0;
     for j=1:k-1
             s=s+l(k,j)*l(k,j);
     end
         l(k,k)=sqrt(a(k,k)-s);
     end
 end
 
 