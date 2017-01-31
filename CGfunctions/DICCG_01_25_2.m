
function[x,flag,relres,iter,resvec]=DICCG_01_25_1(a1,b1,z,tol,iter)
xi=zeros(size(z,1),1);
b=b1(1:size(z,1),1);
h2=0;
h21=0;
size(z)
size(a1)


a=a1(1:size(z,1),1:size(z,1));

%  h2=0;
% h21=0;
% if size(a) == size(z)
%     Z1=z;
% else
%     Z1=eye(size(a));
%     Z1(1:size(z,1),1:size(z,2))=z;
% end
%   z=Z1;
  eig=1;
n=size(a,2);
e=z'*a*z;
ei=sparse(inv(e));
l=ichol(a);
xi=zeros(size(b));
%size(ez)
%  file='E';
%    B=[dir  file def '.fig'];
%   saveas(hd(1),B)
%    B=[dir  file def '.jpg'];
%   saveas(hd(1),B)
%  fileID = fopen(text,'a');
% fprintf(fileID,'%6s %6s %6.2e %6s\n','DCGCh',' &',condadeff, ' \\');
% fclose(fileID);
if eig==1
    if n<3000
     ez=ei*z';
q=z*ez;
pd=eye(n)-a*q;
  pda=pd*a;   
fprintf('DICCG')
[Va,Da] = eigs(inv(l)*pda*inv(l'),n);
 Da=diag(Da);
 Da=abs(Da);
 conda=condest(inv(l)*pda*inv(l'),2)
 condadeff=max(Da)/min(Da)
 

figure(30)
subplot(2,1,2)
h2=semilogy(Da,'*');
axis('tight')
title('DICCG')
ylabel('log scale)')
xlabel('Eigenvalue')
hold on
figure(40)
subplot(2,1,2)
Da=Da;
h21=semilogy(Da,'*');
axis([1, n, min(Da), max(Da)])
title('DICCG')
ylabel('log scale')
xlabel('Eigenvalue')
hold on
    end
else h2=0; h21=0;
end

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
      resvec=relres;
     figure(1235)
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
b1(size(z,1)+1)
a1(size(z,1)+1,size(z,1)+1)
x(1,size(z,1)+1)=b1(size(z,1)+1)/a1(size(z,1)+1,size(z,1)+1);
x(1,size(z,1)+2)=b1(size(z,1)+2)/a1(size(z,1)+2,size(z,1)+2);