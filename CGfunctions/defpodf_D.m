%This program compute the dedaltion vectors from POD from D=x'*x

function [U,S,hd]=defpodf_D(x,dir)
x=normc(x);
lx=size(x,1);
ly=size(x,2);
D=x'*x;
[V,L] = eigs(D,ly);
S=zeros(lx,ly);
g=diag(L);
g1=sqrt(g);
g1=1./g1;
S(1:ly,1:ly)=diag(g1);
S=sparse(S);
V=sparse(V);
x=sparse(x);
U=x*V*S';
%[V,D] = eigs(X,n);
figure(3000)
hd=semilogy(abs(g),'*b'); 
%title('Eigenvalues X');
ylabel('log(Value) ','FontSize',16)
xlabel('Eigenvalue','FontSize',16)
axis('tight')
set(gca,'FontSize',13,'Xdir','reverse')
file='eig_pod';
B=[dir  file '.fig'];
saveas(hd,B)
B=[dir  file '.jpg'];
saveas(hd,B)
