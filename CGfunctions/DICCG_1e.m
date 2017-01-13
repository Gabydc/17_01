
function[xf,iter,ee,hl1,h2,h21]=DICCG(a,b,xi,iteration,tol,z,l,eig,text,nn,dir,def,fc)
h2=0;
h21=0;
n=size(a,2);
e=z'*a*z
ei=sparse(inv(e));
%size(ez)
file=[dir def 'e'];
  save(file,'e')
  de=det(e)
  file=[dir def 'de'];
  save(file,'de')
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
 
 fileID = fopen(text,'a');
fprintf(fileID,'%6s %6s %6.2e %6s\n','DCGCh',' &',condadeff, ' \\');
fclose(fileID);
 text1 = [dir def 'eigs_DCGCh_1.txt'];
 fileID = fopen(text1,'w');
fprintf(fileID,'%6.2e \n', Da);
fclose(fileID);
file=[dir def 'eigs_DCGCh_1'];
 save(file,'Da')
figure(80)

subplot(2,1,1)
h2=semilogy(Da,'*');
axis('tight')
title('DICCG_1')
ylabel('log scale)')
xlabel('Eigenvalue')
hold on

[Ve,De] = eigs(e,nn);
 De=diag(De);
 De=abs(De);
  conde=max(De)/min(De)
figure(90)
h21=semilogy(De,'p');
axis('tight')
title('eigenvalores E')
ylabel('log scale')
xlabel('Eigenvalue')
hold on
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
for iter=1:iteration
         [ap]=deflatevect(z,ei,a,a*p0);
         [apt]=tdeflatevect(z,ei,a,p0);
     alpha=(r0'*r0)/((ap)'*p0);   
     xf=xi+alpha*p0;
     r=r0-alpha*(l\ap);
     beta=(r'*r)/(r0'*r0);
     p=l'\r+beta*p0;  
     p0=p;
     r0=r;
     color=[0.7 0.5 0.7];

      ee=abs(r'*r)/nor;
     figure(fc)
     hl1=semilogy(iter,ee,'p','Color',color);
     hold on
       flag=0;

     if (ee>=tol)
         flag=1;
     end
     if flag==0
         break
     end     
     xi=xf;
end
[xf]=tdeflatevect(z,ei,a,xf);
qb=z'*b;
qb=ei*qb;
qb=z*qb;
xf=qb+xf;