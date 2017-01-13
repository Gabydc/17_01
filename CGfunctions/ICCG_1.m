
function[xf,iter,ee,hl1,h2,h21]=ICCG_1(a,b,xi,iteration,tol,l,eig,text,nn,dir,def,fc) 
n=size(a,2);
if eig==1
n=size(a,2);
fprintf('CG')
[V,D] = eigs(a,n);
 D=diag(D);
 D=abs(D);
 conda1=condest(a,2)
 max(D)
 min(D)
 conda2=max(D)/min(D)
 fprintf('ICCG')
[Va,Da] = eigs(inv(l)*a*inv(l'),n);
 Da=diag(Da);
 Da=abs(Da);
 condal1=condest(inv(l)*a*inv(l'),2)
 max(Da)
 min(Da)
 condadl2=max(Da)/min(Da)

 
 fileID = fopen(text,'a');
 fprintf(fileID,'%6s %6s %6.2e %6s\n','CG',' &',conda2, ' \\');
fprintf(fileID,'%6s %6s %6.2e %6s\n','ICCG',' &',condadl2, ' \\');

fclose(fileID);
 text1 = [dir def 'eigs_CGCh.txt'];
 fileID = fopen(text1,'w');
fprintf(fileID,'%6.2e \n', Da);
fclose(fileID);
file=[dir def 'eigs_CGCh'];
 save(file,'Da')

figure(30)
subplot(2,1,1)
h2=semilogy(D,'*');
axis('tight')
title('CG')
ylabel('log scale')
xlabel('Eigenvalue')
figure(30)
subplot(2,1,2)
h2=semilogy(Da,'*');
axis('tight')
title('ICCG')
ylabel('log scale')
xlabel('Eigenvalue')
hold on
figure(40)
subplot(2,1,1)
h21=semilogy(D,'*');
axis([1, n, min(D), max(D)])
title('CG')
ylabel('log scale')
xlabel('Eigenvalue')
figure(40)
subplot(2,1,2)
h21=semilogy(Da,'*');
axis([1, n, min(Da), max(Da)])
title('ICCG')
ylabel('log scale')
xlabel('Eigenvalue')

hold on
else h2=0; h21=0;
end
r0=b-a*xi;
lb=l\b;
r0=l\r0;
p0=l'\r0;
nor=abs(lb'*lb);
%nor=1;
for iter=1:iteration
    ap=a*p0; 
     alpha=(r0'*r0)/((ap)'*p0);   
     xf=xi+alpha*p0;
     r=r0-alpha*(l\ap);
     beta=(r'*r)/(r0'*r0);
     p=l'\r+beta*p0;  
     p0=p;
     r0=r;
      color=[0.2 0.8 0.6];
%     if xf==0
%        e=0;
%      else
%           e=(abs(xi-xf)./abs(xf))*100;
%      ee=sqrt(e'*a*e);
%      
%     
%      figure(1)
%      hl=semilogy(iter,(ee),'o','Color',color);
%      hold on
%     end
%     
%       fra=abs(a*xf-b);
%      b1=abs(b);
%       for rr=1:n
%           if b1(rr)==0
%               e1(rr)=0;
%           else             
%               e1(rr)=fra(rr)/b1(rr);
%           end
%       end
%      ee1=sqrt(e1*a*e1');
    %rr=inv(l)*r;
     figure(fc)
     ee=abs(r'*r)/nor;
     hl1=semilogy(iter,ee,'o','Color',color);
     hold on
     
%       if xf==0
%          e11=0;
%      else
%           e11=(sqrt((xi-xf)'*a*(xi-xf))./sqrt((xf')*a*(xf)))*100;
%      ee11=sqrt(e'*a*e);
%          figure(6)
%      hl11=semilogy(iter,(ee11),'o','Color',color);
%      hold on
%      end
     
       flag=0;
       
     if (ee>=tol)
         flag=1;
     end
     if flag==0
         break
     end     
     xi=xf;
 end
