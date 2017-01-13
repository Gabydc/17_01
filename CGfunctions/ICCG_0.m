function [xf,h,h2,h21]=ICCG_0(a,b,xi,iteration,tol,l,eig)
%a=A;
n=length(a);
if eig==1

fprintf('ICCG')
[Va,Da] = eigs(a,n);
 Da=diag(Da);
 Da=abs(Da);
 conda=condest(a,2)
 condadeff=max(Da)/min(Da)
 
 %fileID = fopen(text,'a');
%fprintf(fileID,'%6s %6s %6.2e %6s\n','CG',' &',condadeff, ' \\');
%fclose(fileID);
 %text1 = [dir 'eigs_CGCh.txt'];
 %fileID = fopen(text1,'w');
%fprintf(fileID,'%6.2e \n', Da);
%fclose(fileID);
%file=[dir  'eigs_CGCh'];
% save(file,'Da')
figure(80)
subplot(2,1,2)
h2=semilogy(Da,'*');
axis('tight')
title('a')
ylabel('log scale')
xlabel('Eigenvalue')
hold on
h21=0;
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

     figure(500)
     ee=abs(r'*r)/nor;
     h=semilogy(iter,ee,'o','Color',color);
     hold on
ylabel('log(||M^{-1}r^k||_2/||M^{-1}b||_2)','FontSize',16)
xlabel('Iteration','FontSize',16)

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
