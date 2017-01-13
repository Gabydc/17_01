%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function
function solve_DICCG_f_e(num,nn,e,iteration,dir,k)
%This is the matrix that we are going to use, usually is A12
tol=10^-k;
fc=100;
%num='5';  

%The number of deflation verctors to use either 2 or 3
% nn=5;
% e=0;
% pw=0;
% size=32*2^pw;
%  iteration=500;
% k=5;
%  
  
%  dir='/home/wagm/cortes/Localdisk/proofs/';
%  folder=[problem num2str(size)  ];
%  dir = [dir folder '/'];
%  nxi=1;
%  nyi=1;
%  nx=60;
%  ny=220;
%   name=['accuracy_' num2str(nx) '_' num2str(ny) ];
%   mkdir([dir], name )
%  dir=[dir name '/'];

def=['tol-' num2str(k) '_' num2str(nn) ];
text = [dir 'iter' def '.txt'];
text1 = [dir 'cond' def '.txt'];
file=[dir 'xi'];
load(file)
file=[dir 'A' num];
load(file)
file=[dir 'b' num];
load(file)
file=[dir 'p' num];
load(file)
file=[dir 'W' num];
load(file)
file=[dir 'G'];
load(file)
file=[dir 'l'];
load(file)
% for i=1:nn
% file=[dir 'x' num2str(i) ];
% x=load(file);
% 
% end
file=[dir 'x'];
load(file)
file=[dir 'x1'];
load(file)
n=length(b);

%[V,D]=defpodf_w_2(x,nn);

a=A;
p=A\b;
figure(1000)
ht=plotingsolution_1(G,W,'backslash',p);
view(0,90)
%clear p
clear A
%%1
%   [V,D]=POD(x);
%  z(:,1)=x(:,5);
    s=0;
     for i=1:nn
        s=s+1;
     z(:,s)=x(:,i);
      z1(:,s)=x1(:,i);
     
% % % %   %z(:,s)=V(:,s);
    end


fileID = fopen(text1,'w');
fprintf(fileID,'\n%6s %6s %6s %6s\n\n','Method ','& ', 'Cond',' \\');
fclose(fileID);
%[x1,iter1,e1,hl1,hl11,hl111,h2,h21]=CGF(a,b,xi,iteration,tol,e,text,nn,dir,def); 

[x2,iter2,e2,hl2,h2,h21]=ICCG(a,b,xi,iteration,tol,l,e,text1,nn,dir,def,fc);


%[x3,iter3,e3,hl3,hl13,hl113,h2,h21]=DCGF(a,b,xi,iteration,tol,z,e,text1,nn,dir,def);

[x4,iter4,e4,hl4,h2,h21]=DICCG_e(a,b,xi,iteration,tol,z,l,1,text1,nn,dir,def,fc);

[x5,iter5,e5,hl5,h25,h215]=DICCG_1e(a,b,xi,iteration,tol,z1,l,1,text1,nn,dir,def,fc);


 % x1(:,1)=x1;
  x2(:,1)=x2;
  %x3(:,1)=x3;
  x4(:,1)=x4;
  x5(:,1)=x5;
% % filename=[dir '/results/qfs_1/xch' num];
% % save(filename,'xch')
% figure(1)
% legend([hl1,hl2,hl3,hl4],'CG','PCG','DCG','DPCG')
% ylabel('log(Error)')
% xlabel('Iteration')
% title('||x_i-x_f||_2/||x_f||_2')
figure(fc)
legend([hl2,hl4,hl5],'ICCG','DICCG','DICCG_1')
ylabel('log(||M^{-1}r^k||_2/||M^{-1}b||_2)','FontSize',16)
xlabel('Iteration','FontSize',16)


% figure(6)
% legend([hl111,hl112,hl113,hl114],'CG','PCG','DCG','DPCG')
% ylabel('log(Error)')
% xlabel('Iteration')
% title('||x_i-x_f||_A/||x_f||_A')
%txt file with data
fileID = fopen(text,'w');
fprintf(fileID,'\n%6s %6s %6s %6s %6s\n\n','Method &','Iter ', '&','error', ' \\');
%fprintf(fileID,'%6s %d %6s %6.2e %6s\n','CG &',iter1,' &',max(e1), ' \\');
fprintf(fileID,'%6s %d %6s %6.2e %6s\n','ICCG &',iter2,' &',max(e2), ' \\');
%fprintf(fileID,'%6s %d %6s %6.2e %6s\n','DCG &',iter3,' &',max(e3), ' \\');
fprintf(fileID,'%6s %d %6s %6.2e %6s\n','DICCG &',iter4,' &',max(e4), ' \\');
fprintf(fileID,'%6s %d %6s %6.2e %6s\n','DICCG1 &',iter5,' &',max(e5), ' \\');
fclose(fileID);
 fprintf('\n  Method      Iteration #    error   \n');
 % fprintf('\n CG      %8d      %10.0d\n',iter1, max(e1));
  fprintf('\n ICCG %8d      %10.0d\n',iter2, max(e2));
  %fprintf('\n DCG          %8d           %1.0d\n',iter3, max(e3));
  fprintf('\n DICCG         %8d           %1.0d\n',iter4, max(e4));
  fprintf('\n DICCG1         %8d           %1.0d\n',iter5, max(e5));

% 
 figure(fc+1)
%h1=plotingsolution(G,W,'CG',x1,1);
h1=plotingsolution(G,W,'BS',p,1);
colorbar

h1=plotingsolution(G,W,'ICCG',x2,2);
colorbar
%h1=plotingsolution(G,W,'DCG',x3,3);

figure(fc+2)
h1d=plotingsolution(G,W,'DICCG',x4,1);
colorbar
h1d=plotingsolution(G,W,'DICCG1',x5,2);
colorbar
figure(fc+3)
%h1=plotingsolution(G,W,'CG',x1,1);
he=plotingsolution(G,W,'error BS',abs(a*p-b),1);
colorbar
he=plotingsolution(G,W,'error ICCG',abs(a*x2-b),2);
colorbar
%h1=plotingsolution(G,W,'DCG',x3,3);
figure(fc+4)
hed=plotingsolution(G,W,'error DICCG',abs(a*x4-b),1);
colorbar
hed=plotingsolution(G,W,'error DICCG1',abs(a*x5-b),2);
colorbar
figure(fc+5)
%h1=plotingsolution(G,W,'CG',x1,1);
hd=plotingsolution(G,W,'diff BS',abs(p-p),1);
colorbar
hd=plotingsolution(G,W,'diff ICCG',abs(x2-p),2);
colorbar
%h1=plotingsolution(G,W,'DCG',x3,3);
figure(fc+6)
hdd=plotingsolution(G,W,'diff DICCG',abs(x4-p),1);
colorbar
hdd=plotingsolution(G,W,'diff DICCG1',abs(x5-p),2);
colorbar
%mkdir([dir 'solx_5_' num2str(num)])
%dir=[dir 'solx_5_' num2str(num) '/'];
mkdir([dir 'sol' ])
dir=[dir 'sol'  '/'];
% Se guarda la grafica en el directorio dir
% file=['x2'];
% filename=[dir file];
% save(filename,file)
% file=['x4'];
% filename=[dir file];
% file=['x5'];
% filename=[dir file];
% save(filename,file) 
%  file='sol';
%   B=[dir  file def '.fig'];
%   saveas(h1(1),B)
%    B=[dir  file def '.jpg'];
%      saveas(h1(1),B)
%    file='soldef';
%      B=[dir  file def '.fig'];
%   saveas(h1d(1),B)
%    B=[dir  file def '.jpg'];
%   saveas(h1d(1),B)
%    file='error';
%   B=[dir  file def '.fig'];
%   saveas(he(1),B)
%    B=[dir  file def '.jpg'];
%   saveas(he(1),B)
%    file='errordef';
%   B=[dir  file def '.fig'];
%   saveas(hed(1),B)
%    B=[dir  file def '.jpg'];
%   saveas(hed(1),B)
%      file='diffbs';
%   B=[dir  file def '.fig'];
%   saveas(hd(1),B)
%    B=[dir  file def '.jpg'];
%   saveas(hd(1),B)
%    file='diffdef';
%   B=[dir  file def '.fig'];
%   saveas(hdd(1),B)
%    B=[dir  file def '.jpg'];
%   saveas(hdd(1),B)
%    file='sol_bs';
%   B=[dir  file def '.fig'];
%   saveas(ht,B)
%   
%    B=[dir  file def '.jpg'];
%   saveas(ht,B)
  
% crear las carpetas para guardar los resultados



% file='eig_pod';
%  B=[dir  file def '.fig'];
%  saveas(h11,B)
%  B=[dir  file def '.jpg'];
%  saveas(h11,B)
   

% file='conv';
%  B=[dir   file def '.fig'];
%  saveas(hl4,B)
%  B=[dir   file def '.jpg'];
%  saveas(hl4,B) 
 
 file='conv_def';
 B=[dir   file def '.fig'];
 saveas(hl5,B)
 B=[dir   file def '.jpg'];
 saveas(hl5,B) 
 if e==1
 file='eigs_mat';
 B=[dir   file def '.fig'];
 saveas(h2(1),B)
 B=[dir  file def '.jpg'];
 saveas(h2(1),B) 
 file='eigs_mat_zoom'; 
  B=[dir   file def '.fig'];
 saveas(h21(1),B)
 B=[dir  file def '.jpg'];
 saveas(h21(1),B)
  file='eigs_mat_1';
 B=[dir   file def '.fig'];
 saveas(h25(1),B)
 B=[dir  file def '.jpg'];
 saveas(h25(1),B) 
 end
%   clearvars('-except','num','e','size1','iteration','tol', 'j')

