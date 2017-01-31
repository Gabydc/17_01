clear all
close all
%condition to compute the eigenvalues, ife=1, the eigenvalues will be
%computed
e=0;
%The directory where we want to store the matrices, for each size a new
%directory will be created
%dir='/home/wagm/cortes/Localdisk/Research_programs/2016_1_report/1_18_results/';
 dir='/home/wagm/cortes/Localdisk/Results/16_04/prueba/';
 nxi=1;
 nyi=1;
 nx=60;
 ny=220;

  %define the layer of SPE10 that will be used
 layers = 1:85;
 iteration=2000;
 for k=[1]
k
tol=10^-k;
  name=['SPE10'  '_85layers_5w_tol-' num2str(k) ];
  mkdir([dir], name )
  dir=[dir name '/'];
%num is the number of snapshot
 for s=1:5
     %It is neccesary to create a random vector as initial guess and save
     %it for the iteration
    if s==1
        xi=rand(nx*ny*max(layers),1);
        file=[dir 'xi'];
        save(file,'xi')
    end 
    %The grid and the permeability are obtained from the MRST
    cartDims = [  60,  220,   numel(layers)];
    %physDims = [1200, 2200, 2*cartDims(end)] ;   % ft -> m
    G = cartGrid([nx ny numel(layers)]);

    G = computeGeometry(G);
    rock = SPE10_rock(layers);
   % max(rock.perm)/min(rock.perm)
    %The permeability is 
    perm = rock.perm(:,1);
%     perm = reshape(perm,[60 220 1]);
%     permupscaled = sampleFromBox(G,perm);
%     rock.perm=permupscaled;
%     max(rock.perm)
%     min(rock.perm)
%     contrast=max(rock.perm)/min(rock.perm)
%     
%     figure
%     plotCellData(G,log(perm));
%     pause
%% Define well locations
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [100,   100,   100,  100,  300] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  nxi,   nx,     nxi,   nx,   nx/2;
              nyi,   ny,     ny,   nyi,   ny/2];
wname    = {'P1', 'P2', 'P3', 'P4', 'I1'};
sgn      = [ -1 ,  -1 ,  -1 ,  -1 ,   1 ];


wtarget(s)=0;
if s==5
   wtarget(s)=400*barsa();
end



W = [];


for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), [], ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'Sign', sgn(w), 'InnerProduct', 'ip_tpf');
end
% figure
% h = plotCellData(G,log(rock.perm)); 
% h = plotWell(G, W);
% axis equal tight; colormap(jet(128));
% view(60,30)
% file='perm_layer_';
% B=[dir   file  '.fig'];
% saveas(h(1),B)
% B=[dir   file  '.jpg'];
% saveas(h(1),B)
 



%Boundary conditions
%  bc  = pside([], G, 'YMin',10.*barsa());
%  bc  = pside(bc, G, 'YMax',0.*barsa());
sol = initState(G, W, 0);
%% Compute half transmissibilities
% All we need to know to develop the spatial discretization is the reservoir
% geometry and the petrophysical properties. This means that we can compute
% the half transmissibilities without knowing any details about the fluid
% properties and the boundary conditions and/or sources/sinks that will
% drive the global flow:
hT = simpleComputeTrans(G, rock);

%% Fluid model
% When gravity forces are absent, the only fluid property we need in the
% incompressible, single-phase flow equation is the viscosity. However, the
% flow solver is written for general incompressible flow and requires the
% evaluation of a fluid object that can be expanded to represent more
% advanced fluid models. Here, however, we only use a simple fluid object
% that requires a viscosity and a density (the latter is needed when gravity
% is present)
gravity reset off
fluid = initSingleFluid('mu' , 1*centi*poise, ...
                        'rho', 1014*kilogram/meter^3);
display(fluid)

psolve = @(state) incompTPFA(state, G, hT, fluid, 'wells', W,'MatrixOutput',true);

 sol= psolve(sol);
 num=num2str(s);
 
 
 
 
 figure 
plotingsolution_1(G,W,'s',sol.pressure);
h=plotWell(G, W);
p=sol.pressure;
 A=sol.A(1:G.cells.num,1:G.cells.num);
 b=sol.rhs(1:G.cells.num);

if s==1
l=ichol(A);
 file=[dir 'l'];
 save(file,'l')
else
    file=[dir 'l'];
load(file)
end
fc=200;
nn=0;
def=0;
def=['tol-' num2str(k) '_' num2str(s) 'snapshots'];
text = [dir 'iter' def '.txt'];
text1 = [dir 'cond' def '.txt'];

%[x11]=ICCG_1(A,b,xi,iteration,tol,l); 
[x11,iter2,e2,hl2,h2,h21]=ICCG_1(A,b,xi,iteration,tol,l,e,text1,nn,dir,def,fc);
 x(:,s)=x11;
 figure
 [h1]=plotingsolution(G,W,'CG',x(:,s),2);

% colorbar
  
% [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
% clf;
% % Plot the grid
% plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
% plotCellData(G, )
% colorbar
% Plot the wells
% figure

% view(30,50)
%  figure
% plotingsolution(G,W,'s',sol.pressure,1);

 file=[dir 'G'];
 save(file,'G')
 file=[dir 'W' num];
 save(file,'W')
 file=[dir 'b' num];
 save(file,'b')
 file=[dir 'A' num];
 save(file,'A')
 file=[dir 'p' num];
 save(file,'p')

file=['x'];
filename=[dir file];
save(filename,file)
if e==1
 file='eigs_mat';
 B=[dir   file def '.fig'];
 saveas(h2,B)
 B=[dir  file def '.jpg'];
 saveas(h2,B) 
 file='eigs_mat_zoom'; 
  B=[dir   file def '.fig'];
 saveas(h21,B)
 B=[dir  file def '.jpg'];
 saveas(h21,B)
end

 end

%%
%break
 close all
 nn=4;
  iteration=2000;
  for k1=[ 1]
      k1
close all
  fc=100;
  num='5';
  solve_DICCG_f_pod(num,nn,e,iteration,dir,k1)
 
  end
   end
 