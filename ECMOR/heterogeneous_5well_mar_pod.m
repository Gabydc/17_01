close all
clear all
clc

%% Set up grid and petrophysical data
% We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
% data: permeability of 100 mD and porosity of 0.2.

pw=1;
size=32*2^pw;
% w1=10*2^pw;
% w2=2*w1;
nxi=1;
nyi=1;
nx=size;
ny=size;
iteration=500;
e=0;
num='5';

for k=11
    k
tol=10^-k;

for per=[3]
   
close all
per
%Create the directory
dir='/home/wagm/cortes/Localdisk/Results/16_03/POD/';
%dir='/dev/media/Sphinx/Doctorado_Delft/marzo16/2016_03/Heterogeneous/POD/';
%dir='/dev/media/Sphinx/Doctorado_Delft/Research_programs/2016_1_report/2016_02/heterogeneous/matrixE/';
folder=['size_'   num2str(size) 'perm_' num2str(per) '_5wells' num2str(k)  ];
mkdir([dir], folder)
dir = [dir folder '/'];
G = cartGrid([size, size, 1]);
G = computeGeometry(G);
% Disable gravity
gravity off

% Set up uniform permeability and constant porosity
rock.perm = ones(G.cells.num, 1)*1*milli*darcy;
%inhomogeneus permeability
 lsize=round(size*size/8); 
 
% for i=2:2:6
%  %for i=2
%  if i==2
%      init1=lsize*(i-i/4);
%   rock.perm(1+init1:init1+lsize)  = repmat(10^(-per)*milli*darcy(), [lsize, 1]);
%  else
%      init=init1+lsize*(i-2);
%      rock.perm(init+1:init+lsize)  = repmat(10^(-per)*milli*darcy(), [lsize, 1]);
%  end
%  end

for i=1:2:8
 rock.perm(1+lsize*(i-1):lsize*i)  = repmat(10^(-per)*milli*darcy(), [lsize, 1]);
end


 %Visualize the plot
%  hp=plotCellData(G, rock.perm);
%  
%  view(0,90) 
%  file='perm';
%  axis equal tight
%  B=[dir   file  '.fig'];
%  saveas(hp,B)
%  B=[dir   file   '.jpg'];
%  saveas(hp,B) 
% return
 %break
% rock.perm(1:G.cells.num)  = repmat(10*milli*darcy(), [128, 1]);
rock.poro = ones(G.cells.num, 1)*0.2;
% rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
% rock.poro = repmat(0.5            , [G.cells.num, 1]);



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
% is present)3 bar
gravity reset off
fluid = initSingleFluid('mu' , 1*centi*poise, ...
                        'rho', 1014*kilogram/meter^3);
display(fluid);
%Boundary conditions
% bc  = pside([], G, 'YMin',3.*barsa());
% bc  = pside(bc, G, 'YMax',0.*barsa());

%% Source terms
% The number of wells can be changes.
% pv  = sum(poreVolume(G,rock));
%% Change number of solution, 1=complete problem, 1,2,... different wells



%% Define well locations, 4 wells

for s=1:5
    well(1:4)=-1;
    well(5)=3;
num=num2str(s);
    if s==1
    xi=rand(nx*ny,1);
         file=[dir 'xi'];
        save(file,'xi')
        well(1)=0;
    else if s==2
        well(2)=0;
    else if s==3
        well(3)=0;
        else if s==4
               well(4)=0;
            else
        well(5)=4;
            end
        end
        end
    end
%        well
%        pause
%        
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [well(1),   well(2),   well(3),   well(4), well(5)] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  nxi,   nxi,     nx,   nx, nx/2;
              nyi,   ny,     nyi,   ny, ny/2];
wname    = {'W1', 'W2', 'W3', 'W4', 'W1'};
sgn      = [ 1 ,  1 ,  1 ,  1 ,1 ];
W = [];        
for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), [], ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'Sign', sgn(w), 'InnerProduct', 'ip_tpf');
end

% 
%   h1=plotCellData(G, rock.perm);
%  h1=plotWell(G, W);
% axis('tight')
%   view(0,90) 
% colormap('winter')
% file='perm';
%  B=[dir   file '.fig'];
%  saveas(h1(1),B)
%  file='perm';
%  B=[dir   file  '.jpg'];
%  saveas(h1(1),B)
 

% Create a initialized state and set initial saturation to phase 1.

sol = initState(G, W, 0);

% % Find transmissibility.
% T = computeTrans(G, rock);

% Reference TPFA
psolve = @(state) incompTPFA(state, G, hT, fluid, 'wells', W,'MatrixOutput',true);

sol= psolve(sol);
p=sol.pressure;
 A=sol.A(1:G.cells.num,1:G.cells.num);
 b=sol.rhs(1:G.cells.num);
% clf;



%plotWell(G, W);
view(0,90)  
%%Sol ICCG
if s==1
l=ichol(A);
 file=[dir 'l'];
 save(file,'l')
else
    file=[dir 'l'];
load(file)
end

[x11,hc0,hc01,hc02]=ICCG_0(A,b,xi,iteration,tol,l,e); 


if e==1
     file='eigs_a';
  B=[dir  file '.fig'];
  saveas(hc01(1),B)
   B=[dir  file '.jpg'];
  saveas(hc01(1),B)
end


 x(:,s)=x11;
 x1(:,s)=sol.pressure;
 res(:,s)=x(:,s)-x1(:,s);

figure(s+per+200)
[ht]=plotingsolution(G,W,'s', sol.pressure,1) ;
colorbar
title('backslash')
 [h1]=plotingsolution(G,W,'ICCG',x(:,s),2);
 colorbar
title('ICCG')
figure(s+per+300)
 [h1]=plotingsolution(G,W,'bs-ICGCG',res(:,s),2);
 colorbar
title('bs-ICCG')

% well
% pause
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

end

file=['x'];
filename=[dir file];
save(filename,file)
file=['x1'];
filename=[dir file];
save(filename,file) 
% return
close all


 iteration=500;
 for k1=11
     k1
 fc=100;
 %e=1;
 for nn=[4]
     
 
 solve_DICCG_f_pod(num,nn,e,iteration,dir,k1);

 end
 end
end
end