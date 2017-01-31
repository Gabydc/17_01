close all
clear all
clc

%% Set up grid and petrophysical data
% We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
% data: permeability of 100 mD and porosity of 0.2.
sz=32;
G = cartGrid([sz, sz, 1]);
G = computeGeometry(G);
tol=10^-11;
% Disable gravity
gravity off
% Set up uniform permeability and constant porosity
rock.perm = ones(G.cells.num, 1)*100*milli*darcy;

%inhomogeneus permeability
for i=1:2:8
 rock.perm(1+128*(i-1):128*i)  = repmat(0.01*milli*darcy(), [128, 1]);
end

plotCellData(G, rock.perm)
colorbar
view(0,90) 
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
% is present)
gravity reset off
fluid = initSingleFluid('mu' , 1*centi*poise, ...
                        'rho', 1014*kilogram/meter^3);
display(fluid)

%% Source terms
% The number of wells can be changes.
% pv  = sum(poreVolume(G,rock));
%% Change number of solution, 0=complete problem, 1,2,... different wells
num='4';
W = struct([]);
% % 
W = verticalWell(W, G, rock, 10, 10, [], ...
            'Type', 'bhp' , 'Val', -5*barsa(), ...
            'Radius', 0.1, 'InnerProduct', 'ip_tpf', ...
            'Comp_i', 1);
% % 
W = verticalWell(W, G, rock, 20, 20, [], ...
            'Type', 'bhp' , 'Val', -5*barsa(), ...
            'Radius', 0.1, 'InnerProduct', 'ip_tpf', ...
            'Comp_i', 1);
% 
W = verticalWell(W, G, rock, 20, 10, [], ...
            'Type', 'bhp' , 'Val', 5*barsa(), ...
            'Radius', 0.1, 'InnerProduct', 'ip_tpf', ...
            'Comp_i', 1);

W = verticalWell(W, G, rock, 10, 20, [], ...
            'Type', 'bhp' , 'Val', 5*barsa(), ...
            'Radius', 0.1, 'InnerProduct', 'ip_tpf', ...
            'Comp_i', 1);
% W = verticalWell(W, G, rock, 1, 25, [], ...
%             'Type', 'bhp' , 'Val', -10*barsa(), ...
%             'Radius', 0.1, 'InnerProduct', 'ip_tpf', ...
%             'Comp_i', 1);

% W = verticalWell(W, G, rock, 50, 25, [], ...
%             'Type', 'bhp' , 'Val', 10*barsa(), ...
%             'Radius', 0.1, 'InnerProduct', 'ip_tpf', ...
%             'Comp_i', 1);

%Boundary conditions
bc  = pside([], G, 'YMin',3.*barsa());
bc  = pside(bc, G, 'YMax',0.*barsa());
% Create a initialized state and set initial saturation to phase 1.
 
sol = initState(G, W, 0);
% 
% % Find transmissibility.
% T = computeTrans(G, rock);

for i=1:sz

    Z((1:sz)+sz*(i-1),i)=1;
end
%Z=eye(sz*sz);


% Reference TPFA
solver = DPCG_ICSolverAD('tolerance', tol,'maxIterations', 1000, 'Z',Z);
 %solver = BackslashSolverAD();
% pressureSolver = BackslashSolverAD();
% linsolve = LinearSolverAD('ellipticSolver', pressureSolver);
%linsolver = LinearSolverAD('ellipticSolver', pressureSolver);

fn = @(A, b) solver.solveLinearSystem(A, b);
psolve = @(state) incompTPFA(state, G, hT, fluid, 'wells', W,'MatrixOutput',true,'bc',bc,'LinSolve', fn);

% Implicit transport solver
% tsolve   = @(state, dT) implicitTransport(state, G, dT, rock, ...
%                                                 fluid, 'wells', W);
%                                             
% clf;
% plotGrid(G)
% view(30,50)    
 sol= psolve(sol);
% clf;
% plotCellData(G, sol.pressure)
% colorbar
% view(30,50)    
% [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
% clf;
% % Plot the grid
% plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
% plotCellData(G, )
% colorbar
% % Plot the wells
% plotWell(G, W);
% view(30,50)



 p=sol.pressure;
 A=sol.A(1:G.cells.num,1:G.cells.num);
 b=sol.rhs(1:G.cells.num);
 solb=A\b; 
 figure
plotingsolution(G,W,'s',sol.pressure,1);
colorbar
plotingsolution(G,W,'bs',solb,2);
colorbar
figure
plotingsolution(G,W,'diff',sol.pressure-solb,1);
colorbar
 %%Change dir
 
%  dir='/run/media/taurlin/Sphinx/Doctorado_Delft/Research_programs/15_07_laplace2d_MRST/qfs_1/';
   dir='/home/wagm/cortes/Localdisk/Results/17_01/pruebas/';
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

