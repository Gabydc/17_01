close all
clear all
clc
%Linear Solvers
try
    mrstModule add agmg
    1
    %pressureSolver = BackslashSolverAD();
    %pressureSolver = AGMGSolverAD('tolerance', 1e-4);
    %pressureSolver = GMRES_ILUSolverAD('tolerance', 1e-4,'maxIterations', 1000);
     pressureSolver = PCG_ICSolverAD('tolerance', 1e-4,'maxIterations', 1000);
catch
  %  pressureSolver = AGMGSolverAD('tolerance', 1e-4);
    %pressureSolver = BackslashSolverAD();
end

%% Grid and rock properties
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-props ad-core % AD framework
mrstModule add ad-blackoil      % Three phase simulator

[nx,ny,nz] = deal( 20, 20, 20);
[Lx,Ly,Lz] = deal(200, 200, 50);
G = cartGrid([nx, ny, nz], [Lx, Ly, Lz ]);
G = computeGeometry(G);
nzz = 8;
Gt = cartGrid([1, 1, 1], [200, 200, nzz*Lz/nz ]);
Gt = computeGeometry(Gt); Gb = Gt;
Gb.nodes.coords(:,3) = [30; 30; 30; 30; 50; 50; 50; 50];

%% Layers

bottomlayer = G.cells.num-nx*ny*nzz+1:G.cells.num;
toplayer = 1:nx*ny*nzz;
middlelayer = setdiff(1:G.cells.num,[toplayer,bottomlayer]);
rock.perm = repmat(30*milli*darcy, [G.cells.num, 1]);
rock.perm(toplayer) = 150*milli*darcy;
rock.perm(bottomlayer) = 150*milli*darcy;
rock.poro = repmat(0.2 , [G.cells.num, 1]);
rock.poro(toplayer) = 0.3;
rock.poro(bottomlayer) = 0.3;
G.rock = rock;
figure; plotToolbar(G,rock); 
view(-45,30)
plotGrid(Gb, 'facealpha',0); plotGrid(Gt,'facealpha',0);
axis equal tight
gravity on


%% Define fluid properties
% Define a two-phase fluid model without capillarity. The fluid model has
% the values:
%
% * densities: [rho_w, rho_o] = [1000 700 250] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 5 0.2] cP.
% * corey-coefficient: [2, 2] = [2 2 2].

fluid = initSimpleADIFluid('phases', 'OW', ...
                            'mu',      [1, 5]*centi*poise, ...
                            'rho',     [1000, 700]*kilogram/meter^3 , ...
                            'n',       [2, 2]);

fluid0 = fluid; % Incompressible fluid

 pRef = 100*barsa;
 c_w = 1e-8/barsa;
 c_o = 1e-5/barsa;
% c_g = 1e-3/barsa;
% 
% % Add compressibility to the fluid
% fluid.bW = @(p) exp((p - pRef)*c_w);
% fluid.bO = @(p) exp((p - pRef)*c_o);
% fluid.bG = @(p) exp((p - pRef)*c_g);

%% Define three-phase compressible flow model

% We define a three phase blackoil model without dissolved gas or vaporized
% oil. This is done by first instansiating the blackoil model, and then
% manually passing in the internal transmissibilities and the topological
% neighborship from the embedded fracture grid.

% model = ThreePhaseBlackOilModel(G, [], fluid, 'disgas', false, 'vapoil', false);
 model.operators = setupOperatorsTPFA(G, G.rock);

%% Add wells
% We have a single horizontal gas injector that injects a total of one pore
% volume of gas at surface conditions over a period of 'totTime' years. A
% horizontal producer is also added to the top and set to bottom hole
% pressure controls of 100 bar.

totTime = 1*year;
tpv = sum(model.operators.pv);

% inj = 1:nx*ny:nx*ny*(nz-1)+1;
% prod = nx*ny:nx*ny:nx*ny*nz;

inj = nx*ny*(nz-1)+1:nx:nx*ny*nz;
prod = nx:nx:nx*ny;

W = addWell([], G, rock, inj, 'Type', 'rate', 'Val', tpv/totTime, 'Name', 'I', ...
    'Sign', 1, 'Comp_i', [0, 1]);

W = addWell(W, G, rock, prod, 'Type', 'bhp', 'Val', 100*barsa, 'Name', 'P', ...
    'Sign', -1, 'Comp_i', [1, 1]/2);

clf; plotCellData(G, rock.perm, 'facealpha', 0.3, 'edgealpha', 0.05); plotWell(G,W);
axis tight equal; view(-135,30);

%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water and oil initially. We also set up a simple-time step strategy that
% ramps up gradually towards 30 day time-steps.

s0 = [0.8, 0.2];
state  = initResSol(G, pRef, s0);
dt = rampupTimesteps(totTime, 30*day, 5);
% schedule = simpleSchedule(dt, 'W', W);
% 
% %% Simulate problem
% 
% [ws, states, report] = simulateScheduleAD(state, model, schedule, ...
%    'afterStepFn', getPlotAfterStep(state, model, schedule));

%% Create and solve the same problem without any compressibility
% Fractures can have a large impact on fluid displacement, but are often
% associated with great uncertainty in terms of locations, orientation and
% permeability. To demonstrate the effect of fractures on the transport, we
% create another problem where the fractures themselves have been omitted.


% Same model, but incompressible fluid
model0 = TwoPhaseOilWaterModel(G, rock, fluid0);
%model0 = ThreePhaseBlackOilModel(G, rock, fluid0, 'disgas', false, 'vapoil', false);
state0  = initResSol(G, pRef, s0);

% Make a copy of the wells in the new grid
W0 = [];
for i = 1:numel(W)
    w = W(i);
    W0 = addWell(W0, G, rock, w.cells, 'type', w.type, 'val', w.val, ...
                 'wi', w.WI, 'comp_i', w.compi);
end
% Set up and simulate the schedule
schedule0 = simpleSchedule(dt, 'W', W0);
[ws0, states0, report0] = simulateScheduleAD(state0, model0, schedule0, ...
   'afterStepFn', getPlotAfterStep(state0, model0, schedule0),'LinearSolver',pressureSolver);

%% Plot the results
% We plot the production curves for the two wells and compare between the
% two cases, as well as the final gas saturation.
% 
% tm = report.ReservoirTime/day;
% flds = {'qGs', 'qWs', 'qOs'};
% names = {'Gas production', 'Water Production', 'Oil production'};
% sgn = vertcat(W.sign);
% inj = find(sgn > 0);
% prod = find(sgn <= 0);
% colors = lines(nnz(prod));
% 
% for i = 1:numel(flds)
%     figure(i); clf
%     q0 = -getWellOutput(ws0, flds{i}, prod)*day;
%     q =  -getWellOutput(ws, flds{i}, prod)*day;
%     l = {};
%     for j = 1:numel(prod)
%         c = colors(j, :);
%         
%         plot(tm, q(:, j), 'color', c);
%         hold on
%         plot(tm, q0(:, j), '--', 'color', c, 'linewidth', 2);
%         title(names{i})
%         ylabel('Production [m^3/day]');
%         xlabel('Time [days]')
%         n = W(prod(j)).name;
%         l = [l, {[n, ' - Incompressible'], [n, ' - Compressible']}]; %#ok
%     end
%     legend(l)
% end
% 
% figure;
% subplot(2, 1, 1)
% plotCellData(G, states{end}.s(:, 3), 'EdgeColor', 'none')
% view(30,10);
% axis equal tight
% caxis([0, 1])
% title('Compressible')
% subplot(2, 1, 2)
% plotCellData(G, states0{end}.s(:, 3), 'EdgeColor', 'none')
% view(30,10);
% title('Incompressible')
% axis equal tight
% caxis([0, 1])
% % 
% %% Interactive plotting
% % Compare the well cruves
% 
% plotWellSols({ws, ws0}, report.ReservoirTime, 'datasetnames', {'Compressible', 'Incompressible'})
% 
% figure;
% plotToolbar(G, states)
% axis equal tight, view(30,10);
% plotWell(G,W);
% plotGrid(Gb, 'facealpha',0); plotGrid(Gt,'facealpha',0);
% title('Compressible')
% 
% figure;
% plotToolbar(G, states0)
% axis equal tight, view(30,10);
% plotGrid(Gb, 'facealpha',0); plotGrid(Gt,'facealpha',0);
% plotWell(G,W);
% title('Incompressible')