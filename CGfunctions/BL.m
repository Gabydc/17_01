close all
clear all
clc
G = computeGeometry(cartGrid([100,1]));
rock = makeRock(G, 100*milli*darcy, 0.2);
fluid = initSimpleFluid('mu' , [1, 1].*centi*poise, ...
'rho', [1000, 1000].*kilogram/meter^3, 'n', [2,2]);
bc = fluxside([], G, 'Left' , 1, ' sat' , [1 0]);
bc = fluxside(bc, G, 'Right', -1, 'sat' , [0 1]);
hT = computeTrans(G, rock);
rSol = initState(G, [], 0, [0 1]);
rSol = incompTPFA(rSol, G, hT, fluid, 'bc', bc);
rSole = explicitTransport(rSol, G, 10, rock, fluid, 'bc', bc, 'verbose',true);
[rSoli, report] = ...
implicitTransport(rSol, G, 10, rock, fluid, 'bc', bc, 'Verbose', true);