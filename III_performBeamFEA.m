clear
close all

hmax        = 2e-2;
fmodelType  = 'exponential';
tfinal      = 0.1;
dt          = 1e-4;
eps         = 1e-3;
beta        = 1e-2;

performStickSlipBeamFEA(hmax, fmodelType, tfinal, dt, eps, beta)