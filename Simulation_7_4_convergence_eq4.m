nom = '7_4_convergence_eq4';
Ninter_stop = 10000;
xL = 0;
xR = 1000000;
hplage = 15;
hocean = 6000;
xocean = 500000;
Ninter = 1000;
equation = 4;
question = 2;
CFL = 0.25;
nscan = 100;
tfinal = 7000;
bc_l = 3;
bc_r = 3;

Simulation([nom, '_dx'], Ninter, xL, xR, equation, question, 1, 1, hocean, xocean, hplage, CFL, tfinal, 0, 0, bc_l, bc_r, 0, 0, 3, nscan, 0, 0, Ninter_stop, 1)

Ninter = 10000;
CFL = 0.01;
CFL_stop = 1.1;

Simulation([nom, '_dt'], Ninter, xL, xR, equation, question, 1, 1, hocean, xocean, hplage, CFL, tfinal, 0, 0, bc_l, bc_r, 0, 0, 2, nscan, 0, CFL_stop, 0, 1)


donnee1 = Analyse_donne_convergence([nom, '_dx']);
donnee2 = Analyse_donne_convergence([nom, '_dt']);