nom = '7_2_convergence_eq3';
xL = 0;
xR = 3;
u = 1;
Ninter = 10;
equation = 3;
question = 0;
CFL = 0.25;
tfinal = 100;
bc_l = 2;
A = 0.1;
omega = 7*pi()/6;
bc_r = 1;
nscan = 100;
Ninter_stop = 1000;

Simulation([nom, '_dx'], Ninter, xL, xR, equation, question, u, u, 1, 1, 1, CFL, tfinal, 0, 0, bc_l, bc_r, A, omega, 3, nscan, 0, 0, Ninter_stop, 1)

Ninter = 100;
CFL = 0.01;
CFL_stop = 0.9;

Simulation([nom, '_dt'], Ninter, xL, xR, equation, question, u, u, 1, 1, 1, CFL, tfinal, 0, 0, bc_l, bc_r, A, omega, 2, nscan, 0, CFL_stop, 0, 1)

donnee1 = Analyse_donne_convergence([nom, '_dx']);
donnee2 = Analyse_donne_convergence([nom, '_dt']);