nom = '7_2';
xL = 0;
xR = 3;
u = 1;
Ninter = 1000;
equation = 0;
CFL = 0.25;
tfinal = 100;
bc_l = 2;
A = 0.1;
omega = 2*pi();
bc_r = 1;
nscan = 1;

Simulation(nom, Ninter, xL, xR, equation, question, u, u, 1, 1, 1, CFL, tfinal, 0, 0, bc_l, bc_r, A, omega, nscan, 0)

Analyse_donne(nom)