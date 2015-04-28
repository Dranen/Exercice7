nom = '7_3_frequence';
xL = 0;
xR = 3;
ul = sqrt(0.3);
ur = sqrt(1.7);
Ninter = 1000;
equation = 1;
question = 1;
CFL = 0.25;
tfinal = 100;
bc_l = 2;
A = 0.1;
omega = 0;
bc_r = 1;
nscan = 10000;
omega_stop = 1000;

Simulation(nom, Ninter, xL, xR, equation, question, ul, ur, 1, 1, 1, CFL, tfinal, 0, 0, bc_l, bc_r, A, omega, 1, nscan, omega_stop, 0, 0, 1)

donnee = Analyse_donne_frequence(nom);