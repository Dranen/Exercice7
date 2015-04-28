nom = '7_2_frequence';
xL = 0;
xR = 3;
u = 1;
Ninter = 1000;
equation = 1;
question = 0;
CFL = 0.25;
tfinal = 100;
bc_l = 2;
A = 0.1;
omega = 0;
bc_r = 1;
nscan = 10000;
omega_stop = 50;

Simulation(nom, Ninter, xL, xR, equation, question, u, u, 1, 1, 1, CFL, tfinal, 0, 0, bc_l, bc_r, A, omega, 1, nscan, omega_stop, 0, 0, 1)

donnee = Analyse_donne_frequence(nom);