nom = '7_2_frequence_t140';
xL = 0;
xR = 3;
u = 1;
Ninter = 100;
equation = 1;
question = 0;
CFL = 0.25;
tfinal = 140;
bc_l = 2;
A = 1;
omega = 0;
bc_r = 1;
nscan = 2000;
omega_stop = 20;

Simulation(nom, Ninter, xL, xR, equation, question, u, u, 1, 1, 1, CFL, tfinal, 0, 0, bc_l, bc_r, A, omega, 1, nscan, omega_stop, 0, 0, 1)

donnee = Analyse_donne_frequence(nom);