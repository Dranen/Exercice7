nom = '7_2';
xL = 0;
xR = 3;
ul = sqrt(0.3);
ur = sqrt(1.7);
Ninter = 100;
equation = 2;
question = 1;
CFL = 0.25;
tfinal = 100;
bc_l = 2;
A = 0.1;
omega = 2*pi();
bc_r = 1;

Simulation(nom, Ninter, xL, xR, equation, question, ul, ur, 1, 1, 1, CFL, tfinal, 0, 0, bc_l, bc_r, A, omega, 0, 0, 0, 0, 0, 1)

donnee = Analyse_donne_unique(nom);