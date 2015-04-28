nom = '7_4';
xL = 0;
xR = 1000000;
hplage = 15;
hocean = 6000;
xocean = 500000;
Ninter = 100000;
equation = 2;
question = 2;
CFL = 0.25;
tfinal = 7000;
bc_l = 3;
bc_r = 3;

Graphique_fond_ocean(xocean, xR, hocean, hplage, nom);

Simulation(nom, Ninter, xL, xR, equation, question, 1, 1, hocean, xocean, hplage, CFL, tfinal, 0, 0, bc_l, bc_r, 0, 0, 0, 0, 0, 0, 0, 1000)

donnee = Analyse_donne_unique(nom, 1, 10);