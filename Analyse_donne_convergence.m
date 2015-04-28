function donnee = Analyse_donne_convergence(nom)

donnee = lecture_fichier_convergence(nom);
x =  donnee{1};
u = donnee{2};
dt = donnee{3};
dx = donnee{4};
f = donnee{5};
energy = donnee{6};
max = donnee{7};

Graphique_u(x,u, nom);
Graphique_dt_max_energy(dt, max, nom);
Graphique_dx_max_energy(dx, max, nom);
Graphique_dt_maxf(dt, f, nom);
Graphique_dx_maxf(dx, f, nom);
Graphique_dt_energie(dt, energy, nom);
Graphique_dx_energie(dx, energy, nom);

end
