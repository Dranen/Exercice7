function donnee = Analyse_donne_convergence(nom)

donnee = lecture_fichier_convergence(nom);
x =  donnee{1};
u = donnee{2};
dt = donnee{3};
dx = donnee{4};
maxf_fin = donnee{5};
energy = donnee{6};
max = donnee{7};
maxf = donnee{8};

Graphique_u;
Graphique_dt_max_energy(dt, max, nom);
Graphique_dx_max_energy(dx, max, nom);
Graphique_dt_maxf(dt, maxf, nom);
Graphique_dx_maxf(dx, maxf, nom);
Graphique_dt_maxf(dx, maxf_fin, nom);
Graphique_dx_maxf(dx, maxf_fin, nom);
Graphique_dt_energie(dt, energy, nom);
Graphique_dx_energie(dx, energy, nom);

end
