function donnee = Analyse_donne_unique(nom)

donnee = lecture_fichier_unique(nom);

x =  donnee{1};
u = donnee{2};
t = donnee{3};
f = donnee{4};
energy = donnee{5};

Graphique_u(x,u, nom);
Graphique_maxf(x,f, nom);
Graphique_xtf(x,t,f, nom);
Graphique_xtf3D(x,t,f, nom);
Graphique_E(t,energy, nom);
Graph_ftx(x, f, nom);

end
