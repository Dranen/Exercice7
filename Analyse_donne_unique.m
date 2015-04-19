function donnee = Analyse_donne_unique(nom)

donnee = lecture_fichier(nom);

x =  donnee{1};
u = donnee{2};
t = donnee{3};
f = donnee{4};
energy = donnee{5};
omega = donnee{6};
maxenergy = donnee{7};

Graphique_u(x,u);
Graphique_maxf(x,f);
Graphique_xtf(x,t,f);
Graphique_xtf3D(x,t,f);
Graphique_omega(omega, maxenergy);
Graphique_E(t,energy);
Graph_ftx(x, f);

end
