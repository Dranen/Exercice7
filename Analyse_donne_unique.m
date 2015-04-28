function donnee = Analyse_donne_unique(nom, parse_t, parse_x)

donnee = lecture_fichier_unique(nom, parse_t, parse_x);

x =  donnee{1};
u = donnee{2};
t = donnee{3};
f = donnee{4};
E = donnee{5};

Graphique_u
Graphique_maxf
Graphique_xtf
Graphique_xtf3D
Graphique_E
Graphique_ft

end
