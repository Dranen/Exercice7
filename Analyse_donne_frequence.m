function donnee = Analyse_donne_frequence(nom)

donnee = lecture_fichier_frequence(nom);

x =  donnee{1};
u = donnee{2};
omega = donnee{3};
maxenergy = donnee{4};

Graphique_u
Graphique_omega(omega, maxenergy, nom);

end
