function donnee = lecture_fichier_unique(nom)

a = load([nom '_pos_usquared.dat']);
x = a(:,1);
u = sqrt(a(:,2));


a = load([nom '_wave.dat']);
t = a(:,1);
f = a(:,2:end);


a = load([nom '_energy.dat']);
energy = a(:,2:end);


donnee = {x, u, t, f, energy};
end

