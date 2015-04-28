function donnee = lecture_fichier_unique(nom, parse_t, parse_x)

a = load([nom '_pos_usquared.dat']);
x = a(1:parse_x:end,1);
u = sqrt(a(1:parse_x:end,2));


a = load([nom '_wave.dat']);
t = a(1:parse_t:end,1);
f = a(1:parse_t:end,2:parse_x:end);


a = load([nom '_energy.dat']);
energy = a(1:parse_t:end,2);


donnee = {x, u, t, f, energy};
end

