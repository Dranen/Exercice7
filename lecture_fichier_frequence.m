function donnee = lecture_fichier_frequence(nom)

a = load([nom '_pos_usquared.dat']);
x = a(:,1);
u = sqrt(a(:,2));

a = load([nom '_maxenergy.dat']);
omega = a(:,1);
max = a(:,2);

temp = sortrows([omega max], 1);
omega = temp(:,1);
max = temp(:,2);

donnee = {x, u, omega, max};
end

