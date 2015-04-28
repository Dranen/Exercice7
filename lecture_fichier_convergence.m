function donnee = lecture_fichier_convergence(nom)

a = load([nom '_pos_usquared.dat']);
x = a(:,1);
u = sqrt(a(:,2));

a = load([nom '_wave.dat']);
dt = a(:,1);
dx = a(:,2); 
maxf = a(:,3);
maxf_fin = a(:,4);

a = load([nom '_energy.dat']);
energy = a(:,3:end);

a = load([nom '_maxenergy.dat']);
max = a(:,3:end);

donnee = {x, u, dt, dx, maxf, energy, max, maxf_fin};
end

