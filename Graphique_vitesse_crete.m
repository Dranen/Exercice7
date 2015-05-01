function Graphique_vitesse_crete(x, t, f, nom)

imax = 1;
maxf = max(f(1,:));
    while(f(1,imax)~=maxf)
       imax = imax+1;
    end
    imaxp = imax;

for j = 2:max(size(t))
    maxf = max(f(j,:));
    while(f(j,imax)~=maxf)
       imax = imax+1;
    end
    v(j-1) = (x(imax)-x(imaxp))/(t(imax) - t(imaxp));
    tg(j-1) = t(j);
    xg(j-1) = x(imax);
    imaxp = imax;
end

figure
plot(xg,v);
grid on;
xlabel('x [m]', 'FontSize', 20);
ylabel('Vitesse de la crete [m/s]', 'FontSize', 20);
saveas(gcf, [nom, '_vitesse_crete_x.fig'])
saveas(gcf, [nom, '_vitesse_crete_x.eps'], 'epsc')

figure
plot(tg,v);
grid on;
xlabel('t [s]', 'FontSize', 20);
ylabel('Vitesse de la crete [m/s]', 'FontSize', 20);
saveas(gcf, [nom, '_vitesse_crete_t.fig'])
saveas(gcf, [nom, '_vitesse_crete_t.eps'], 'epsc')

end

