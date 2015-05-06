function Graphique_vitesse_crete(x, t, f, nom)

for j = 2:(max(size(x))-1)
    maxf = max(f(:,j));
    imax = 2;
    while  imax < (max(size(t))-1) && f(imax,j)~=maxf
       imax = imax+1;
    end

        p = polyfit([t(imax-1) t(imax) t(imax+1)], [f(imax-1, j) f(imax, j) f(imax+1, j)], 2);
        xg(j-1) = x(j);
        tg(j-1) = roots(polyder(p));

end

for j = 1:(max(size(xg))-1)
v(j) = abs((xg(j+1)-xg(j))/(tg(j+1)-tg(j)));
end

figure
plot(xg(1:end-1),v);
grid on;
xlabel('x [m]', 'FontSize', 20);
ylabel('Vitesse de la crete [m/s]', 'FontSize', 20);
saveas(gcf, [nom, '_vitesse_crete_x.fig'])
saveas(gcf, [nom, '_vitesse_crete_x.eps'], 'epsc')

end