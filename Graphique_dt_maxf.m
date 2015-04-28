function Graphique_dt_maxf(dt,f, nom)
f
figure
plot(dt,f);
grid on;
xlabel('dt [s]', 'FontSize', 20);
ylabel('Maximum de f [m]', 'FontSize', 20);
saveas(gcf, [nom, '_dt_maxf.fig'])
saveas(gcf, [nom, '_dt_maxf.eps'], 'epsc')
end

