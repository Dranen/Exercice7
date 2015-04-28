function Graphique_dt_maxf_fin(dt,f, nom)


figure
plot(dt,f);
grid on;
xlabel('dt [s]', 'FontSize', 20);
ylabel('Maximum de f_{fin} [m]', 'FontSize', 20);
saveas(gcf, [nom, '_dt_maxf_fin.fig'])
saveas(gcf, [nom, '_dt_maxf_fin.eps'], 'epsc')
end

