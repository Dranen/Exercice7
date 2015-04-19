function Graphique_dt_max_energy(dt, max)
figure
loglog(dt,max);
grid on;
xlabel('dt [s]', 'FontSize', 20);
ylabel('Maximum de l energie [J]', 'FontSize', 20);
saveas(gcf, [nom, '_dt_maxenergy.fig'])
saveas(gcf, [nom, '_dt_maxenergy.eps'], 'epsc')
end

