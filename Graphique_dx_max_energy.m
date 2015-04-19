function Graphique_dx_max_energy(dx, max)
figure
loglog(dx,max);
grid on;
xlabel('dx [s]', 'FontSize', 20);
ylabel('Maximum de l energie [J]', 'FontSize', 20);
saveas(gcf, [nom, '_dx_maxenergy.fig'])
saveas(gcf, [nom, '_dx_maxenergy.eps'], 'epsc')
end

