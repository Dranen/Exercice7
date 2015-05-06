function Graphique_dx_energie(dx,E, nom)
figure
plot(dx,E);
grid on;
xlabel('dx [m]', 'FontSize', 20);
ylabel('E [J]', 'FontSize', 20);
saveas(gcf, [nom, '_dx_E.fig'])
saveas(gcf, [nom, '_dx_E.eps'], 'epsc')
end

