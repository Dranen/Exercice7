function Graphique_E(t,E, nom)
figure
plot(t,E);
grid on;
xlabel('x [m]', 'FontSize', 20);
ylabel('E [J]', 'FontSize', 20);
saveas(gcf, [nom, '_E.fig'])
saveas(gcf, [nom, '_E.eps'], 'epsc')
end

