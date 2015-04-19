function Graphique_dt_energie(dt,E)
figure
plot(dt,E);
grid on;
xlabel('dt [s]', 'FontSize', 20);
ylabel('E [J]', 'FontSize', 20);
saveas(gcf, [nom, '_dt_E.fig'])
saveas(gcf, [nom, '_dt_E.eps'], 'epsc')
end

