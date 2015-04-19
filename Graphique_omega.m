function Graphique_omega(omega, maxenergy)
figure
plot(omega,maxenergy);
grid on;
xlabel('\omega', 'FontSize', 20);
ylabel('Energie maximale [J]', 'FontSize', 20);
saveas(gcf, [nom, '_omega.fig'])
saveas(gcf, [nom, '_omega.eps'], 'epsc')
end

