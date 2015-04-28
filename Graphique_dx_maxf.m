function Graphique_dx_maxf(dx,f, nom)
f
figure
plot(dx,f);
grid on;
xlabel('dx [m]', 'FontSize', 20);
ylabel('Maximum de f [m]', 'FontSize', 20);
saveas(gcf, [nom, '_dx_maxf.fig'])
saveas(gcf, [nom, '_dx_maxf.eps'], 'epsc')
end

