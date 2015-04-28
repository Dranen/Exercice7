function Graphique_dx_maxf_fin(dx,f, nom)

figure
plot(dx,f);
grid on;
xlabel('dx [m]', 'FontSize', 20);
ylabel('Maximum de f_{fin}[m]', 'FontSize', 20);
saveas(gcf, [nom, '_dx_maxf_fin.fig'])
saveas(gcf, [nom, '_dx_maxf_fin.eps'], 'epsc')
end

