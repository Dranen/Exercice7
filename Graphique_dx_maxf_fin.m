function Graphique_dx_maxf_fin(dx,f, nom)
s = size(f);
maxf = zeros(s(2));
for i = 1:s(2)
maxf(i) = max(f(:,i));
end

figure
plot(dx,maxf);
grid on;
xlabel('dx [m]', 'FontSize', 20);
ylabel('Maximum de f [m]', 'FontSize', 20);
saveas(gcf, [nom, '_dx_maxf_fin.fig'])
saveas(gcf, [nom, '_dx_maxf_fin.eps'], 'epsc')
end

