function Graphique_dt_maxf(dt,f, nom)
s = size(f);
maxf = zeros(s(2));
for i = 1:s(2)
maxf(i) = max(f(:,i));
end

figure
plot(dt,maxf);
grid on;
xlabel('dt [s]', 'FontSize', 20);
ylabel('Maximum de f [m]', 'FontSize', 20);
saveas(gcf, [nom, '_dt_maxf.fig'])
saveas(gcf, [nom, '_dt_maxf.eps'], 'epsc')
end

