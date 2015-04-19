function Graphique_u(x, u)
figure
plot(x,u);
grid on;
xlabel('x [m]', 'FontSize', 20);
ylabel('u [m.s^{-1}]', 'FontSize', 20);
saveas(gcf, [nom, '_u.fig'])
saveas(gcf, [nom, '_u.eps'], 'epsc')
end

