figure
contourf(x,t,f, 'LineStyle', 'none');
grid on;
xlabel('x [m]', 'FontSize', 20);
ylabel('t [s]', 'FontSize', 20);
saveas(gcf, [nom, '_xtf.fig'])
saveas(gcf, [nom, '_xtf.eps'], 'epsc')

