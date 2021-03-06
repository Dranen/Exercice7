figure
surf(x,t,f, 'LineStyle', 'none');
grid on;
colorbar;
xlabel('x [m]', 'FontSize', 20);
ylabel('t [s]', 'FontSize', 20);
zlabel('f [m]', 'FontSize', 20);
saveas(gcf, [nom, '_xtf3D.fig'])
saveas(gcf, [nom, '_xtf3D.eps'], 'epsc')

