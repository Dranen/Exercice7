figure
plot(t,E);
grid on;
xlabel('t [s]', 'FontSize', 20);
ylabel('E [J]', 'FontSize', 20);
saveas(gcf, [nom, '_E.fig'])
saveas(gcf, [nom, '_E.eps'], 'epsc')

