function Graphique_xtf(x,t,f)
figure
contourf(x,t,f, 'LineStyle', 'none');
grid on;
xlabel('x [m]', 'FontSize', 20);
ylabel('t [s]', 'FontSize', 20);
end

