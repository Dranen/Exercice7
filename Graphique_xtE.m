function Graphique_xtE(x,t,E)
figure
contourf(x,t,E);
grid on;
xlabel('x [m]', 'FontSize', 20);
ylabel('t [s]', 'FontSize', 20);
end

