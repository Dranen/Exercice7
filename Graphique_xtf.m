function Graphique_xtf(x,t,f)
figure
contourf(x,t,f);
grid on;
xlabel('x [m]', 'FontSize', 20);
ylabel('t [s]', 'FontSize', 20);
end

