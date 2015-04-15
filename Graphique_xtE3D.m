function Graphique_xtE3D(x,t,E)
figure
surf(x,t,E);
grid on;
xlabel('x [m]', 'FontSize', 20);
ylabel('t [s]', 'FontSize', 20);
zlabel('E [J]', 'FontSize', 20);
end

