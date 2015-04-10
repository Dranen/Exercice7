function Graphique_xtf3D(x,t,f)
figure
surf(x,t,f);
grid on;
xlabel('x [m]', 'FontSize', 20);
ylabel('t [s]', 'FontSize', 20);
zlabel('f [m]', 'FontSize', 20);
end

