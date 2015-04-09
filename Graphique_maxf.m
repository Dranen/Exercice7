function Graphique_maxf(x,f)

maxf = zeros(max(size(f)));
for i = 1:max(size(f))
maxf(i) = max(f(:,i));
end

figure
plot(x,maxf);
grid on;
xlabel('x [m]', 'FontSize', 20);
ylabel('Maximum de f [m]', 'FontSize', 20);
end

