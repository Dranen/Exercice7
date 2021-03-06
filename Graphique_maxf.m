s = size(f);
%maxf = zeros(s(2));
for i = 1:s(2)
maxf(i) = max(f(:,i));
end

figure
plot(x,maxf)
grid on
xlabel('x [m]', 'FontSize', 20)
ylabel('Maximum de f [m]', 'FontSize', 20)
saveas(gcf, [nom, '_maxf.fig'])
saveas(gcf, [nom, '_maxf.eps'], 'epsc')

