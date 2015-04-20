function Graphique_fond_ocean(xocean, xR, hocean, hplage, nom)

x = [0:1000000];
h = zeros(size(x));

for i = x
    if(x(i+1) > xocean)
        h(i+1) = hocean + (hplage-hocean)*((sin(pi() * (x(i+1) - xocean) / (2.0 * (xR - xocean))))^2);
    else
        h(i+1) = hocean;
    end
end
h =-h;
figure
grid on;
plot(x, h)
xlabel('x [m]')
ylabel('h [m]')

saveas(gcf, [nom, '_fond_ocean.fig'])
saveas(gcf, [nom, '_fond_ocean.eps'], 'epsc')

end

