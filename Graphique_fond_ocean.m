function Graphique_fond_ocean(xocean, xR)

x = [0:1000000];

for i = x
    hocean = sin(pi() * (x(i) - xocean) / (2.0 * (xR - xocean)));
end

end

