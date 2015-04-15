function Graph_ftx(x, f)
s = [-abs(min(min(f))*1.2) abs(max(max(f))*1.2)];

for i = 1:max(size(f))
    figure(101)
    ax = gca;
    ax.YLim = s;
    ax.YLimMode = 'manual';
    plot(ax, x,f(i,:))
    ax.YLim = s;
    ax.YLimMode = 'manual';
    pause(0.01)
end

end

