for i = 1:max(size(x))
    if(x(i) > xocean)
        Ath(i) = 0.5*(hocean^(0.25))/((hocean + (hplage-hocean)*((sin(pi() * (x(i) - xocean) / (2.0 * (xR - xocean))))^2))^0.25);
    else
        Ath(i) = 0.5*(hocean^(0.25))/(hocean^(0.25));
    end
    Anum(i) = (max(f(:,i))-min(f(:,i)));
end

figure
hold all
grid on
plot(x,Ath)
plot(x,Anum)
xlabel('X [m]')
ylabel('A [m]')
legend('Amplitude théorique','Amplitude numérique')
