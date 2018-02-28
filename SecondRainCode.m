clear all
% radar frequency (GHz);
freq = linspace(1, 1000, 1000);
nf = length(freq); % nl = the number of elements in lambda, the radar wavelengths
% fill with zeros
att = zeros(nf, 6); 
trans = zeros(nf, 6);
% declare all values from tables needed to make calculations
mk1 = -0.18961;
ck1 = 0.71147;
mk2 = -0.16398;
ck2 = 0.63297;
ma1 = 0.67849;
ca1 = -1.95537;
ma2 = -0.053739;
ca2 = 0.83433;
a1 = [-5.33980 -0.35351 -0.23789 -0.94158];
b1 = [-0.10008 1.26970 0.86036 0.64552];
c1 = [1.13098 0.454 0.15354 0.16817];
a2 = [-3.80595 -3.44965 -0.39902 0.50167];
b2 = [0.56934 -0.22911 0.73042 1.07319];
c2 = [0.81061 0.51059 0.11899 0.27195];
a3 = [-0.14318 0.29591 0.32177 -5.3761 16.1721];
b3 = [1.82442 0.77564 0.63773 -0.96230 -3.29980];
c3 = [-0.55187 0.19822 0.13164 1.47828 3.4399];
a4 = [-0.07771 0.56727 -0.20238 -48.2991 48.5833];
b4 = [2.3384 0.95545 1.1452 0.791669 0.791459];
c4 = [-0.76284 0.54039 0.26809 0.116226 0.116479];
R = [5, 10, 15, 20, 25, 30]; % Rain rates in mm/h
theta = 0; % can be 0 or 90
tilt = 90; % can be 0 or 90
for i = 1 : 6 % for each rain rate
    for f = 1 : nf % for each frequency 
        r = R(i); 
        khsum = 0;
        kvsum = 0;
            for j = 1 : 4
                kh = a1(j)*exp(-1*((log10(freq(f))-b1(j))/c1(j))^2)+mk1*log10(freq(f))+ck1;
                kv = a2(j)*exp(-1*((log10(freq(f))-b2(j))/c2(j))^2)+mk2*log10(freq(f))+ck2;
                khsum = khsum + kh;
                kvsum = kvsum + kv;
            end
        khsum = 10^khsum;
        kvsum = 10^kvsum;
        ahsum = 0;
        avsum = 0;
            for j = 1 : 5
                ah = a3(j)*exp(-1*((log10(freq(f))-b3(j))/c3(j))^2)+ma1*log10(freq(f))+ca1;
                av = a4(j)*exp(-1*((log10(freq(f))-b4(j))/c4(j))^2)+ma2*log10(freq(f))+ca2;
                ahsum = ahsum + ah;
                avsum = avsum + av;
            end
        k = (khsum + kvsum + (khsum - kvsum) * (cos(theta)^2)* cos(2*tilt))/2;
        a = ((khsum*ahsum) + (kvsum*avsum) + ((khsum*ahsum) - (kvsum*avsum))*(cos(theta)^2)*cos(2*tilt))/(2*k);
        att(f, i) = k*r^a; % calculates attenuation
        trans(f, i) = 1 - att(f, i); % calculates transmission
    end
end    
% Print Attenuation Graph
fntsz = 14;
figure(1)
clf
plot (freq, att(:, 1), 'b', 'DisplayName', '5mm')
hold on
plot (freq, att(:, 2), 'r', 'DisplayName', '10mm')
plot (freq, att(:, 3), 'k', 'DisplayName', '15mm')
plot (freq, att(:, 4), 'g', 'DisplayName', '20mm')
plot (freq, att(:, 5), 'r--', 'DisplayName', '25mm')
plot (freq, att(:, 6), 'b--', 'DisplayName', '30mm')
ylim ([2 10])
lgd = legend('location', 'northwest', 'orientation', 'vertical')
title(lgd, 'Rain Rate')
set(gca, 'fontsize', fntsz)
xlabel('Frequency (GHz)')
ylabel('Attenuation (dB / Km)')
hold off
% Print Transmission Graph
figure(2)
clf

plot (freq, trans(:, 1), 'b', 'DisplayName', '5mm')
hold on
plot (freq, trans(:, 2), 'r', 'DisplayName', '10mm')
plot (freq, trans(:, 3), 'k', 'DisplayName', '15mm')
plot (freq, trans(:, 4), 'g', 'DisplayName', '20mm')
plot (freq, trans(:, 5), 'r--', 'DisplayName', '25mm')
plot (freq, trans(:, 6), 'b--', 'DisplayName', '30mm')
ylim ([-10 -2])
lgd = legend('location', 'southwest', 'orientation', 'vertical')
title(lgd, 'Rain Rate')
set(gca, 'fontsize', fntsz)
xlabel('Frequency (GHz)')
ylabel ('Transmission (dB / Km)')
hold off