%clear all
% radar frequency (Hz);
freq = linspace(5e9, 50e9, 501);
c = 3e8;
% radar wavelength (mm);
lambda = c./freq;
% temperature of water in celsius
T_w = [-10 0 10 20];
% n from Battan Table 4.1: not sure where this is from
n_w = [8.88 8.14 6.15 4.44; ...
       9.02 7.80 5.45 3.94; ...
       8.99 7.14 4.75 3.45; ...
       NaN  6.48 4.15 3.10];...
% k from Battan Table 4.1
k_w = [0.63 2.00 2.86 2.59; ...
       0.90 2.44 2.90 2.37; ...
       1.47 2.89 2.77 2.04; ...
       NaN  NaN  2.55 1.77];
% Km^2 from Battan Table 4.1
Km2_w = [0.928  0.9275 0.9193 0.8926; ...
     0.9313 0.9282 0.9152 0.8726; ...
     0.9340 0.9300 0.9055 0.8312; ...
     NaN    NaN    0.8902 0.7921];
% Im(-Km) from Battan Table 4.1
ImKm_w = [0.00474 0.01883 0.0471 0.0915; ...
      0.00688 0.0247  0.0615 0.1142; ...
      0.01102 0.0335  0.0807 0.1441; ...
      NaN     NaN     0.1036 0.1713];
% Range of normalize diameters (pi D / lambda) to plot
D_n = linspace(0.01, 10, 201);
deltaD = 10;

% Defining variables
nt = length(T_w); % nt = the number of elements in the temperature of the water
nl = length(lambda); % nl = the number of elements in lambda, the radar wavelengths
nd = length(D_n); % nd = the number of elements in the range of normalized diameters
 
% set up arrays filling them with zeros
sig_a = zeros(nt, nl, nd); 
sig_s = zeros(nt, nl, nd);
sig_an = zeros(nt, nl, nd);
sig_sn = zeros(nt, nl, nd);
sig_ed = zeros(nt, nl, nd);
att = zeros(nt, nl, nd);
Nd = zeros(nt, nl, nd);
trans = zeros(nt, nl, nd);

R = input('Enter the rainfall rate in mm per hour: ');

    for ti = 1: nt % loop to calculate for range of ti (element of temperature) values
        for lambi = 1: nl % start of loop to calculate for a given lambi (element of lambda)
            D = lambda(lambi)*reshape(D_n, 1, 1, nd)/pi;
            sig_a(ti, lambi, :) = pi^2*D.^3/lambda(lambi)*ImKm_w(ti, lambi); % calculates absorption cross section
            sig_s(ti, lambi, :) = 2*pi^5*D.^6/3/lambda(lambi)^4*Km2_w(ti, lambi); % calculates scattering cross section
            sig_an(ti, lambi, :) = 4*sig_a(ti, lambi, :)./(pi*D.^2); % calculates normalized absorption cross section
            sig_sn(ti, lambi, :) = 4*sig_s(ti, lambi, :)./(pi*D.^2); % calculates normalized scattering cross section
            Nd (ti, lambi, :) = (0.08^-1)*exp(-41*R^(-0.21).*D);
        end
    end
sig_ed = sig_an + sig_sn; % calculates normalized extinction cross section
norm = 4.34e3 * 3e8;
att = norm * Nd .* sig_ed .* deltaD; % calculates attenuation
trans = 1 - att;

fprintf(1, 'Possible temperature values\n')
fprintf(1, '1: -10C\n')
fprintf(1, '2:  0C\n')
fprintf(1, '3: 10C\n')
fprintf(1, '4: 20C\n')
ti = input('Enter the temperature index (1 - 4): ');
switch ti
 case 1
  ptit = 'Temperature = -10^{\circ}C';
 case 2
  ptit = 'Temperature = 0^{\circ}C';
 case 3
  ptit = 'Temperature = 10^{\circ}C';
 case 4
  ptit = 'Temperature = 20^{\circ}C';
 otherwise
  fprintf('*** Invalid Selection ... exiting!\n');
  return
end


fntsz = 14;
figure(1)
clf
lambi = 1;
sig = squeeze(att(ti, lambi, :));
loglog(D_n, sig, 'r')
hold on
lambi = 2;
sig = squeeze(att(ti, lambi, :));
loglog(D_n, sig, 'b')
lambi = 3;
sig = squeeze(att(ti, lambi, :));
loglog(D_n, sig, 'k')
lambi = 4;
sig = squeeze(att(ti, lambi, :));
loglog(D_n, sig, 'g')
hold off
set(gca, 'xlim', [1e-2 10])
set(gca, 'ylim', [1e-3 10])
set(gca, 'fontsize', fntsz)
legend('10 cm', '3.21 cm', '1.24 cm', '0.62 cm', ...
       'Location', 'NorthWest')
xlabel('Normalized Drop Diameter (\pi D / \lambda)')
ylabel('Attenuation (dB/Km)')
title(ptit)
  
  fntsz = 14;
figure(2)
clf
lambi = 1;
sig = squeeze(trans(ti, lambi, :));
loglog(D_n, sig, 'r')
hold on
lambi = 2;
sig = squeeze(trans(ti, lambi, :));
loglog(D_n, sig, 'b')
lambi = 3;
sig = squeeze(trans(ti, lambi, :));
loglog(D_n, sig, 'k')
lambi = 4;
sig = squeeze(trans(ti, lambi, :));
loglog(D_n, sig, 'g')
hold off
set(gca, 'xlim', [1e-2 10])
set(gca, 'ylim', [1e-3 10])
set(gca, 'fontsize', fntsz)
legend('10 cm', '3.21 cm', '1.24 cm', '0.62 cm', ...
       'Location', 'NorthWest')
xlabel('Normalized Drop Diameter (\pi D / \lambda)')
ylabel('Transmission (dB/Km)')
title(ptit)