%clear all
% radar wavelength (mm);
lambda = [10 3.21 1.24 0.62]*10; % lambda = radar wavelength in mm 
% temperature of water
T_w = [-10 0 10 20]; % T_w = temperature of water in celsius
% n from Battan Table 4.1 
% WHERE IS THIS TABLE FROM!?!
% n_w = ?
n_w = [8.88 8.14 6.15 4.44; ...
       9.02 7.80 5.45 3.94; ...
       8.99 7.14 4.75 3.45; ...
       NaN  6.48 4.15 3.10];...
% k from Battan Table 4.1
% k_w = ?
k_w = [0.63 2.00 2.86 2.59; ...
       0.90 2.44 2.90 2.37; ...
       1.47 2.89 2.77 2.04; ...
       NaN  NaN  2.55 1.77];
% Km^2 from Battan Table 4.1
% Km2_w = ?
Km2_w = [0.928  0.9275 0.9193 0.8926; ...
     0.9313 0.9282 0.9152 0.8726; ...
     0.9340 0.9300 0.9055 0.8312; ...
     NaN    NaN    0.8902 0.7921];
% Im(-Km) from Battan Table 4.1
% ImKm_w = ?
ImKm_w = [0.00474 0.01883 0.0471 0.0915; ...
      0.00688 0.0247  0.0615 0.1142; ...
      0.01102 0.0335  0.0807 0.1441; ...
      NaN     NaN     0.1036 0.1713];
% Range of normalize diameters (pi D / lambda) to plot
D_n = linspace(0.01, 10, 201); % D_n = normalized diameter range I think
 % Defining variables
nt = length(T_w); % nt = the number of elements in the temperature of the water, don't really need this since we are assuming only one possible value, 0
nl = length(lambda); % nl = the number of elements in lambda, the radar wavelengths
nd = length(D_n); % nd = the number of elements in the range of normalized diameters
 
% set up the arrays for sigma a and s and sigma a and s normalized, filling
% them with all zeros
sig_a = zeros(nt, nl, nd); 
sig_s = zeros(nt, nl, nd);
sig_an = zeros(nt, nl, nd);
sig_sn = zeros(nt, nl, nd);
sig_ed = zeros(nt, nl, nd);
att_d = zeros(nt, nl, nd);
att = zeros(nt, nl, nd);

%lambi = 1; %just to see if it'll work
d = 0.001; %ditto

deltaD = 1; % can change or ask for input, went with 1mm because of what Chilson said in one of our meetings
R = input('Enter the rainfall rate in mm per hour: ');

while (d < deltaD)
    for ti = 1: nt % loop to calculate for range of ti (element of temperature) values
        for lambi = 1: nl % start of loop to calculate cross sections for a given lambi (element of lambda)
        sig_a(ti, lambi, :) = pi^2*d^3/lambda(lambi)*ImKm_w(ti, lambi); % calculates absorption cross section
        sig_s(ti, lambi, :) = 2*pi^5*d^6/3/lambda(lambi)^4*Km2_w(ti, lambi); % calculates scattering cross section
        sig_an(ti, lambi, :) = 4*sig_a(ti, lambi, :)/(pi*d^2); % calculates normalized absorption cross section
        sig_sn(ti, lambi, :) = 4*sig_s(ti, lambi, :)/(pi*d^2); % calculates normalized scattering cross section
        end
    end
Nd = (0.08^-1)*exp(-41*R^(-0.21)*d);
sig_ed = sig_an + sig_sn; % calculates normalized extinction cross section
att_d = Nd * sig_ed * deltaD; % calculates normalized attenuation?
att = att + att_d; % calculates total attenuation by adding to what we already had from previous circuits of the loop
d = d + 0.001; % d++ for loop purposes
end
  
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
sig = squeeze(att_d(ti, lambi, :));
loglog(D_n, sig, 'r')
hold on
lambi = 2;
sig = squeeze(att_d(ti, lambi, :));
loglog(D_n, sig, 'b')
lambi = 3;
sig = squeeze(att_d(ti, lambi, :));
loglog(D_n, sig, 'k')
lambi = 4;
sig = squeeze(att_d(ti, lambi, :));
loglog(D_n, sig, 'g')
hold off
set(gca, 'xlim', [1e-2 10])
set(gca, 'ylim', [1e-3 10])
set(gca, 'fontsize', fntsz)
legend('10 cm', '3.21 cm', '1.24 cm', '0.62 cm', ...
       'Location', 'NorthWest')
xlabel('Normalized Drop Diameter (\pi D / \lambda)')
ylabel('Attenuation (dBkm^-1)')
title(ptit)
 
reply = input('Create pdf figure? (n/y): ', 's');
if isempty(reply), reply = ' '; end
if upper(reply) == 'Y'
  fname1 = 'prob_1_fig1.pdf';
  fprintf(1, 'Creating file %s\n', fname1);
  orient(1, 'landscape')
  print(1, '-dpdf', fname1)
end