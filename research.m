% radar frequency (Hz);
freq = linspace(5e9, 100e9, 1001);
c = 3e8;
% radar wavelength (mm);
lambda = c./freq;
% temperature of water in celsius
T_w = [-10 0 10 20];

% Range of normalize diameters (pi D / lambda) to plot
D_n = linspace(0.01, 10, 501);
deltaD = 10;

% Defining variables
nt = length(T_w); % nt = the number of elements in the temperature of the water
nl = length(lambda); % nl = the number of elements in lambda, the radar wavelengths
nd = length(D_n); % nd = the number of elements in the range of normalized diameters

%Fill n_w with values from the interpolated excel equations for each
%temperature
n_w = zeros (nl, nt);
for n = 1: nl
n_w(n, 1) = 9.7979 * exp (-2e-11 .* freq(n));
end
for n = 1: nl
n_w(n, 2) = 9.6233 * exp (-2e-11 .* freq(n));
end
for n = 1: nl
n_w(n, 3) = 9.6649 * exp (-3e-11 .* freq(n));
end
for n = 1: nl
n_w(n, 4) = 8.9104 * exp (-3e-11 .* freq(n));
end

%Fill k_w with values from the interpolated excel equations for each
%temperature
k_w = zeros (nl, nt);
for k = 1: nl
k_w(k, 1) = 0.2666 * reallog((freq(k))) - 3.5211;
end
for k = 1: nl
k_w(k, 2) = -0.3 * reallog((freq(k))) + 9.7807;
end
for k = 1: nl
k_w(k, 3) = 0.7645 * reallog((freq(k))) - 15.309;
end
for k = 1: nl
k_w(k, 4) = 1.0353 * reallog((freq(k))) - 21.881;
end

%Calcualate other components of the complex refractive index needed
m_w = n_w + (i*k_w);
Km = ((m_w.^2)+1)./((m_w.^2)+2);
Km2_w = abs(Km.^2);
ImKm_w = imag(-1 * Km);
 
% set up arrays filling them with zeros
sig_a = zeros(nl, nd); 
sig_s = zeros(nl, nd);
sig_an = zeros(nl, nd);
sig_sn = zeros(nl, nd);
sig_ed = zeros(nl, nd);
att = zeros(nl, nd);
Nd = zeros(nl, nd);
trans = zeros(nl, nd);

%Ask user which temperature index they want to use
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

R = input('Enter the rainfall rate in mm per hour: ');

%Loop through values of diameter and lambda to fill in sigma's and Nd
    for D_n_ind = 1: nd
        for lambi = 1: nl % start of loop to calculate for a given lambi (element of lambda)
            D=lambda(lambi)*D_n(D_n_ind)/pi;
            sig_a(lambi, D_n_ind) = pi^2*D.^3/lambda(lambi)*ImKm_w(lambi, ti); % calculates absorption cross section
            sig_s(lambi, D_n_ind) = 2*pi^5*D.^6/3/lambda(lambi)^4*Km2_w(lambi, ti); % calculates scattering cross section
            sig_an(lambi, D_n_ind) = 4*sig_a(lambi, D_n_ind)./(pi*D.^2); % calculates normalized absorption cross section
            sig_sn(lambi, D_n_ind) = 4*sig_s(lambi, D_n_ind)./(pi*D.^2); % calculates normalized scattering cross section
            Nd (lambi, D_n_ind) = (0.08^-1)*exp(-41*R^(-0.21).*D);
        end     
    end   
    
sig_ed = sig_an + sig_sn; % calculates normalized extinction cross section
norm = 4.34e3 * c; % this number will convert our output to dB/Km
att = norm * Nd .* sig_ed .* deltaD; % calculates attenuation
trans = 1 - att; % calculates transmission


% Print Attenuation Graph
fntsz = 14;
figure(1)
clf
mesh (D_n, freq, att)
set(gca, 'fontsize', fntsz)
zlim ([2 4])
xlabel('Normalized Drop Diameter (\pi D / \lambda)')
ylabel('Frequency (Hz)')
zlabel('Attenuation (Db / Km)')
title(ptit)
  
% Print Transmission Graph
fntsz = 14;
figure(2)
clf
mesh (D_n, freq, trans)
set(gca, 'fontsize', fntsz)
zlim ([-5 -2])
xlabel('Normalized Drop Diameter (\pi D / \lambda)')
ylabel ('Frequency (Hz)')
zlabel('Transmission (dB/Km)')
title(ptit)