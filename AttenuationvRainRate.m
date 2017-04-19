clear all
% radar frequency (Hz);
freq = linspace(0e9, 100e9, 10000);
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
att1 = zeros(nl, 1);
att2 = zeros(nl, 1);
att3 = zeros(nl, 1);
att4 = zeros(nl, 1);
att5 = zeros(nl, 1);
att6 = zeros(nl, 1);
Nd = zeros(nl, nd);
trans1 = zeros(nl, 1);
trans2 = zeros(nl, 1);
trans3 = zeros(nl, 1);
trans4 = zeros(nl, 1);
trans5 = zeros(nl, 1);
trans6 = zeros(nl, 1);

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

   R = [5, 10, 15, 20, 25, 30];
   
    for lambi = 1 : nl
        for D_n_ind = 1: nd
            r = R(1);
            D=lambda(lambi)*D_n(D_n_ind)/pi;
            sig_a = pi^2*D^3/lambda(lambi)*ImKm_w(lambi, ti); % calculates absorption cross section
            sig_s = 2*pi^5*D^6/3/lambda(lambi)^4*Km2_w(lambi, ti); % calculates scattering cross section
            sig_an = 4*sig_a/(pi*D^2); % calculates normalized absorption cross section
            sig_sn = 4*sig_s/(pi*D^2); % calculates normalized scattering cross section
            Nd = (0.08^-1)*exp(-41*r^(-0.21)*D); 
            sig_ed = sig_an + sig_sn; % calculates normalized extinction cross section
            att1(lambi) = att1(lambi) + (4.34e3 * Nd * sig_ed * deltaD); % calculates attenuation
            trans1(lambi) = 1 - att1(lambi); % calculates transmission
        end
    end
    
  for lambi = 1 : nl
        for D_n_ind = 1: nd
            r = R(2);
            D=lambda(lambi)*D_n(D_n_ind)/pi;
            sig_a = pi^2*D^3/lambda(lambi)*ImKm_w(lambi, ti); % calculates absorption cross section
            sig_s = 2*pi^5*D^6/3/lambda(lambi)^4*Km2_w(lambi, ti); % calculates scattering cross section
            sig_an = 4*sig_a/(pi*D^2); % calculates normalized absorption cross section
            sig_sn = 4*sig_s/(pi*D^2); % calculates normalized scattering cross section
            Nd = (0.08^-1)*exp(-41*r^(-0.21)*D); 
            sig_ed = sig_an + sig_sn; % calculates normalized extinction cross section
            att2(lambi) = att2(lambi) + (4.34e3 * Nd * sig_ed * deltaD); % calculates attenuation
            trans2(lambi) = 1 - att2(lambi); % calculates transmission
        end
  end
    
  for lambi = 1 : nl
        for D_n_ind = 1: nd
            r = R(3);
            D=lambda(lambi)*D_n(D_n_ind)/pi;
            sig_a = pi^2*D^3/lambda(lambi)*ImKm_w(lambi, ti); % calculates absorption cross section
            sig_s = 2*pi^5*D^6/3/lambda(lambi)^4*Km2_w(lambi, ti); % calculates scattering cross section
            sig_an = 4*sig_a/(pi*D^2); % calculates normalized absorption cross section
            sig_sn = 4*sig_s/(pi*D^2); % calculates normalized scattering cross section
            Nd = (0.08^-1)*exp(-41*r^(-0.21)*D); 
            sig_ed = sig_an + sig_sn; % calculates normalized extinction cross section
            att3(lambi) = att3(lambi) + (4.34e3 * Nd * sig_ed * deltaD); % calculates attenuation
            trans3(lambi) = 1 - att3(lambi); % calculates transmission
        end
  end
    
  for lambi = 1 : nl
        for D_n_ind = 1: nd
            r = R(4);
            D=lambda(lambi)*D_n(D_n_ind)/pi;
            sig_a = pi^2*D^3/lambda(lambi)*ImKm_w(lambi, ti); % calculates absorption cross section
            sig_s = 2*pi^5*D^6/3/lambda(lambi)^4*Km2_w(lambi, ti); % calculates scattering cross section
            sig_an = 4*sig_a/(pi*D^2); % calculates normalized absorption cross section
            sig_sn = 4*sig_s/(pi*D^2); % calculates normalized scattering cross section
            Nd = (0.08^-1)*exp(-41*r^(-0.21)*D); 
            sig_ed = sig_an + sig_sn; % calculates normalized extinction cross section
            att4(lambi) = att4(lambi) + (4.34e3 * Nd * sig_ed * deltaD); % calculates attenuation
            trans4(lambi) = 1 - att4(lambi); % calculates transmission
        end
  end
    
  for lambi = 1 : nl
        for D_n_ind = 1: nd
            r = R(5);
            D=lambda(lambi)*D_n(D_n_ind)/pi;
            sig_a = pi^2*D^3/lambda(lambi)*ImKm_w(lambi, ti); % calculates absorption cross section
            sig_s = 2*pi^5*D^6/3/lambda(lambi)^4*Km2_w(lambi, ti); % calculates scattering cross section
            sig_an = 4*sig_a/(pi*D^2); % calculates normalized absorption cross section
            sig_sn = 4*sig_s/(pi*D^2); % calculates normalized scattering cross section
            Nd = (0.08^-1)*exp(-41*r^(-0.21)*D); 
            sig_ed = sig_an + sig_sn; % calculates normalized extinction cross section
            att5(lambi) = att5(lambi) + (4.34e3 * Nd * sig_ed * deltaD); % calculates attenuation
            trans5(lambi) = 1 - att5(lambi); % calculates transmission
        end
  end
    
  for lambi = 1 : nl
        for D_n_ind = 1: nd
            r = R(6);
            D=lambda(lambi)*D_n(D_n_ind)/pi;
            sig_a = pi^2*D^3/lambda(lambi)*ImKm_w(lambi, ti); % calculates absorption cross section
            sig_s = 2*pi^5*D^6/3/lambda(lambi)^4*Km2_w(lambi, ti); % calculates scattering cross section
            sig_an = 4*sig_a/(pi*D^2); % calculates normalized absorption cross section
            sig_sn = 4*sig_s/(pi*D^2); % calculates normalized scattering cross section
            Nd = (0.08^-1)*exp(-41*r^(-0.21)*D); 
            sig_ed = sig_an + sig_sn; % calculates normalized extinction cross section
            att6(lambi) = att6(lambi) + (4.34e3 * Nd * sig_ed * deltaD); % calculates attenuation
            trans6(lambi) = 1 - att6(lambi); % calculates transmission
        end
    end
     

    
% Print Attenuation Graph
fntsz = 14;
figure(1)
clf
plot (freq/1e9, att1, 'b')
hold on
plot (freq/1e9, att2, 'r')
plot (freq/1e9, att3, 'k')
plot (freq/1e9, att4, 'g')
plot (freq/1e9, att5, 'r--')
plot (freq/1e9, att6, 'b--')
set(gca, 'fontsize', fntsz)
xlabel('Frequency (GHz)')
ylabel('Attenuation (dB / Km)')
title(ptit)
hold off
  
% Print Transmission Graph
fntsz = 14;
figure(2)
clf
plot (freq, trans1, 'b')
plot (freq, trans2, 'r')
plot (freq, trans3, 'k')
plot (freq, trans4, 'g')
plot (freq, trans5, 'r--')
plot (freq, trans6, 'b--')
set(gca, 'fontsize', fntsz)
xlabel('Frequency (Hz)')
ylabel ('Transmission (dB / Km)')
title(ptit)