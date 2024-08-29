% Main code to calculate eigenvalues for the parametric instability of the
% interface between two fluids using the approach to Kumar & Tuckerman (1994)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_dim = 64*pi*1e-3;                           % width, m
g_dim = 9.81;                                 % gravity, m/s²
rho_dim = [519.933, 415.667];                 % densities (bottom/top), kg/m³
mu_dim = [3.908e-5, 3.124e-5];                % dyn. visc. Pa s
nu_dim = mu_dim./rho_dim;                     % kin. visc. m²/s
h_dim = [1000, 1000];                         % layer heights, m
gamma_dim = 2.181e-6;                         % interfacial tension, N/m²
f_dim = 100;                                  % frequency, Hz
omega_dim = 2*pi*f_dim;                       % frequency, rad/s
m_max = 16000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 20;                                       % Trunctation Fourier series
atwood = -diff(rho_dim)/sum(rho_dim);         % Atwood number


k0 = [1:m_max] * pi / w_dim;
omegan = sqrt(atwood*g_dim*k0);               % Nat. frequency

figure(1)
plot(k0, omegan)
hold on 
plot(k0,k0*0+omega_dim)
hold off


[a_harmonic, a_subharmonic] = sweep_over_k(w_dim, N, omega_dim, nu_dim, mu_dim, rho_dim, g_dim, gamma_dim, h_dim, m_max);


figure(2)
k = a_harmonic(:,1);
plot(k, a_harmonic(:,3), '.k')
hold on 
plot(k, a_subharmonic(:,3), '.r')
for i = 4:12
  plot(k, a_harmonic(:,i), '.k')
  plot(k, a_subharmonic(:,i), '.r')
end
hold off
xlabel('$k$', 'Interpreter', 'latex')
ylabel('$a/g$', 'Interpreter', 'latex')
xlim([0,1.5e5])

% Set figure properties
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4 3]); % Set the figure size (width x height)
set(gcf, 'PaperSize', [4 3]); % Ensure the paper size matches the position

% Set font size and line width
set(gcf, 'DefaultAxesFontSize', 12);
set(gcf, 'DefaultTextFontSize', 12);
set(gca, 'LineWidth', 1.5);

% Set the resolution (DPI)
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', '-r300', 'KT94_fig1b.eps'); % Save as EPS with 300 DPI

save('KT94_fig1b.mat', "a_harmonic", "a_subharmonic")