%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble the full hydrodynamic system matrix 
% Taken from the appendix B of Kumar & Tuckerman (1994).
% 
% The only one of the equations which couples the different Fourier modes is the
% pressure jump condition (3.9)
%
% Δ[ρ{σ + i(ɑ + nω)} + 3μk²]∂_z(w_n) 
% - [Δμ ∂_zzz(w_n)] + (Δ[ρ]g - γk²) ζ_n = Δ[ρ] k² (1/2) a (ζ_{n+1} + ζ_{n-1})
%
% where : 
% - ρ is the density 
% - μ is the viscosity
% - k is the wavenumber
% - g is the gravity
% - γ is the surface tension coefficient
% - σ + i(ɑ + nω) are the Floquet exponents
% 
% Note that ∂_z(w_n) and ∂_zzz(w_n) are also functions of ζ_n through 
% the kinematic condition
%
%   w_{1n}(z=z1) = w_{2n}(z=z2) = [σ + i(ɑ + nω)] ζ_n
%
% the entire left-hand side of (3.9) can be written as
%
%   (1/2) Δ[ρ] k² A_n ζ_n = a (ζ_{n+1} + ζ_{n-1})
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dn] = matrix_An(sigma, alpha, n, omega, k, nu, mu, rho, g, gamma, h)
  % Compute the complex frequency term
  fexp = sigma + 1i * (alpha + n * omega);

  % Compute the wave numbers q1n and q2n
  q1n = sqrt(k.^2 + fexp / nu(1));
  q2n = sqrt(k.^2 + fexp / nu(2));

  % Define the z-coordinate for the interface surface
  interface_surface = 0;

  % Compute exponential terms for the wave number k at the interface surface
  exp_k_z1_plus = exp(k * interface_surface);
  exp_k_z1_minus = exp(-k * interface_surface);

  eps = 1e-3;
  if fexp ~= 0 % general case
    % Assemble Full Hydrodynamic System
    lhs = fhs_general(k, mu, h, q1n, q2n);

    if (rcond(lhs) > eps)
			rhs = zeros(8, 1);
      rhs(1) = fexp;

			% Solve Coefficients a1n, b1n, c1n, d1n, a2n, b2n, c2n, d2n, 
      coefficients = lhs \ rhs;

      % Compute derivatives
      exp_q1n_z1_plus = exp( q1n * interface_surface);
      exp_q2n_z1_plus = exp( q2n * interface_surface);
      exp_q1n_z1_minus = exp(-q1n * interface_surface);
      exp_q2n_z1_minus = exp(-q2n * interface_surface);
      exp_terms_z1_q1 = [exp_k_z1_plus, exp_k_z1_minus, exp_q1n_z1_plus, exp_q1n_z1_minus];
      exp_terms_z1_q2 = [exp_k_z1_plus, exp_k_z1_minus, exp_q2n_z1_plus, exp_q2n_z1_minus];
      dz1 = ([k, -k, q1n, -q1n] .* exp_terms_z1_q1) * coefficients(1:4);
      dz2 = ([k, -k, q2n, -q2n] .* exp_terms_z1_q2) * coefficients(5:8);
      dzzz1 = ([k^3, -k^3, q1n^3, -q1n^3] .* exp_terms_z1_q1) * coefficients(1:4);
      dzzz2 = ([k^3, -k^3, q2n^3, -q2n^3] .* exp_terms_z1_q2) * coefficients(5:8);
		else 
      lhs = fhs_general_deep(k, mu, q1n, q2n);

      % limit to infinite layer
      rhs = zeros(4, 1);
      rhs(1) = fexp;
			
      % Solve Coefficients a1n, c1n, b2n, and d2n
			coefficients = lhs \ rhs;
			
      % Compute derivatives
      dz1 = [ k,  q1n] * coefficients(1:2);
      dz2 = [-k, -q2n] * coefficients(3:4);
      dzzz1 = [ k^3, q1n^3,] * coefficients(1:2);
      dzzz2 = [-k^3, -q2n^3] * coefficients(3:4);
	  end

  else % special case
    % Assemble Full Hydrodynamic System for the special case
    lhs = fhs_special(k, mu, h);
    
    if (rcond(lhs) > eps)
			rhs = zeros(8, 1);
      rhs(1) = fexp;

			% Solve Coefficients a1n, b1n, c1n, d1n, a2n, b2n, c2n, d2n, 
      coefficients = lhs \ rhs;

      % Compute derivatives dw1/dz, and d³w1/dz³
      exp_terms_z1 = [exp_k_z1_plus, exp_k_z1_minus, exp_k_z1_plus, exp_k_z1_minus];
      dz1 = ([k, -k, (1 + k * interface_surface), (1 - k * interface_surface)] .* exp_terms_z1) * coefficients(1:4);
      dz2 = ([k, -k, (1 + k * interface_surface), (1 - k * interface_surface)] .* exp_terms_z1) * coefficients(5:8);
      dzzz1 = ([k^3, -k^3, 3 * (k^2) + (k^3) * interface_surface, 3 * (k^2) - (k^3) * interface_surface] .* exp_terms_z1) * coefficients(1:4);
      dzzz2 = ([k^3, -k^3, 3 * (k^2) + (k^3) * interface_surface, 3 * (k^2) - (k^3) * interface_surface] .* exp_terms_z1) * coefficients(5:8);

		else 
      lhs = fhs_special_deep(k, mu);

      % limit to infinite layer
      rhs = zeros(4, 1);
      rhs(1) = fexp;
			
      % Solve Coefficients a1n, c1n, b2n, and d2n
			coefficients = lhs \ rhs;

      % Compute derivatives dw1/dz, and d³w1/dz³
      dz1 = [ k, (1 + k * interface_surface)] * coefficients(1:2);
      dz2 = [-k, (1 - k * interface_surface)] * coefficients(3:4);
      dzzz1 = [ k^3, 3 * (k^2) + (k^3) * interface_surface] * coefficients(1:2);
      dzzz2 = [-k^3, 3 * (k^2) - (k^3) * interface_surface] * coefficients(3:4);
      
	  end

    
  end

  % Assemble Matrix
  Dn = ((rho(2) * fexp + 3 * mu(2) * k^2) * dz2 - mu(2) * dzzz2 - ...
        (rho(1) * fexp + 3 * mu(1) * k^2) * dz1 + mu(1) * dzzz1 + ...
        ((rho(2) - rho(1)) * g - gamma * k^2) * k^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
