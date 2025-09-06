%==========================================================================
% MATLAB SCRIPT TO SOLVE FOR TORPEDO ANCHOR PENETRATION DEPTH IN SAND
% according to the formulation by Li et. al (2022)
% which is based on the cavity expansion theory
%==========================================================================

clear; clc; close all;

% --- Input Parameters ---

% Anchor Properties
a = 0.0125;       % Radius of the anchor shaft (m)
L = 0.175;        % Length of the anchor (m)
m = 0.235;        % Mass of the anchor (kg)

% Lenghth of the anchor tip = (2a)/(2tan(15 deg))  
Lc = (2*a)/(2*tand(15));

% Radius of curvature of the anchor tip arc = (Lc^2 + a^2)/(2a)
s = (Lc^2 + a^2)/(2*a);     

% Physical Constants
g = 9.81;       % Acceleration due to gravity (m/s^2)

% Soil properties

rho_0 = 1834;     % Initial bulk density of soil (kg/m^3)
gamma_soil = rho_0 * g; % Unit weight of soil (N/m^3)

rho_star = 1930;  % Locked (compacted) density of soil (kg/m^3)
                      % Used as maximum dry density
                  
E = 100e6;       % Young's modulus of the soil (Pa)


%%%%%%%%%%%%%%%%%%%%%
tau_0 = 31549;     % Soil cohesion or yield shear stress (Pa)
lambda = 1.47;    % Yield parameter
%%%%%%%%%%%%%%%%%%%%%

eta_1 = 0.30;      % Friction coefficient (anchor tip-soil interface)

K = 2.74;         % Earth pressure coefficient

delta = 19.29;    % Interface friction angle (degrees)
phi = 27;         % Internal friction angle of soil (degrees)


% Initial Conditions
V_0 = 22;       % Initial velocity of the anchor at seabed (m/s)

% --- Calculations ---

fprintf('Calculating intermediate parameters...\n');

% 1. Geometric and material ratios
Psi = s / (2*a);
phi_0 = asin((2*Psi - 1) / (2*Psi));
eta_2 = tan(deg2rad(delta)) / tan(deg2rad(phi));
eta_star = 1 - rho_0 / rho_star;
alpha_param = 3 * lambda / (3 + 2 * lambda);


%%
% Cavity Expansion Parameters
%  V is the cavity expansion velocity (m/s)
%  c is the elastic-plastic interface velocity (m/s)

% gamma_ratio = V / c; 
gamma_ratio = 0.385;

%%

% 2. Cavity expansion parameters A, B, C1, C2
C1 = (3 * tau_0 / E) + eta_star * (1 - (3 * tau_0) / (2 * E))^2;
C2_num = -gamma_ratio^3 * (2 * (1 - eta_star) * (2 - alpha_param) + 3 * gamma_ratio^3);
C2_den = (1 - eta_star) * (1 - 2 * alpha_param) * (2 - alpha_param) * (1 + tau_0 / (2*E))^4;
C2 = C2_num / C2_den;

A = (1/alpha_param) * ( (1 + tau_0 / (2*E)) / gamma_ratio ) - (1/lambda);
B_term1 = 3 / ((1 - eta_star) * (1 - 2 * alpha_param) * (2 - alpha_param));
B_term2 = (1 / gamma_ratio^2) * (( (1 + tau_0 / (2*E)) / gamma_ratio )^(2 * alpha_param));
B = B_term1 + B_term2 * (C1 + C2);

% 3. Static and dynamic resistance terms (alpha_s, beta_s)
term1_alpha = 1 + 4 * eta_1 * Psi^2 * (pi/2 - phi_0);
term2_alpha = eta_1 * (2*Psi - 1) * sqrt(4*Psi - 1);
alpha_s = pi * a^2 * tau_0 * A * (term1_alpha - term2_alpha);

term1_beta = (8*Psi - 1) / (24 * Psi^2);
term2_beta = eta_1 * Psi^2 * (pi/2 - phi_0);
term3_beta = (eta_1 * (2*Psi - 1) * (6*Psi^2 + 4*Psi - 1) * sqrt(4*Psi - 1)) / (24*Psi^2);
beta_s = pi * a^2 * rho_0 * B * (term1_beta + term2_beta - term3_beta);

% 4. Coefficients for the partial embedment equation (a1, a2, a3)
a2_factor = (m / (2 * beta_s))^3;
a2 = (2 * pi * a * eta_2 * K * gamma_soil / m) * a2_factor;

a1_term1 = 2 * a2; 
a1_term2 = -(m*g - alpha_s) / beta_s;
a1 = a1_term1 + a1_term2 + V_0^2;

a3 = ((m*g - alpha_s) / beta_s) - 2*a2;

% --- Solve for Penetration Depth (Z_max) ---
% DESCRIPTION: The script first attempts to solve the equation for partial
% embedment. If no valid solution is found, it proceeds to solve the
% equation for full embedment.

fprintf('Solving for maximum penetration depth (Z_max)...\n');
syms Z_max

% Equation for partial embedment (Z_max <= L)
x = (2 * beta_s * Z_max) / m;
eqn1 = a1 * exp(-x) - a2 * x^2 + 2 * a2 * x + a3 == 0;

% Use vpasolve to find a numerical solution within the valid range [0, L]
sol_Z1 = vpasolve(eqn1, Z_max, [0, L]);

if ~isempty(sol_Z1)
    % A solution was found for partial embedment
    Z_final = double(sol_Z1);
    fprintf('\n--- SOLUTION FOUND ---\n');
    fprintf('Penetration corresponds to PARTIAL EMBEDMENT.\n');
    fprintf('Final Penetration Depth (Z_max): %.4f m\n', Z_final);
    is_full_embedment = false;
else
    % No solution for partial embedment, so it must be full embedment
    fprintf('No partial embedment solution found. Checking for full embedment...\n');

    % 5. Calculate R, W, and a4 for the full embedment case
    R = 2*g - (2*alpha_s/m) + (2*pi*a*K*eta_2*gamma_soil*L^2 / m);
    W = (4*pi*a*K*eta_2*gamma_soil*L) / m;
    
    a4_term1 = exp((2*beta_s*L)/m);
    a4_term2 = a1 * exp((-2*beta_s*L)/m);
    a4_term3 = -a2 * ((2*beta_s*L)/m)^2;
    a4_term4 = 2 * a2 * ((2*beta_s*L)/m);
    a4_term5 = W * (m/(2*beta_s))^2 * (((2*beta_s*L)/m) - 1);
    a4_term6 = -(m*R)/(2*beta_s);
    a4 = a4_term1 * (a4_term2 + a4_term3 + a4_term4 + a3 + a4_term5 + a4_term6);
    
    % Equation for full embedment (Z_max > L)
    eqn2 = (2 * a4 * beta_s / m) * exp(-(2 * beta_s * Z_max) / m) - W * Z_max + R + (W * m) / (2 * beta_s) == 0;
    
    % Solve the second equation for depths greater than L
    sol_Z2 = vpasolve(eqn2, Z_max, [L, Inf]);
    
    if ~isempty(sol_Z2)
        Z_final = double(sol_Z2);
        fprintf('\n--- SOLUTION FOUND ---\n');
        fprintf('Penetration corresponds to FULL EMBEDMENT.\n');
        fprintf('Final Penetration Depth (Z_max): %.4f m\n', Z_final);
        is_full_embedment = true;
    else
        Z_final = NaN;
        fprintf('\n--- ERROR ---\n');
        fprintf('Could not find a valid solution for either partial or full embedment.\n');
        is_full_embedment = false;
    end
end

% --- Calculate and Plot Instantaneous Velocity ---
% DESCRIPTION: If a valid solution for Z_max was found, this section
% calculates the velocity at each point and generates a plot.

if ~isnan(Z_final)
    fprintf('\nGenerating velocity profile plot...\n');
    z_vals = linspace(0, Z_final, 500);
    Vz_sq = zeros(size(z_vals)); % Vz squared

    for i = 1:length(z_vals)
        z = z_vals(i);
        x_z = (2 * beta_s * z) / m;
        
        if z <= L
            % Partial embedment velocity equation
            Vz_sq(i) = a1 * exp(-x_z) - a2 * x_z^2 + 2 * a2 * x_z + a3;
        else
            % Full embedment velocity equation
            % Note: a4, R, and W must be calculated as they are only defined
            % for the full embedment case.
            R = 2*g - (2*alpha_s/m) + (2*pi*a*K*eta_2*gamma_soil*L^2 / m);
            W = (4*pi*a*K*eta_2*gamma_soil*L) / m;
            a4_term1 = exp((2*beta_s*L)/m);
            a4_term2 = a1 * exp((-2*beta_s*L)/m);
            a4_term3 = -a2 * ((2*beta_s*L)/m)^2;
            a4_term4 = 2 * a2 * ((2*beta_s*L)/m);
            a4_term5 = W * (m/(2*beta_s))^2 * (((2*beta_s*L)/m) - 1);
            a4_term6 = -(m*R)/(2*beta_s);
            a4 = a4_term1 * (a4_term2 + a4_term3 + a4_term4 + a3 + a4_term5 + a4_term6);
            
            Vz_sq(i) = a4 * exp(-x_z) + (m*R)/(2*beta_s) - W * (m/(2*beta_s))^2 * (x_z - 1);
        end
    end
    
    % Ensure Vz^2 is not negative due to numerical approximations near zero
    Vz_sq(Vz_sq < 0) = 0;
    Vz = sqrt(Vz_sq);
    
    % Plotting
    figure('Name', 'Anchor Velocity Profile');
    plot(Vz, 100*z_vals, 'LineWidth', 2, 'Color', [0, 0, 0]);  
    set(gca, 'YDir','reverse'); % Set depth (z) to increase downwards
    title('Instantaneous Velocity vs. Penetration Depth');
    xlabel('Instantaneous Velocity, V_z (m/s)');
    ylabel('Penetration Depth, z (cm)');
    axis tight;
    
    % Add a line indicating the anchor length L if full embedment occurred
    if is_full_embedment
        hold on;
        yline(100*L, '--r', {'Anchor Length (L)'}, 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
        hold off;
    end
end



