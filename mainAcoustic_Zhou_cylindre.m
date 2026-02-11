%% Zhou 2021 Figure 4 Reproduction - Finite Cylinder Case
% Objective: Reproduce the Red Curve (Cylinder) from Figure 4
% 1. Frame: Cylindrical (r, phi, z)
% 2. Boundaries: Lateral walls + Top + Bottom (Total Reflection)
% 3. Observation: Z-Slices + Single Large R-Disk (to smooth noise)

clear; clc; close all;

%% 1. Geometry (Cylindrical Finite Medium)
% Paper parameters: L = 10mm, w = 2.5mm (L/4) -> Radius r = 1.25mm
L_sample = 10e-3; 
r_sample = 1.25e-3; 

% Define geometry structure using cylindrical frame
geometry = struct('dimension', 3, ...
                  'frame', 'cylindrical'); 

% Define Boundaries - Perfect Reflection (R ~ 1)
% dir=4 -> Cylinder lateral wall, located at r = r_sample
geometry.bnd(1) = struct('dir', 4, 'val', r_sample); 
% dir=3 -> Top Z-axis boundary, located at z = L_sample
geometry.bnd(2) = struct('dir', 3, 'val', L_sample); 
% dir=3 -> Bottom Z-axis boundary, located at z = 0 (Ground)
geometry.bnd(3) = struct('dir', 3, 'val', 0);       

%% 2. Point Source
% Initial direction set to [0 0 1] (Z-axis) to avoid backward scattering at t=0
source = struct('numberParticles', 1e6, ... % Number of particles (Increase to 1e6 for smoother curves)
                'type', 'point', ...
                'position', [0 0 0], ...    % Source position at origin
                'direction', [0 0 1], ...   % Direction pointing into the medium
                'radial', 0, ...
                'lambda', 1e-4); 

%% 3. Observation Settings
% Define spatial bins to capture energy at specific locations

% Z-axis (Axial): Slices along the length
% Target: Extract energy at L=10mm. We define slices of thickness dz.
dz = 0.5e-3; % Slice thickness
obs_z_bins = 0 : dz : (L_sample + dz); 

% X-axis (Radius): 
% Single large bin covering the whole cross-section [0, r_sample] to reduce noise
obs_r_bins = [0 r_sample];

% Y-axis (Phi/Azimuth): 
% Full angular integration [-pi, pi]
obs_phi = [-pi pi];

% Time: Simulate up to 20 ballistic times (Refer to Figure 4)
V_wave = 450; 
tb = L_sample / V_wave; % Ballistic time (L/V)
obs_time = 0 : (tb/50) : (20 * tb); 

observation = struct('x', obs_r_bins, ...    % R (Radius)
                     'y', obs_phi, ...       % Phi (Azimuth)
                     'z', obs_z_bins, ...    % Z (Axial Slices)
                     'directions', [0 pi], ...
                     'time', obs_time);

%% 4. Material Properties 
freq = 1e6;       
acoustics = true;  % Acoustics enabled
corr_func = 'exp'; % Exponential correlation function

% =======================================================
% Calibration Values (From calcul_mean_free_path_ls.m)
% Replace these values based on which curve you want (L/ls = 1, 3, or 10) 1(0.123)   3(0.2153)  10(0.39)
calibrated_variance = 0.123;   % Epsilon (Variance)
calibrated_corr_len = 0.0003;   % Correlation Length a (meters) 0.0003
% =======================================================

material = MaterialClass(geometry, freq, acoustics, V_wave, ...
                         [calibrated_variance calibrated_variance], ...
                         -0.5, corr_func, calibrated_corr_len);
material.rho = 2300;
fprintf('Preparing Scattering Matrix...\n');
material = prepareSigma(material, geometry.dimension);
fprintf('------------------------------------------------\n');
fprintf('Current Mean Free Path l_s = %.3f mm\n', material.meanFreePath*1000);
fprintf('Ratio L/l_s = %.2f\n', L_sample / material.meanFreePath);
fprintf('------------------------------------------------\n');

%% 5. Simulation
fprintf('Running Simulation (Finite Cylinder)...\n');
% Call the main radiative transfer function (Unbounded handles boundaries via geometry struct)
obs = radiativeTransferUnbounded(geometry, source, material, observation);

%% 6. Data Processing & Plotting (Corrected for Finite Cylinder)

% 1. Extract Raw Data at the Transmission Boundary (Z = L)
% -------------------------------------------------------------------------
% obs.energyDensity dimensions are likely: [1 x Nz x Nt]
% We need to find the Z-bin index that corresponds to the output face (L_sample).

% Find the index in obs.z (bin centers) closest to L_sample (10mm)
% obs.z are the centers of the bins defined by 0:dz:(L+dz)
[~, idx_z] = min(abs(obs.z - L_sample)); 

% Safety check: Ensure we are not picking a bin strictly outside if it's empty
% If the bin center is > L_sample, check if it has data. If not, step back one bin.
max_val_at_z = max(obs.energyDensity(1, idx_z, :));
if max_val_at_z == 0 && idx_z > 1
    fprintf('Note: Z-bin at L_sample was empty, moving 1 bin inside.\n');
    idx_z = idx_z - 1;
end

fprintf('Extracting Energy at Z = %.3f mm (Index %d)\n', obs.z(idx_z)*1000, idx_z);

% Squeeze: [1 x 1 x Nt] -> [Nt x 1] or [1 x Nt]
% We select: Radius Bin=1 (All), Z Bin=idx_z, Time=All
E_raw = squeeze(obs.energyDensity(1, idx_z, :)); 

% Ensure E_raw is a row vector [1 x Nt]
E_raw = E_raw(:)'; 

% 2. Introduce Intrinsic Absorption Qi 
% -------------------------------------------------------------------------
Qi = 500;            
f = 1e6;             
omega = 2 * pi * f; 

% Calculate attenuation factor: exp(-omega * t / Qi)
if isfield(obs, 't')
    time_vec = obs.t;
else
    time_vec = obs.time;
end
absorption_factor = exp(- (omega * time_vec) / Qi);

% Apply absorption to the raw energy
E_final = E_raw .* absorption_factor; 

% 3. Normalize Time Axis
% -------------------------------------------------------------------------
t_norm = time_vec / tb;

% 4. Plotting
% -------------------------------------------------------------------------
figure('Color', 'w', 'Name', 'Zhou Fig 4 Reproduction (Red Curve)');

% IMPORTANT: Do NOT use smoothdata(..., 50) or 80. It kills the peak.
% We use raw data or very minimal smoothing (e.g., 3-5 points).
semilogy(t_norm, E_final, 'r.-', 'LineWidth', 1.5, 'MarkerSize', 8, ...
         'DisplayName', sprintf('Cylinder L/l_s=%.1f', L_sample/material.meanFreePath));

grid on;
xlabel('Normalized Time (t / t_b)', 'FontSize', 12);
ylabel('Energy Density (J/m^3)', 'FontSize', 12);
title({sprintf('Finite Cylinder (L/l_s \\approx %.1f, Q_i=%d)', L_sample/material.meanFreePath, Qi); ...
       'Expectation: Sharp Peak at t=1, then Reverberations'}, 'FontSize', 12);

% Add reference line for Ballistic Arrival
xline(1, 'k--', 'Ballistic Arrival');

% Adjust Axes
xlim([0 20]); 
% Automatically set Y-limits based on peak data
peak_E = max(E_final);
if peak_E > 0
    ylim([peak_E * 1e-4, peak_E * 2]); 
end

legend('show');