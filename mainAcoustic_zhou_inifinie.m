

%% 1. Geometry
% Infinite medium (spherical frame)
geometry = struct( 'dimension', 3 , ...   
                   'frame', 'spherical' ); 
               
%% 2. Point Source
% Small lambda to avoid numerical singularity at origin
source = struct( 'numberParticles', 2e5 , ... 
                 'type', 'point', ...         
                 'position', [0 0 0], ...     
                 'direction', 'uniform',...%[0 0 1], ...    
                 'lambda', eps ); 

%% 3. factor
factor = 10;
L = 10; dr = 0.25; m = ceil(L/dr);
% Spatial bins (Radius R)
if rem(m,2)==0
    obs_radius_bins = L-((m-1):-2:-1)*dr;
else
    obs_radius_bins = L-((m-2):-2:-1)*dr;
end
%% 4. Material Properties
freq = 1e4;          
material.v = 450000; 
material.acoustics = true; 
material.sigma{1} = @(th) 1/4/pi*ones(size(th))*material.v*factor/L;

material = prepareSigma(material,geometry.dimension);
% Prepare scattering matrices
fprintf('Preparing scattering cross-sections...\n');
material = prepareSigma(material, geometry.dimension);
fprintf('Mean Free Path l_s = %.2f mm\n', material.meanFreePath);
% Fix: Z is elevation in [-pi/2, pi/2], not colatitude [0, pi]

%% 5. observation
tb = L/material.v;
dt = material.meanFreeTime;
obs_time = 0:dt:20*tb; 
observation = struct('x', obs_radius_bins, ...      
                     'y', [-pi pi], ...           % Azimuth (Full circle)
                     'z', [-pi/2 pi/2], ...       % Elevation (Full sphere)
                     'directions', [0 pi], ...    
                     'time', obs_time );     


     
%% 6. Solution
fprintf('Running Simulation (Infinite Medium)...\n');
tic ; obs = radiativeTransferUnbounded( geometry, source, material, observation ); toc;

% Paasschens 
%[EP, Ediff] = Comparison.Paasschens_RTE_Unbounded(source, material, observation, geometry);

% %% 6. Results Processing and Comparison
% % --- Step 1: Replicate Paasschens' Internal Grid Logic ---
% Lmax = material.v * max(observation.time);
% Rmax = Lmax; % For infinite medium
% r = linspace(0, Rmax, ceil(Rmax/observation.dr_paa));
% dr = mean(diff(r));
% 
% % --- Step 2: Extract Analytical Results at Target Distance ---
% target_dist = 10; 
% [~, paa_idx] = min(abs(r - target_dist));
% r_actual = r(paa_idx);
% 
% % step 2 Extract shell energy (Note: Paasschens returns energy integrated over the shell)
% y_paa_shell = EP(paa_idx, :); 
% y__paa_diff_shell = Ediff(paa_idx, :);
% 
% % --- Step 3: Density and Unit Conversion ---
% % 1. Convert Shell Energy to Energy Density (1/mm^3)
% % 3D Formula: Density = Energy / (4 * pi * r^2 * dr)
% V_shell = 4 * pi * r_actual^2 * dr;
% y_paa_density = y_paa_shell / V_shell;
% y_paa_diff_density = y__paa_diff_shell / V_shell;

% % 2. Apply the 1e9 factor to convert J/mm^3 to J/m^3 
% y_paa_final = y_paa_density * 1e9;
% y__paa_diff_final = y_paa_diff_density * 1e9;

% --- Step 4: Monte Carlo Processing ---
% Extract MC energy density and apply your 1e9 scaling
en_mc = squeeze(obs.energyDensity);
y_sim_mc = en_mc(end, :) * 1e9; 
x_sim = obs.t / tb;  
% --- Step 5: Intrinsic Attenuation (Qi) ---
% Apply the same exponential decay to both results
Qi = 500; f = 1e6; omega = 2 * pi * f; 
att = exp(- (omega * observation.time) / Qi);

y_sim_mc = y_sim_mc .* att;
% y_paa_final = y_paa_final .* att;
% y__paa_diff_final = y__paa_diff_final .* att;

% --- Step 6: Amplitude Scaling (Normalization) ---
% [~, t_ballistic_idx] = min(abs(observation.time - r_actual/material.v));
% scale_factor = y_sim_mc(t_ballistic_idx) / y_paa_final(t_ballistic_idx);

% %% 7. Plotting
% --- 1. Plot the Monte Carlo result ---
h1 = semilogy(x_sim, y_sim_mc, 'k-', 'LineWidth', 2); 
hold on; % Keep the current plot to add more data

% --- 2. Define and plot the second dataset (h2) ---
% Extract x from the 1st column and y from the 2nd column
h2_10_data = [
    1.1958041958041954, 254043.43448706894;  
    2.3379953379953378, 44072.915332119934;
    3.3822843822843813, 14643.475730989523;  
    4.002331002331002,  7865.123048451399;  
    5.568764568764569,  2540.434344870684; 
    7.6573426573426575, 584.6166306860788;  
    9.386946386946388,  194.2421865872216;  
    11.051282051282051, 78.65123048451399; 
];

% Plot h2 as red circles or dashed line (adjust style as needed)
h2 = semilogy(h2_10_data(:,1), h2_10_data(:,2), 'ro--', 'LineWidth', 1.5, 'MarkerSize', 8);

% --- 3. Figure Decorations ---
grid on;
xlabel('Normalized Time (t / t_b)');
ylabel('Energy Density (J/m^3)');
title(['Infinite Medium (L = 10mm, l_s \approx 10mm)']);
% Add the ballistic arrival reference line
xline(1, 'r--', 'Ballistic Arrival', 'LabelVerticalAlignment', 'bottom');
% Set axes limits
xlim([0 20]);
% Add a legend to distinguish the two lines
legend([h1, h2], {'Monte Carlo (y\_sim\_mc)', 'Reference Data (h2)'}, 'Location', 'northeast');
hold off; 
% Release the plot
% figure;
% semilogy(observation.time/tb, y_sim_mc, '-k', 'LineWidth', 1.5, 'DisplayName', 'Monte Carlo'); 
% hold on; grid on;
% semilogy(observation.time/tb, y_paa_final * scale_factor, '--r', 'LineWidth', 1, 'DisplayName', 'Paasschens (Aligned)');
% semilogy(observation.time/tb, y__paa_diff_final * scale_factor, ':g', 'DisplayName', 'Diffusion Approx');
% 
% xlabel('Normalized Time (t / t_b)');
% ylabel('Energy Density (J/m^3)');
% title(['Comparison at r = ' num26str(r_actual) ' mm (l_s = 1 mm)']);
% legend('Location', 'northeast');
% xlim([0 20]);