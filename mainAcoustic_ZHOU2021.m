
geometry = struct( 'dimension', 3 , ...   
                   'frame', 'spherical' );
%dimension (integer 2 or 3): dimensionality of the propagation space.
%frame (chain of characters spherical (default) or cartesian): frame in which coordinates are expressed on output.
%bnd (structured array(s)): description of the boundaries through their normal (field dir, with values 1=x, 2=y, 3=z) and position (indicated by field pos).
% boundaries with normal 'dir' (1='x', 2='y', 3='z', 4='r_cylindrical') and position 'val'
% For dir=4, the symmetry axis of the cylinder is by default 'z' 
%geometry.bnd(1) = struct('dir',3,'val',0);

% Point source
source = struct( 'numberParticles', 1e4 , ...
                 'type', 'point', ...         % 'point' (default) or 'plane'
                 'position', [0 0 0], ...    % always in cartesian frame
                 'direction',[0 0 1],     ... % with point source: 'uniform' (default) or 'outgoing'; with plane source: 1='x', 2='y', 3='z' [x y z]
                 'radial', 0, ...             % only used with point sources. Initial condition function of radius (always spherical coordinates)
                 'extent', [10 10], ...       % only used with plane sources. Extent of the source around 'position' (always cartesian coordinates)
                 'lambda', 1e-4 );    

%% 3. Observation (??)
% ????: 0.5mm ? 20.5mm
obs_radius_bins = (0.5 : 1.0 : 20.5) * 1e-3; 

% ??: 0 ? 500 us (?? 22.2us ?????)
obs_time = 0 : 0.5e-6 : 500e-6;               % <--- ?????????

observation = struct('x', obs_radius_bins, ...      
                     'y', [0 2*pi], ...           % Phi ??
                     'z', [0 2*pi], ...             % Theta ??
                     'directions', [0 pi], ...    
                     'time', obs_time );          

%% 4. Material (??????? calibration)
freq = 1e6;       
V_paper = 450;    
variance = 0.08;       
corr_length = 6.25e-4;  

material = MaterialClass( geometry, ...
                          freq, ...
                          true, ...          
                          V_paper, ...             
                          [variance variance], ...     
                          -0.5, ...          
                          'exp', ...         
                          corr_length);            
                          
%% 5. Solution
fprintf('Preparing scattering cross-sections (Sigma)...\n');
material = prepareSigma(material, geometry.dimension);
fprintf('Running Simulation (Infinite Medium)...\n');
% ?? Unbounded ??
obs = radiativeTransferUnbounded( geometry, source, material, observation );

%% 6. Plotting (?? Figure 1)
target_L = 10e-3; % ???? 10mm

% ??????? bin
bin_centers = (obs_radius_bins(1:end-1) + obs_radius_bins(2:end)) / 2;
[~, idx_L] = min(abs(bin_centers - target_L)); 
fprintf('Extracting energy at R = %.2f mm\n', bin_centers(idx_L)*1000);

% ???? (??????)
E_t = squeeze(obs.energyDensity(idx_L, :, :, :));
E_t = E_t(:);

% ??????
tb = target_L / V_paper;
t_normalized = obs_time / tb;

% ???????????
min_len = min(length(E_t), length(t_normalized));
E_t = E_t(1:min_len);
t_normalized = t_normalized(1:min_len);
% ... ? E_t = E_t(1:min_len); ??? semilogy ???? ...

% --- DEBUG START ---
fprintf('Data Check:\n');
fprintf('  Max Energy: %e\n', max(E_t));
fprintf('  Min Energy: %e\n', min(E_t));
if max(E_t) == 0
    warning('??????????? 0???? material ??? source ???');
    return; % ???????????
end
% --- DEBUG END ---

figure('Name', 'Reproduction Figure 1');

% ??
figure('Name', 'Reproduction Figure 1');
semilogy(t_normalized, E_t, 'k-', 'LineWidth', 2); % ?????(k-)??????
grid on;
xlabel('Normalized Time (t / t_b)');
ylabel('Energy Density (J/m^3)');
title(['Infinite Medium (L = 10mm), t_b = ' num2str(tb*1e6, '%.1f') ' \mus']);
xline(1, 'r--', 'Ballistic Arrival'); % ??????
xlim([0 20]); % ??X???????



