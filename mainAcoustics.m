
% This code works in 2D/3D. The initial condition is modeled as gaussian 
% with standard deviation source.lambda and uniformly-distributed angle(s).
% The initial direction is uniform.

% geometry
geometry = struct( 'dimension', 3 , ...   
                   'frame', 'cartesian' );
%dimension (integer 2 or 3): dimensionality of the propagation space.
%frame (chain of characters spherical (default) or cartesian): frame in which coordinates are expressed on output.
%bnd (structured array(s)): description of the boundaries through their normal (field dir, with values 1=x, 2=y, 3=z) and position (indicated by field pos).
% boundaries with normal 'dir' (1='x', 2='y', 3='z', 4='r_cylindrical') and position 'val'
% For dir=4, the symmetry axis of the cylinder is by default 'z' 
geometry.bnd(1) = struct('dir',3,'val',100);

% Point source
source = struct( 'numberParticles', 1e3 , ...
                 'type', 'point', ...         % 'point' (default) or 'plane'
                 'position', [0 0 0], ...    % always in cartesian frame
                 'direction',[1 0 0],     ... % with point source: 'uniform' (default) or 'outgoing'; with plane source: 1='x', 2='y', 3='z' [x y z]
                 'radial', 0.1, ...             % only used with point sources. Initial condition function of radius (always spherical coordinates)
                 'extent', [10 10], ...       % only used with plane sources. Extent of the source around 'position' (always cartesian coordinates)
                 'lambda', 0.1 );    

% observations
% choose 2 variables only to perform histograms (among x, y, z, directions)
% if geometry.frame = 'spherical', (x,y,z) correspond respectively to
%  x=r, y=azimuth (in [-pi pi]), z=elevation (in [-pi/2 pi/2])
% if geometry.frame = 'cylindrical', (x,y,z) correspond respectively to
%  x=r, y=azimuth (in [-pi pi]), z
observation = struct('x', -5:.1:5, ...                  % bins in space
                     'y', linspace(-5,5,100), ...                 
                     'z', [-Inf Inf], ...                 % unused in 2D
                     'directions', [0 pi], ...          % bins for directions [0 pi]         
                     'time', 0:0.05:5 );                % observation times

% material properties - more examples can be found in Example folder
% here is a basic example
freq = 1e6; % in Hz
material = MaterialClass( geometry, ...
                          freq, ...
                          true, ...          % true for acoustics
                          450, ...             % average wave velocity
                          [0.08 0.08], ...     % coefficients of variation of kappa and rho.
                          -0.5, ...          % correlation coefficient of kappa/rho
                          'exp', ...         % autocorrelation function
                           6.25e-4);            % correlation length
                          
% radiative transfer solution - acoustic with boundaries
obs = radiativeTransferUnbounded( geometry, source, material, observation );

% plotting output
sensors = [0 0 0; 
           0 0 5];

plotting = struct( 'checkEnergy', false, ...
                   'movieTotalEnergy', true, ...
                   'timehistory', false, ...
                   'sensors', sensors );

plotEnergies( plotting, obs, material, source.lambda )