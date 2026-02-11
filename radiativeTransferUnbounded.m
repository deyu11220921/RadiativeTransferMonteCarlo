function obs = radiativeTransferUnbounded( geometry, source, material, observation )

% physics
d = geometry.dimension;
acoustics = material.acoustics;

if ~isfield(material,'timeSteps')
    material.timeSteps = 0; % default value
end

if ~isfield(geometry,'frame')
    geometry.frame = 'spherical';
end

% discretization in packets of particles  (for optimal vectorization)
Npk = 1e5; % size of packets 
Np = ceil(source.numberParticles/Npk); % number of packets

% initialize observation structure
[ obs, E, bins, ibins, vals, Nt, t , d1, d2 ] = ...
         initializeObservation( geometry, acoustics, observation, Np*Npk );

% prepare scattering cross sections
material = prepareSigma( material, d );      

% loop on packages of particles
for ip = 1:Np

    % initialize particles
    P = initializeParticle( Npk, d, acoustics, source );

    % loop on time
    for it = 2:Nt

        % propagate particles
        if material.timeSteps==0
            P = propagateParticleSmallDt( material, geometry, P, t(it) );
        elseif material.timeSteps==1
            P = propagateParticle( material, P, t(it) );
        end

        % observe energies (as a function of [Psi x t])
        E(:,:,it,:) = E(:,:,it,:) + observeTime( geometry, acoustics, ...
                                      P.x, P.p, P.dir, bins, ibins, vals );

    % end of loop on time
    end

% end of loop on packages
end

% energy density as a function of [x1 x2 t]
obs.energyDensity = (1./(d1'*(d2*obs.N))).*double(E);

% energy as a function of [t]
%obs.energy = squeeze(tensorprod(d1,reshape(tensorprod(d2,obs.energyDensity,2,2),[length(d1) Nt 1+~acoustics]),2,1));
% --- MATLAB R2019b Compatibility Workaround Start ---

% Determine the number of time steps (Nt) from the observation structure
% Based on your provided image, this should be 61
Nt = length(observation.time); 

% Check if the simulation is for acoustics to determine the number of energy components
% This corresponds to the (1 + ~acoustics) logic in the original code
if exist('acoustics', 'var') && ~acoustics
    numComponents = 2;
else
    numComponents = 1;
end

% Step 1: Perform the inner contraction (equivalent to inner tensorprod)
% This weighs the energyDensity by the directional vector d2
% We reshape d2 to align with the 2nd dimension of obs.energyDensity
d2_weight = reshape(d2, [1, length(d2), 1]);
inner_result = sum(obs.energyDensity .* d2_weight, 2);

% Step 2: Reshape the result to match the expected dimensions [Spatial x Time x Components]
% The teacher mentioned counting particles in specific volumes/shells
intermediate_matrix = reshape(inner_result, [length(d1), Nt, numComponents]);

% Step 3: Perform the outer contraction (equivalent to outer tensorprod)
% This weighs the result by d1 across the spatial/basis dimension
% We use squeeze to remove singleton dimensions and produce the final energy curve
d1_weight = reshape(d1, [length(d1), 1, 1]);
obs.energy = squeeze(sum(intermediate_matrix .* d1_weight, 1));

% --- Compatibility Workaround End ---