clear
%% %% Parameter definition %%%%%

%%% Input files
population_file = 'Population_AB_Ext_3km.mat';
Vs30_file = 'Vs30_AB_Ext_3km_lit.mat';


%%% Parameters for grids
% !! Note that source grid nodes are offset from impact nodes to avoid
% instances where distance between source and impact node are zero !!

% EQ source grid
source_x0 = -121.55; % Minimum longitude
source_y0 = 48.05; % Minimum latitude
source_xfin = -108.5; % Maximum longitude
source_yfin = 61; % Maximum latitude
source_dx = 0.1; % Grid spacing (degrees longitude)
source_dy = 0.1; % Grid spacing (degrees latitude)

% Impact grid
% Note: impact grid extent and spacing is defined by input Population grid file


%%% Simulation parameters

Mag_rng = (2 : 0.5 : 5); % Magnitudes for simulation (minimum : magnitude step : maximum)
Nr = 80; % Number of realizations

% Trailing seismicity:
trailing_seismicity = true; % true or false
% Parameters for calculation of trailing seismicity
bm=1; % b-value from Gutenberg-Richter distribution
db=0.05; % For perturbation of b-value
Rf=1.2; % Ratio between earthquakes occurring during the stimulation and trailing seismicity periods

% Depth Variation
% Elevation convention: units are km with positive values below surface
z_d=[3 5 7]; % Possible Z coordinates for earthquake
wz=[30 50 20]; % Relative probability of each depth as a percentage

% Extent of consideration for nuisance and damage calculations
nuisance_distance_limit = 400; % Nuisance distance limit in km
damage_distance_limit = 40; % Damage distance limit in km

% Specify whether random state is reset to a defined state for each EQ
% source node.
random_type = 'shuffle'; % must be 'consistent' or 'shuffle'

% Check for expected type:
if ~or(strcmp(random_type, 'consistent'), strcmp(random_type, 'shuffle'))
    error("Random_type must be defined as 'consistent' or 'shuffle'")
end


% Simulation file save location
save_location = 'simulations.mat';


%% %% Import data grids %%%%%
Pop = importdata(population_file); %%% read Population data %%%
disp("Reading " + Vs30_file)
Vs30 = importdata(Vs30_file); %%% read Vs30 data %%%
disp(Vs30_file + " read")


%% %% Create the structure array grid %%%%%

%%% Create a source grid

grid.source.dx = source_dx;
grid.source.dy = source_dy;
grid.source.x0 = source_x0;
grid.source.y0 = source_y0;
grid.source.xfin = source_xfin;
grid.source.yfin = source_yfin;

grid.source.xcoords = [grid.source.x0 : grid.source.dx : grid.source.xfin];
grid.source.ycoords = [grid.source.y0 : grid.source.dy : grid.source.yfin];

% Initialize a cell array
grid.source.coordinates = {};

% Populate the cell array
for i = length(grid.source.xcoords) : -1 : 1
    for j = length(grid.source.ycoords) : -1 : 1
        grid.source.coordinates{i,j} = [grid.source.xcoords(i), grid.source.ycoords(j)];
    end
end

% Note Magnitudes in source grid
grid.source.Magnitudes = Mag_rng;


%%% Define the impact grids and their properties
% Note: impact grid extent and spacing is defined by input Population grid
% file, which matches the Vs30 grid.

% X, Y coordinates for Impact cells
x=Pop(:,1);y=Pop(:,2);

% Pop, Vs30 (mean) for Impact cells
Pop2=Pop(:,3); Vs30_2=Vs30(:,5);

% Create the grids
grid.impact.xcoords = sort(unique(x));
grid.impact.ycoords = sort(unique(y));

grid.impact.coordinates = {}; % Initialize coordinates cell array
%grid.impact.z = []; % Initialize depth array

% Populate the impact coordinates cell array
for i = length(grid.impact.xcoords) : -1 : 1
    for j = length(grid.impact.ycoords) : -1 : 1
        grid.impact.coordinates{i,j} = [grid.impact.xcoords(i), grid.impact.ycoords(j)];
    end
end

% Populate impact grids
grid.impact.Pop = [];
grid.impact.Vs30 = [];
grid.impact.Vs30_simulations = {};

for i = 1:length(Pop)
    impacti = find(grid.impact.xcoords == Pop(i,1));
    impactj = find(grid.impact.ycoords == Pop(i,2));
    grid.impact.Pop(impacti, impactj) = Pop(i,3);
end
for i = 1:length(Vs30)
    impacti = find(grid.impact.xcoords == Vs30(i,1));
    impactj = find(grid.impact.ycoords == Vs30(i,2));
    grid.impact.Vs30(impacti, impactj) = Vs30(i,5);
    grid.impact.Vs30_simulations{impacti, impactj} = sort(Vs30(i,6:end));
end

% Populate empty Vs30 and Vs30_simulation cells with mean of Vs30
grid.impact.Vs30(grid.impact.Vs30==0) = mean(Vs30(:,3));

for i = length(grid.impact.xcoords) : -1 : 1
    for j = length(grid.impact.ycoords) : -1 : 1
        if length(grid.impact.Vs30_simulations{i,j}) == 0
            grid.impact.Vs30_simulations{i,j} = zeros(size(Vs30(1,6:end))) + mean(Vs30(:,3));
        end
    end
end


%%% Pre-calculate epicentral distances between source and impact cells.
%%% Create a cellular array on the source grid where each cell
%%% contains an array of the epicentral distances to all the impact cells.

% Initiate cell array to hold distance matricies
grid.source.epi_distances = grid.source.coordinates;

% Transform coordinate grid into vectors
temp = cell2mat(reshape(grid.impact.coordinates, [], 1));
temp_x = temp(:,1);
temp_y = temp(:,2);

% This process can take some time. Show progress with tracking bar.
f = waitbar(0, 'Starting Distance Calculation');

for i = 1:length(grid.source.xcoords)
    for j = 1:length(grid.source.ycoords)
        % Calculate distances
        epi_distances=Distance_epicenter_calc(grid.source.coordinates{i,j},temp_x,temp_y);
        % Reshape vectors into grid
        epi_distances = reshape(epi_distances, length(grid.impact.xcoords), length(grid.impact.ycoords));
        % Insert distance grid into cell array
        grid.source.epi_distances{i,j} = epi_distances;
        
        % Update progress bar
        waitbar((((i - 1) * length(grid.source.ycoords)) + j) / (length(grid.source.xcoords) * length(grid.source.ycoords)),...
            f, sprintf('Distance Calculation Progress: %d %%',...
            floor((((i - 1) * length(grid.source.ycoords)) + j) / (length(grid.source.xcoords) * length(grid.source.ycoords))*100)));
    end
     
end

close(f) % closing distance calculation waitbar


%% %% Calculate nuisance for each source cell %%%%%

% Initialize cell arrays to hold results
grid.source.Nuisance_simulations = grid.source.coordinates;
grid.source.Nuisance_median = grid.source.coordinates;

% This process can take some time. Show progress with tracking bar.
f = waitbar(0, 'Starting Simulation');

for sourcei = 1:length(grid.source.xcoords)
    for sourcej = 1:length(grid.source.ycoords)
    
        % Set random state for source cell
        if strcmp(random_type, 'consistent')
            rng(SI);
        elseif strcmp(random_type, 'shuffle')
            rng('shuffle'); % shuffle state for source node
            SI = rng; % record new state
        else  % if random_type is neither 'consistent' or 'shuffle'
            error("Random_type must be defined as 'consistent' or 'shuffle'")
        end
        
        
        % Retrieve values for variables for distances below nuisance distance limit

        Zin = grid.source.epi_distances{sourcei,sourcej} <= nuisance_distance_limit; % create logical array
        
        % Retrieve data for impact cells within relevant distance
        distances2 = grid.source.epi_distances{sourcei,sourcej}(Zin);
        Pop3 = grid.impact.Pop(Zin);
        Vs30_3 = grid.impact.Vs30(Zin);
        Vs30_simulations_3 = grid.impact.Vs30_simulations(Zin);

        % And Retrieve values for variables for distances below damage distance limit
        Zin2 = distances2 <= damage_distance_limit; % create logical array
        distances4 = distances2(Zin2); % Extract distances for damage calculations
        Vs30_5 = Vs30_3(Zin2); % Extract the Vs30s for damage calculations
        % Vs30_std_5 is taken within the simulation loops to ensure a match
        % for a cell's Vs30 perturbation used in nuisance and damage.


        % Copy distances into matrix for all realizations
        distances2_realizations = repmat(distances2, [1,Nr]); % Nuisance distances
        distances4_realizations = repmat(distances4, [1,Nr]); % Damage distances

        % Assign event a random depth for each realization according to
        % depth variation weights
        z = randsample(z_d, Nr, true, wz);

        % Factor EQ depth into distance calculations
        distances3=sqrt(((double(distances2_realizations)).^(2))+((z).^(2)));
        distances5=sqrt(((double(distances4_realizations)).^(2))+((z).^(2)));
        
        % Vary parameter values for each realization
        dE = rand(1,Nr); % Random number for trailing seismicity generation for each realization
        dEp = normrnd(1,0.05,[1,Nr]); % Perturbation of total population for each realization
        dEb = normrnd(bm,db,[1, Nr]); % Perturbation b-value
        % Prepare Vs30 perturbations
        dEv_idx = randsample([1:length(Vs30_simulations_3{1})], Nr, true); % Select random Vs30 simulations
        dEv_interim = cell2mat(Vs30_simulations_3); % Convert cell array to matrix of Vs30 simulation values
        dEv_interim = dEv_interim(:,dEv_idx); % Select Vs30 values for selected simulations
        dEv = Vs30_3 - dEv_interim; % Calculate perturbation from mean Vs30 value - recognize this is a redundant approach
        dEv_damage = dEv(repmat(Zin2,[1,Nr])); % Subset of dEv for Damage

        % Apply the perturbation to Vs30s
        Vs30_4 = repmat(Vs30_3, [1,Nr]) + dEv; % for Nuisance (the back half of a redundant approach)
        Vs30_4(Vs30_4 < 100) = 100; % Ensure velocity is physically reasonable
        Vs30_6 = repmat(Vs30_5, [1, Nr]) + reshape(dEv_damage, length(Vs30_5), []); % for Damage
        Vs30_6(Vs30_6 < 100) = 100; % Ensure velocity is physically reasonable


        % Apply population perturbation
        Pop4 = Pop3 * dEp;

        % Extract the perturbed Populations for damage calculations
        clear Pop5
        for j = 1:Nr
            Pop5(:,j) = Pop4(Zin2,j);
        end


        for i=1:length(Mag_rng) % Loop through range of magnitudes

            % Reset random state for source cell
            rng(SI);
            
            %%% Create input array of Magnitudes, including any trailing seismicity
           
            Mag = ones([1,Nr]) * Mag_rng(i); % Magnitude at earthquake node

            % Perturb Magnitude according to trailing seismicity, if any
            if(trailing_seismicity)
                q = dE .* (1-exp(-(Rf.^dEb))) + exp(-(Rf.^dEb));
                dM = log10(Rf) ./ dEb-log10(-log(q))./dEb;
                Mag = Mag + dM;
            end

            %%% Calculate Ground Motions at impact cells given EQs at a
            % given source cells and magnitude for all realizations

            GM_cat = GM_calc(Mag,distances3,Vs30_4, 'PGV');

            % Convert from logarithmic scale and scale from cm to m
            GM2=power(10,GM_cat)/100;


            %%% Determine number of people nuisanced at impact cells for all realizations
            
            % Use GMs to calculate Nuisance at impact cells for all realizations
            Mu = Nuisance_function(GM2, 'PGV');

            % Scale population to number of households - assume 4 people
            % per household
            People_nuis = (Pop4 * 0.25) .* Mu;

            % Compile in array for simulations [Realization, Magnitude]
            People_nuis_tot(:,i)=sum(People_nuis, 1);
            
            
            %%% Determine number of households damaged at impact cells for
            % all realizations

            GM_dam = GM_calc(Mag,distances5,Vs30_6, 'PGA');
            
            % Convert from logarithmic scale and scale from cm to m
            GM3 = power(10,GM_dam)/100;

            % Calculate damage
            P = Damage_function(GM3);

            % Scale population to number of households - assume 4 people
            % per household
            People_dam=(Pop5*1/4) .* P;

            % Compile in matrix for simulations [Realization, Magnitudes]
            People_dam_tot(:,i)=sum(People_dam, 1);
            
        end % Close the loop of magnitudes

        % Calculate median nuisance and damage
        Nuisance_median = squeeze(median(People_nuis_tot, 1));
        Damage_median = squeeze(median(People_dam_tot, 1));

        % Insert Nuisance and Damage into source grid
        grid.source.Nuisance_simulations{sourcei, sourcej} = People_nuis_tot;
        grid.source.Nuisance_median{sourcei, sourcej} = Nuisance_median;

        grid.source.Damage_simulations{sourcei, sourcej} = People_dam_tot;
        grid.source.Damage_median{sourcei, sourcej} = Damage_median;

        
        % Update wait bar
        waitbar((((sourcei - 1) * length(grid.source.ycoords)) + sourcej) / (length(grid.source.xcoords) * length(grid.source.ycoords)),...
            f, sprintf('Simulation Progress: %d %%',...
            floor((((sourcei - 1) * length(grid.source.ycoords)) + sourcej) / (length(grid.source.xcoords) * length(grid.source.ycoords))*100)));

    end % Close the loop of sourcej

end % Close the loop of sourcei

close(f) % closing simulation waitbar


%% Save the data structure variable named grid

% Clear large distance and Vs30 grids
grid.source = rmfield(grid.source, 'epi_distances');
grid.impact = rmfield(grid.impact, 'Vs30_simulations');

% Save results
save(save_location, 'grid');