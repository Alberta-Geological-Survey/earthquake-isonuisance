clear
%% %% Parameter definition %%%%%

% Input files
simulation_file = 'simulations.mat';
population_file = 'Population_AB_Ext_3km.mat';
provincial_boundary_file = 'Shape_Alberta.txt';

% Iso-Nuisance and Iso-Damage thresholds
% Acceptible number of households affected with 50% probability (via
% median curve)
nuisance_threshold = 30000;
damage_threshold = 3;

%%% Option: apply smoothing to iso-nuisance map
nuisance_smooth = true; % set as true or false
% Define a filter, H, for smoothing
% This replaces a cell's value with the mean of the cell and its four
% immediate neighbours.
H = (1/5) * [0,1,0; 1,1,1; 0,1,0];
% !! Note !! Applying a filter like this will leave an artifact around the
% edge of the map. Lines 81-83 compensate for the original
% (1/5) * [0,1,0; 1,1,1; 0,1,0] filter. Other filters would require
% different compensation.

%% %% Setup %%%%%
% Import grids for visualization
Pop=importdata(population_file); % Input Population data

% Import provincial boundary
AB_Boundary=importdata(provincial_boundary_file);
Lat_AB =AB_Boundary(:,2); % X-Coordinates
Lon_AB =AB_Boundary(:,1) ; % Y-Coordinates

% Load simulations as variable grid
% Remember, array for simulations: [Realization, Magnitudes]
load(simulation_file, 'grid'); % Loads structure array as variable 'grid'

% Retrieve number of simulation realizations
Nr = size(grid.source.Nuisance_simulations{1,1}, 1);

% Extract range of Magnitudes simulated from structure array
Mag_rng = grid.source.Magnitudes;


%% Create an Iso-Nuisance map

% Initialize an iso-nuisance grid
grid.source.Nuisance_threshold_magnitude = zeros(length(grid.source.xcoords), length(grid.source.ycoords)); % default value is 0

% Cycle through iso-nuisance grid to assign magnitudes that reach nuisance
% threshold. Values are interpolated from median nuisance curves.
for i = 1:length(grid.source.xcoords)
    for j = 1:length(grid.source.ycoords)
        if sum(grid.source.Nuisance_median{i,j}) == 0 % Edge case: if no households are nuisanced for largest magnitude...
            grid.source.Nuisance_threshold_magnitude(i,j) = max(Mag_rng); % ... assign largest magnitude as threshold magnitude
        else
            % In general, interpolate threshold magnitude from median curves
            grid.source.Nuisance_threshold_magnitude(i,j) = ...
                interp1(grid.source.Nuisance_median{i,j}, Mag_rng, nuisance_threshold, 'linear', 'extrap');
        end
    end
end

% Apply a minimum and maximum magnitude
Mag_max = max(Mag_rng);
Mag_min = min(Mag_rng);
grid.source.Nuisance_threshold_magnitude(grid.source.Nuisance_threshold_magnitude < Mag_min) = Mag_min;
grid.source.Nuisance_threshold_magnitude(grid.source.Nuisance_threshold_magnitude > Mag_max) = Mag_max;


%%% Apply smoothing, if desired
if nuisance_smooth == true
    % Retain an un-filtered version of the iso-nuisance map
    grid.source.Nuisance_threshold_magnitude_unfiltered = grid.source.Nuisance_threshold_magnitude;
    
    % Apply the filter to the iso-nuisance map
    grid.source.Nuisance_threshold_magnitude = filter2(H, grid.source.Nuisance_threshold_magnitude);
    
    % Compensate for the halo this filter creates around the edge of the map
    grid.source.Nuisance_threshold_magnitude(1:end,[1,end]) = 5/4 * grid.source.Nuisance_threshold_magnitude(1:end,[1,end]); % Top edges
    grid.source.Nuisance_threshold_magnitude([1,end],2:end-1) = 5/4 * grid.source.Nuisance_threshold_magnitude([1,end],2:end-1); % Side edges
    grid.source.Nuisance_threshold_magnitude([1,end],[1,end]) = 4/3 * grid.source.Nuisance_threshold_magnitude([1,end],[1,end]); % Additional compensation for corners

end

%% Create an Iso-Damage Map

% Initialize an iso-damage grid
grid.source.Damage_threshold_magnitude = zeros(length(grid.source.xcoords), length(grid.source.ycoords)); % default value is 0

% Cycle through iso-damage grid to assign magnitudes that reach damage
% threshold. Values are interpolated from median damage curves.
for i = 1:length(grid.source.xcoords)
    for j = 1:length(grid.source.ycoords)
        if sum(grid.source.Damage_median{i,j}) == 0 % Edge case: if no households are nuisanced for largest magnitude...
            grid.source.Damage_threshold_magnitude(i,j) = max(Mag_rng); % ... assign largest magnitude as threshold magnitude
        else
            % In general, interpolate threshold magnitude from median curves
            grid.source.Damage_threshold_magnitude(i,j) = ...
                interp1(grid.source.Damage_median{i,j}, Mag_rng, damage_threshold, 'linear', 'extrap');
        end
    end
end

% Apply a minimum and maximum Magnitude
grid.source.Damage_threshold_magnitude(grid.source.Damage_threshold_magnitude < Mag_min) = Mag_min;
grid.source.Damage_threshold_magnitude(grid.source.Damage_threshold_magnitude > Mag_max) = Mag_max;


%% Create a combined map

% Define the combined map as the minimum of the iso-damage and iso-nuisance
% maps for each source grid cell.
grid.source.Combined_threshold_magnitude = min(grid.source.Nuisance_threshold_magnitude,...
    grid.source.Damage_threshold_magnitude);

%% Plotting maps

% Initialize figure
f = figure('Name', 'Iso-Nuisance and Iso-Damage Mapping');
sgtitle("Iso-Nuisance and Iso-Damage Mapping")
f.Position = [10, 50, 1600, 400];

% Plot population density
subplot(1,4,1);
[C,h]=contourf(grid.impact.xcoords', grid.impact.ycoords', log(grid.impact.Pop+1)', 50);
set(h,'LineColor','none')
pbaspect([1,2,1])
cc = colorbar;
cc.Label.String = 'log(Population Density)';
xlabel('Longitude')
ylabel('Latitude')
title ('Population Denisty Map')
set(gcf,'color','w');
hold on

plot(Lon_AB,Lat_AB,'r','LineWidth',1.5); % Add geographic context


% Plot Iso-Nuisance
subplot(1,4,2)
[C,h]=contourf(grid.source.xcoords', grid.source.ycoords', grid.source.Nuisance_threshold_magnitude', 50);
set(h,'LineColor','none')
pbaspect([1,2,1])
cc = colorbar;
cc.Label.String = 'Magnitude';
xlabel('Longitude')
ylabel('Latitude')
title ('Iso-Nuisance Map')
set(gcf,'color','w');
xlim([min(grid.impact.xcoords), max(grid.impact.xcoords)])
ylim([min(grid.impact.ycoords), max(grid.impact.ycoords)])
caxis([min(Mag_rng) max(Mag_rng)])
hold on

plot(Lon_AB,Lat_AB,'r','LineWidth',1.5); % Add geographic context


% Plot iso-Damage
subplot(1,4,3)
[C,h]=contourf(grid.source.xcoords', grid.source.ycoords', grid.source.Damage_threshold_magnitude', 50);
set(h,'LineColor','none')
pbaspect([1,2,1])
cc = colorbar;
cc.Label.String = 'Magnitude';
xlabel('Longitude')
ylabel('Latitude')
title ('Iso-Damage Map')
set(gcf,'color','w');
xlim([min(grid.impact.xcoords), max(grid.impact.xcoords)])
ylim([min(grid.impact.ycoords), max(grid.impact.ycoords)])
caxis([min(Mag_rng) max(Mag_rng)])

hold on

plot(Lon_AB,Lat_AB,'r','LineWidth',1.5); % Add geographic context


% Plot a combined map
subplot(1,4,4)
[C,h]=contourf(grid.source.xcoords', grid.source.ycoords', grid.source.Combined_threshold_magnitude', 50);
set(h,'LineColor','none')
pbaspect([1,2,1])
cc = colorbar;
cc.Label.String = 'Magnitude';
xlabel('Longitude')
ylabel('Latitude')
title ('Combination Iso-Risk Map')
set(gcf,'color','w');
xlim([min(grid.impact.xcoords), max(grid.impact.xcoords)])
ylim([min(grid.impact.ycoords), max(grid.impact.ycoords)])
caxis([min(Mag_rng) max(Mag_rng)])
hold on

plot(Lon_AB,Lat_AB,'r','LineWidth',1.5); % Add geographic context




%% Plotting curves

% Define coordinates of example point of interest
xy=[-117.3592 54.4184]; % X, Y coordinates for location near Fox Creek, AB

% Find source grid indicies for example point of interest
[val, sourcei] = min(abs(grid.source.xcoords-xy(1)));
[val, sourcej] = min(abs(grid.source.ycoords-xy(2)));

% Initialize figure
f = figure('Name', 'Nuisance and Damage Curves');
sgtitle('Example Nuisance-Magnitude and Damage-Magnitude Curves')
f.Position = [10, 50, 800, 400];

subplot(1,2,1)
for k=1:Nr
    semilogy(Mag_rng,grid.source.Nuisance_simulations{sourcei,sourcej}(k,:),'r');
    hold on 
end
semilogy(Mag_rng,grid.source.Nuisance_median{sourcei,sourcej},'b','LineWidth',2);
title(["Households Nuisanced vs Magnitude", "for Example Location"])
xlabel('Magnitude')
ylabel('Number of Households')
set(gcf,'color','w');


subplot(1,2,2)
for k=1:Nr
    semilogy(Mag_rng,grid.source.Damage_simulations{sourcei,sourcej}(k,:),'r');
    hold on 
end
semilogy(Mag_rng,grid.source.Damage_median{sourcei,sourcej},'b','LineWidth',2);
title(["Households Damaged vs Magnitude", "for Example Location"])
xlabel('Magnitude')
ylabel('Number of Households')
set(gcf,'color','w');