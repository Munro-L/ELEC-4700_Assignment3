%% Question 1
% setup
clearvars
clearvars -GLOBAL
close all

% compute electric field, assuming its uniform in semiconductor
V = 0.1;
E = -V/200;     % in volts per nanometer

% compute force on each electron, as well as acceleration to make life easy
% later
F = E * 1.60217662E-26; % stay in nm, not m
m = 9.10938356E-31;
a = F/m;

% The following section of code is mostly pulled directly from the first
% part of Assignment 1, however with the addition of the force on each
% electron, and computation for current over time

% current is coulombs/s over a cross sectional area, so given the density
% is 10E15/cm^2, current can be measured by how many particles pass over a
% given slice of the semiconductor. In this case, we might as well put
% this in the boundary case for the left/right sides of the box. Every
% electron that crosses the boundary contributes an elementary charge in
% one time step.

num_particles = 10;
colours = ["r", "g", "b", "c", "y", "k", "m"];

kb = 1.38064852;
T = 300;

vth = sqrt(kb * T / m) / 1E15;  % scaled to femtoseconds
mean_time_collision = 0.2;   % measured in picoseconds
timesteps = 1E-3;            % also in picoseconds, so this will be 1 femto second timesteps

% Creates the positions of the particles and assigned a velocity to each,
% chosen from a normal distrobution shifted by the thermal velocity

% particle positions
particles = rand(num_particles, 2);
particles(:, 2) = particles(:, 2)*200;  % x-coordinates
particles(:, 1) = particles(:, 1)*100;  % y-coordinates

% normal distribution shifted by themal velocity
% column 3 is x-velocity, 4 is x-velocity
angles = randn(num_particles, 1) .* 2 * pi;
particles(:, 3) = randn(num_particles, 1) + vth*cos(angles);
particles(:, 4) = randn(num_particles, 1) + vth*sin(angles);

current_tracker = [];
for i = 0:50      % each step is a femtosecond
    previous_particles = particles;
    current_count = 0;
    
    % update positions
    particles(:, 1) = particles(:, 1) + particles(:, 4);
    particles(:, 2) = particles(:, 2) + particles(:, 3);
    
    % check if any particles passed the boundary and deal with them
    x_boundary_changes_right = particles(:, 2) > 200;
    if any(x_boundary_changes_right)
        particles(:, 2) = particles(:, 2) .* ~x_boundary_changes_right;
        current_count = current_count + nnz(x_boundary_changes_right);
    end

    x_boundary_changes_left = particles(:, 2) < 0;
    if any(x_boundary_changes_left)
        particles(:, 2) = particles(:, 2) + 200 * x_boundary_changes_left - abs(particles(:, 2) .* x_boundary_changes_left);
        current_count = current_count + nnz(x_boundary_changes_left);
    end
    
    y_boundary_changes_upper = particles(:, 1) > 100;
    if any(y_boundary_changes_upper)
        particles(:, 4) = particles(:, 4) - (2 * particles(:, 4) .* y_boundary_changes_upper);
        overshoot = (particles(:, 1) - 100) .* y_boundary_changes_upper;
        particles(:, 1) = particles(:, 1) - 2 * overshoot;
    end
    
    y_boundary_changes_lower = particles(:, 1) < 0;
    if any(y_boundary_changes_lower)
        particles(:, 4) = particles(:, 4) - (2 * particles(:, 4) .* y_boundary_changes_lower);
        overshoot = abs(particles(:, 1)) .* y_boundary_changes_lower;
        particles(:, 1) = particles(:, 1) + 2 * overshoot;
    end
    
    % update velocities based on acceleration from elecric field
    % theta = 0 because applied field is only horizontal, code will work
    % for x/y component with a different angle if set here
    theta = 0;
    particles(:,3) = particles(:,3) + a*cos(theta);
    particles(:,4) = particles(:,4) + a*sin(theta);
    
    % add current to array tracking it per timestep
    current_tracker = [current_tracker, current_count*1.602E-19];
    
    % plot position updates of particles, 
    % ignoring those that passed the x-boundary
    x_boundary_affected = x_boundary_changes_right | x_boundary_changes_left;
    temp_avg = mean(((sqrt(particles(:, 3).^2 + particles(:, 4).^2) .* 1E15).^2) .* m ./ kb);
    title(sprintf("Avg Temperature: %s", temp_avg))
    for i = 1:length(particles)
        if ~x_boundary_affected(i)
            plot([previous_particles(i, 2) particles(i, 2)], [previous_particles(i, 1) particles(i, 1)], colours(mod(i, length(colours)) + 1))
        end
    end
    axis([0 200 0 100])
    hold on
    pause(0.1)
end

% plot current vs. time
figure
plot(1:51, current_tracker)
title('Plot of current vs. timestep')

% This section constructs a density map based on the final positions of all
% particles.
density_data = ceil((particles(:, 1:2)./10)).*10;      % fancy way of rounding to the nearest 10
density_data(:,2) = mod(density_data(:,2), 200); 
density_data(:,2) = density_data(:,2) + ~density_data(:,2);  % janky way of making sure there are no zeros, set to 1 instead
grid = zeros(100, 200);
for i = 1:length(density_data)
    grid(density_data(i, 1), density_data(i, 2)) = grid(density_data(i, 1), density_data(i, 2)) + 1;
end
figure
imagesc(grid)
title('Density Plot')

% Temperature map can be done in a similar way as density map, except
% instead of considering particle positions, we look at the velocities and
% find the temperature of each particle to show in the density map with their position.
heat_data = (particles(:, 3:4) .* 1E15 .* m ./ kb).^2;
position_data = ceil((particles(:, 1:2)./10)).*10;
position_data(:,2) = mod(position_data(:,2), 200); 
position_data(:,2) = position_data(:,2) + ~position_data(:,2);  % janky way of making sure there are no zeros, set to 1 instead
grid = zeros(100, 200);
for i = 1:length(position_data)
    grid(position_data(i, 1), position_data(i, 2)) = grid(position_data(i, 1), position_data(i, 2)) + sqrt(heat_data(i, 1)^2 + heat_data(i, 2)^2);
end
figure
imagesc(grid)
title('Temperature Plot')

