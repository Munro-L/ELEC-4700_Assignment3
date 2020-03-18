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
% electron
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

for i = 0:50      % each step is a femtosecond
    previous_particles = particles;
    
    % update positions
    particles(:, 1) = particles(:, 1) + particles(:, 4);
    particles(:, 2) = particles(:, 2) + particles(:, 3);
    
    % check if any particles passed the boundary and deal with them
    x_boundary_changes_right = particles(:, 2) > 200;
    if any(x_boundary_changes_right)
        particles(:, 2) = particles(:, 2) .* ~x_boundary_changes_right;
    end

    x_boundary_changes_left = particles(:, 2) < 0;
    if any(x_boundary_changes_left)
        particles(:, 2) = particles(:, 2) + 200 * x_boundary_changes_left - abs(particles(:, 2) .* x_boundary_changes_left);
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