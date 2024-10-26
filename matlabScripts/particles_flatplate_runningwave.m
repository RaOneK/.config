% MIT License
% 
% Copyright (c) 2024 Robert Schütze
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
% If you use this code for your own animations, please give credit to me.

% 2D Fluid Simulation in MATLAB (with viscous effects)
% Cleaned and commented version
% To increase performance decrease the number of Jacobi iterations, add less inflow
% particles or add them less often (e.g. only every 2nd simulation iteration)

% Simulation parameters
s = 200;                            % Grid size
ar = 1;                             % Aspect ratio (to represent a longer domain for flow past a flat plate)
J = [0 1 0; 1 0 1; 0 1 0]/4;        % Stencil for Jacobi method
% Reynolds number definition
Re = 100;  % Reynolds number
L = 1;  % Characteristic length (e.g., domain length or inlet width)
U = 1;  % Characteristic velocity (e.g., inflow velocity)
nu = U * L / Re;  % Viscosity coefficient derived from Reynolds number                          % Viscosity coefficient
A = 0;
omega = -5;
k = 10;
% Create a grid
[X, Y] = meshgrid(1:s*ar, 1:s);

% Initialize pressure and velocity fields
[p, vx, vy] = deal(zeros(s, s*ar));

% Initial positions of particles
[px, py] = meshgrid(10:15, 1:200);
px = reshape(px, numel(px), 1);
py = reshape(py, numel(py), 1);

% Save these initial positions for the inflow
pxo = px;
pyo = py;

% Create a flat plate along the bottom boundary (y = 1)
plate_y = 1;
vy(plate_y, :) = 0;
vx(plate_y, :) = 0;

f = figure(1); % Create figure to check if it's closed

% Main simulation loop (stops when closing figure window)
ctr = 1;
dt = 1e-1;
while ishandle(f)
    % Set initial velocity in a specific region (inlet flow)
    vx(:, 1:2) = (1+A*cos(omega*ctr*dt))*U;
    vy(:, 1:5) = 0;
    
    % Compute right-hand side for pressure equation
    rhs = -divergence(vx, vy);
    
    % Jacobi iteration to solve for pressure
    % Higher number of iterations yields better solution
    for i = 0:100
        p = conv2(p, J, 'same') + rhs/2;
    end
    
    % Compute velocity gradient and update velocities for non-boundary pixels
    [dx, dy] = gradient(p);
    vx(2:end-1, 2:end-1) = vx(2:end-1, 2:end-1) - dx(2:end-1, 2:end-1);
    vy(2:end-1, 2:end-1) = vy(2:end-1, 2:end-1) - dy(2:end-1, 2:end-1);   
    
    % Apply viscous effects using the diffusion equation
    [vx_diff_x, vx_diff_y] = gradient(vx);
    vx(2:end-1, 2:end-1) = vx(2:end-1, 2:end-1) + nu * (vx_diff_x(2:end-1, 2:end-1) + vx_diff_y(2:end-1, 2:end-1));
    [vy_diff_x, vy_diff_y] = gradient(vy);
    vy(2:end-1, 2:end-1) = vy(2:end-1, 2:end-1) + nu * (vy_diff_x(2:end-1, 2:end-1) + vy_diff_y(2:end-1, 2:end-1));
    
    % Apply open boundary condition at the right boundary
    vx(:, end) = vx(:, end-1);
    vy(:, end) = vy(:, end-1);
    p(:, end) = p(:, end-1);
    
    % Apply prescribed tangential velocity at the top boundary
    vx(end, :) = (1+A*cos(omega*ctr*dt + k*linspace(1,s*ar,s*ar)))*U;
    vy(end, :) = 0;
    
    % Enforce no-slip boundary condition at the flat plate
    vx(plate_y, :) = 0;
    vy(plate_y, :) = 0;
    
    % Advect velocity field using Runge-Kutta 4th order method (-1 = backward)
    [pvx, pvy] = RK4(X, Y, vx, vy, -1);
    vx = interp2(X, Y, vx, pvx, pvy, 'linear', 0);
    vy = interp2(X, Y, vy, pvx, pvy, 'linear', 0);  
    
    % Advect particles using Runge-Kutta 4th order method (1 = forward)
    [px, py] = RK4(px, py, vx, vy, 1);
    
    % Remove particles that leave the domain
    valid_idx = px > 0 & px <= s*ar & py > 0 & py <= s;
    px = px(valid_idx);
    py = py(valid_idx);
    
    % Add the inflow particles
    px = [px; pxo];
    py = [py; pyo];
    
    if mod(ctr,10)==0
        % Visualization of particle positions
        subplot(2,1,1);
        scatter(px, py, 1, 'filled');
        axis equal; 
        axis([0 s*ar 0 s]);
        title('Particle Positions');
        
        % Visualization of kinetic energy and streamlines
        subplot(2,1,2);
        kinetic_energy = 0.5 * (vx.^2 + vy.^2);
        imagesc(kinetic_energy);
        hold on;
        streamslice(X, Y, vx, vy);
        hold off;
        axis equal;
        title('Kinetic Energy and Streamlines');
        colorbar;
        drawnow;
    end

    ctr = ctr + 1;
end

% Function for Runge-Kutta 4th order method for advection
function [x_new, y_new] = RK4(px, py, vx, vy, h)
   k1x = interp2(vx, px, py, 'linear', 0);
   k1y = interp2(vy, px, py, 'linear', 0);
   k2x = interp2(vx, px + h/2 * k1x, py + h/2 * k1y, 'linear', 0);
   k2y = interp2(vy, px + h/2 * k1x, py + h/2 * k1y, 'linear', 0);
   k3x = interp2(vx, px + h/2 * k2x, py + h/2 * k2y, 'linear', 0);
   k3y = interp2(vy, px + h/2 * k2x, py + h/2 * k2y, 'linear', 0);
   k4x = interp2(vx, px + h * k3x, py + h * k3y, 'linear', 0);
   k4y = interp2(vy, px + h * k3x, py + h * k3y, 'linear', 0);
   x_new = px + h/6 * (k1x + 2*k2x + 2*k3x + k4x);
   y_new = py + h/6 * (k1y + 2*k2y + 2*k3y + k4y);
end
