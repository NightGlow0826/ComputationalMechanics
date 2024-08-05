clear;
clc;

% This is a program to solve PDE with Alternative Direction Implicit(ADI) method
% See CFD4_0418 for the problem
%

% U_next
sigma = 1;
dx = 0.1;
dy = dx;
N = 1 / dx;
dt = 0.1;
steps = 2/dt;
k = (sigma * dt) / (2 * dx.^2);

X = linspace(0, 1, N+1);
Y = linspace(0, 1, N+1);
1
analytical_solution = analytical(X, Y, dt*steps);
U_numerical = get_numerical(dx, 1);

get_numerical_err();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function a = analytical(X, Y, t)
sigma = 1;
l = length(X);
a = zeros(l, l);
for i = 1:l
    for j = 1:l
        x = X(j);
        y = Y(i);
        a(i, j) = 20 + 80 * (y - exp(-0.5* sigma * pi^2 * t) * sin(pi * x/ 2) * sin(pi* y/ 2));
    end
end
end

function init_u = get_init(N) % It doesnt fit the analytical
init_u = zeros(N+1, N+1);
dx = 1 / N;
dy = 1 / N;
for iter = 1:N+1
    x = (iter - 1) * dx;
    y = (iter - 1) * dy;
    init_u(1, iter) = y0_b(x, 0);
    init_u(iter, 1) = x0_b(y, 0);
    init_u(end, iter) = y1_b(x, 0);
    init_u(iter, end) = x1_b(y, 0);
end
end


function implicit_mat = get_implicit_mat(N, k)
implicit_mat = zeros(N+1, N+1);
implicit_mat(1, 1) = 1;
implicit_mat(N+1, N+1) = 1;
for row = 2 : N
    implicit_mat(row, [row-1, row, row+1]) = [k, -2*k-1, k];
end
end

function x0_b = x0_b(y, t)
x0_b = 20 + 80*y;
end

function x1_b = x1_b(y, t)
x1_b = 20 + 80*(y - exp(-0.5*1*pi.^2*t)* sin(pi/2*y));
end

function y0_b = y0_b(x, t)
y0_b = 20;
end

function y1_b = y1_b(x, t)
y1_b = 20 + 80*(1 - exp(-0.5*1*pi.^2*t)* sin(pi/2*x));
end

function bc = ustar_bc(f, var, t, dvar, dt, k)
bc = (f(var, t) + f(var, t+dt))/2 - ...,
    1/2 * k * (f(var - dvar, t+dt) - 2 * f(var, t+dt) + f(var + dvar, t+dt)) + ...,
    1/2 * k * (f(var - dvar, t) - 2 * f(var, t) + f(var + dvar, t));
end

function U_numerical = get_numerical(dx, compute_err)
sigma = 1;
dx = dx;
dy = dx;
N = fix(1 / dx);
dt = 0.1;
k = (sigma * dt) / (2 * dx.^2);



X = linspace(0, 1, N+1);
Y = linspace(0, 1, N+1);
% U_N = get_init(N) % Should be initialize to zero to fit IC
U_N = analytical(X, Y, 0);
U_star = zeros(N+1, N+1);
U_next = zeros(N+1, N+1);

err = [];
steps = 2 / dt;
for timestep = 0:steps
    t = timestep * dt;
    % first step: calculate  U star. Firstly, traverse i direction in U star
    % Notice: i, j is reversed from the textbook
    composed = zeros(N+1, N+1);
    for i = 2:N
        composed(i, :) = -k*(U_N(i-1, :) + U_N(i+1, :)) + (2*k- 1) * U_N(i, :);
    end
    y = ((1:N+1) - 1)*dy;
    % Although at x bound, RHS could be composed of U_N, but we have already
    % Used the boundary condition
    
    U_star(:, 1) = ustar_bc(@x0_b, y, t, dy, dt, k);
    composed(:, 1) = ustar_bc(@x0_b, y, t, dy, dt, k);
    U_star(:, N+1) = ustar_bc(@x1_b, y, t, dy, dt, k);
    composed(:, N+1) = ustar_bc(@x1_b, y, t, dy, dt, k);
    % The y = 0, 1 of ustar is not important, since when calculation U_next, y_bound are given.
    implicit_mat = get_implicit_mat(N, k)';
    % Since ii's row combination, multiple it at right side with transpose
    
    % implicit_mat
    U_star(2:N, :) = composed(2:N, :) *inv(implicit_mat);
    
    % U_star to U_next
    composed = zeros(N+1, N+1);
    
    y = ((1:N+1) - 1) * dy;
    U_next(:, 1) = x0_b(y, t+dt);
    U_next(:, end) = x1_b(y, t+dt);
    for j = 2:N
        composed(:, j) = -k*(U_star(:, j - 1) + U_star(:, j+1)) + (2*k - 1) * U_star(:, j);
    end
    
    % need more constrains for y = 0, 1
    for j = 2:N
        x = (j - 1) * dx;
        composed(1, j) = y0_b(x, t+dt);
        composed(end, j) = y1_b(x, t+dt);
    end
    
    implicit_mat = get_implicit_mat(N, k);
    U_next(:, 2:N) = implicit_mat \ composed(:, 2:N);
    err(timestep+1) = max(max(abs(U_next - U_N)));
    U_N = U_next;
    
end % This end for timestep
if compute_err == 1
    plot(0:steps, log(err), '.-')
end

U_numerical = U_next;
end

function err = get_numerical_err()
expo = linspace(0, 2, 100);
meshes = fix(10 * 3 .^ expo);
dx = 1 ./ meshes;

err = [];
for i = 1:length(expo)
    X = linspace(0, 1, meshes(i) + 1);
    err_mat = abs(analytical(X, X, 2) - get_numerical(dx(i), 0));
    err(i) = max(max(abs(err_mat))); %#ok<AGROW>
end
plot(log(dx), log(err), '*-')
ylabel('log err')
xlabel('log h')
end