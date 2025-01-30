function [f_value, df_deta_value, d2f_deta2_value] = solve_blasius_at_eta(eta_input)
    % Initial guess for the solution
    solinit = bvpinit(linspace(0, 10, 100), @init_guess);
    
    % Solve the differential equation using bvp4c
    sol = bvp4c(@blasius_eq, @blasius_bc, solinit);
    
    % Extract the solution at the given eta_input
    f_value = deval(sol, eta_input, 1);          % f(eta_input)
    df_deta_value = deval(sol, eta_input, 2);    % df/deta(eta_input)
    d2f_deta2_value = deval(sol, eta_input, 3);  % d^2f/deta^2(eta_input)
end

% Define the Blasius equation as a system of first-order ODEs
function dydeta = blasius_eq(eta, y)
    dydeta = [y(2);          % dy1/deta = y2
              y(3);          % dy2/deta = y3
              -0.5 * y(1) * y(3)];  % dy3/deta = -(1/2) * y1 * y3
end

% Define the boundary conditions
function res = blasius_bc(ya, yb)
    res = [ya(1);       % f(0) = 0
           ya(2)-1;       % df/deta(0) = 0
           yb(2)];  % df/deta(eta -> infinity) = 1
end

% Initial guess for the solution
function yinit = init_guess(eta)
    yinit = [eta^2/2; eta; 1];  % Initial guess based on a simple profile
end

