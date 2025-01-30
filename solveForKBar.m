function k_bar = solveForKBar(omega_bar)
    % Constants
    U1_bar = 1;    % Defined within the function
    Ug_bar = 0;    % Defined within the function
    rho_bar = 998/1.2;   % Defined within the function
    We_bar = 0.82;    % Defined within the function
    
    syms k

    % Define the equation
    equation = (U1_bar^2 * rho_bar + Ug_bar^2) * k^2 ...
               - 2 * (omega_bar * U1_bar * rho_bar + omega_bar * Ug_bar) * k ...
               + omega_bar^2 * (rho_bar + 1) ...
               - k^3 / We_bar == 0;

    % Solve the cubic equation for k
    k_solutions = solve(equation, k);

    % Evaluate the solutions numerically
    k_vals = double(vpa(k_solutions));

    % Find the solution with the most negative imaginary part
    [~, idx] = min(imag(k_vals));
    k_bar = k_vals(idx);
end
