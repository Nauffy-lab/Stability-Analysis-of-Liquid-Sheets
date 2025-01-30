function c_bar_solution = solve_c_bar_symbolic(v_hat_0, x_bar, k_bar)
    % Constants
    rho = 998 / 1.2; % Density ratio
    U_1 = 1; % Reference velocity (unitless)
    We_g = 0.82; % Weber number for the gas phase
    Re_g = 378; % Reynolds number for the gas phase
    delta = sqrt(x_bar / Re_g); % Boundary layer thickness
    du_BL_dy_tilde_0 = -0.4439; % Derivative of boundary layer velocity profile at y=0

    % Symbolic variable
    syms c_bar;

    % Construct the equation
    U_c_bar_diff = (U_1 - c_bar);
    term1 = v_hat_0 * rho*tanh(k_bar)* U_c_bar_diff^2;
    term2 = U_c_bar_diff * 1i;
    term3 = v_hat_0 * (k_bar / We_g + (2 / (Re_g * delta)) * du_BL_dy_tilde_0 * (1i - 1 / (2 * x_bar * k_bar)));

    equation = term1 - term2 - term3 == 0;

    % Solve the equation symbolically
    c_bar_solutions_sym = solve(equation, c_bar);

    % Evaluate solutions numerically and keep the precision to four decimal places
    c_bar_solutions_numeric = vpa(c_bar_solutions_sym, 8);  % Using higher precision internally

    % Find the solution with the most negative imaginary part
    imag_parts = imag(c_bar_solutions_numeric);
    [~, idx] = min(imag_parts);  % Get the index of the solution with the most negative imaginary part
    selected_solution = c_bar_solutions_numeric(idx);

    % Convert to double and round to four decimal places
    c_bar_solution = round(real(selected_solution), 4) + 1i * round(imag(selected_solution), 4);

end
