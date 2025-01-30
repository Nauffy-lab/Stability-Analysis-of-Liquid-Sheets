function v_hat_at_zero = solve_v_hat_at_zero_user(k_bar, c_bar, Re_g, x_bar, H)
    options = bvpset('RelTol', 1e-3, 'AbsTol', 1e-5);
    solinit_xi = bvpinit(linspace(0, 10, 100), [0.5 0 0]);
    solinit_vhat = bvpinit(linspace(0, 10, 100), [1 0 0]);
    sol_xi = bvp4c(@ode_system_xi, @bc_conditions_xi, solinit_xi, options);
    y_tilde = linspace(0, 10, 100);
    xi_solution = deval(sol_xi, y_tilde);
    u_BL = u_BL_func(xi_solution(2,:));
    v_BL = v_BL_func(xi_solution(1,:), xi_solution(2,:), y_tilde, Re_g, x_bar, H);
    sol_vhat = bvp4c(@(y_tilde, v) ode_system_vhat_scalar(y_tilde, v, u_BL, v_BL, k_bar, c_bar, 1j, @delta_func, @ddelta_dx_func, x_bar, H, Re_g), ...
        @(ya, yb) bc_conditions_vhat(ya, yb, c_bar, k_bar, H, x_bar), solinit_vhat, options);
    v_hat_at_zero = deval(sol_vhat, 0);
    v_hat_at_zero = v_hat_at_zero(1);
end

function dydx = ode_system_xi(y_tilde, xi)
    dydx = [xi(2); 
            xi(3); 
            -xi(1) * xi(3) / 2];
end

function res = bc_conditions_xi(ya, yb)
    res = [ya(1);  
           ya(2) - 1;  
           yb(1)];
end

function dydx = ode_system_vhat_scalar(y_tilde, v, u_BL, v_BL, k_bar, c_bar, i, delta_func, ddelta_dx_func, x_bar, H, Re_g)
    u_BL_val = interp1(linspace(0, 10, length(u_BL)), u_BL, y_tilde, 'linear', 'extrap');
    v_BL_val = interp1(linspace(0, 10, length(v_BL)), v_BL, y_tilde, 'linear', 'extrap');
    du_BL_dy = gradient(u_BL, linspace(0, 10, length(u_BL)));
    du_BL_dy_val = interp1(linspace(0, 10, length(du_BL_dy)), du_BL_dy, y_tilde, 'linear', 'extrap');
    k_coef = k_bar / i;
    epsilon = 1e-8;
    term1 = (-1 * k_coef * (u_BL_val - c_bar) + (-1/2) * y_tilde * (1/H) * (1/x_bar) * du_BL_dy_val) * v(2);
    term2 = sqrt(Re_g * x_bar) * v_BL_val * v(3);
    term3 = k_coef * du_BL_dy_val * v(1);
    term4 = k_bar^2 * (delta_func(x_bar)/H) * exp(-k_bar * y_tilde * (delta_func(x_bar)/H));
    dydx = zeros(3, 1);
    dydx(1) = v(2);
    dydx(2) = v(3);
    dydx(3) = (term1 + term3 + term4 + epsilon) + term2;
end

function res = bc_conditions_vhat(ya, yb, c_bar, k_bar, H, x_bar)
    res = [ya(1) + bval_1_func(c_bar, k_bar, H, x_bar) * ya(2);
           yb(1) - 1i / bval_2_func(c_bar);
           yb(2) + bval_3_func(c_bar, k_bar, H, x_bar) * 1i];
end

function u_BL = u_BL_func(dxi_dy_tilde)
    u_BL = dxi_dy_tilde;
end

function v_BL = v_BL_func(xi, dxi_dy_tilde, y_tilde, Re_g, x_bar, H)
    delta = @(x) H*sqrt(x / Re_g);
    ddelta_dx = @(x) 0.5*H / sqrt(Re_g * x);
    v_BL = ddelta_dx(x_bar) * (dxi_dy_tilde .* y_tilde - xi);
end

function delta_val = delta_func(x_bar)
    Re_g = 378;
    H = 0.672*10^-3;
    delta_val = H * sqrt(x_bar / Re_g);
end

function ddelta_dx_val = ddelta_dx_func(x_bar)
    Re_g = 378;
    H = 0.672*10^-3;
    ddelta_dx_val = 0.5 * H / sqrt(Re_g * x_bar);
end

function bval_1_val = bval_1_func(c_bar, k_bar, H, x_bar)
    bval_1_val = ((1 - c_bar) / ((k_bar * (delta_func(x_bar) / H) * (1 - c_bar)) - 0.4439));
end

function bval_2_val = bval_2_func(c_bar)
    bval_2_val = c_bar;
end

function bval_3_val = bval_3_func(c_bar, k_bar, H, x_bar)
    bval_3_val = k_bar * (delta_func(x_bar) / H) * (1 / c_bar);
end
