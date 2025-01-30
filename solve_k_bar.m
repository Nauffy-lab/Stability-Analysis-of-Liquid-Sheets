function k_bar_solution = solve_k_bar(c_bar)
    U_1 = 1;   
    U_g = 0;    
    rho_bar = 998/1.2;  
    We_g = 0.82; 
    
    term1 = (U_1 - c_bar)^2 * rho_bar;  
    term2 = (U_g - c_bar)^2;            
    
    k_bar_solution = We_g * (term1 + term2);  
end




