 function C = CouplingMatrix(z_fixed_numerical, K,tau0)
    
    N = length(z_fixed_numerical); % number of nodes
    C = zeros(N, N); % initialize matrix C						
    
    % calculate F'(z_fixed_numerical)
    F_prime = @(z, i) -1/sqrt(8 * z(i) - 629.6 / 27); % define F(z) function
    
    for i = 1:N
        
        sum_K_ij = sum(K(i, :)); % calculate Σ K_ij
        sum_K_ii = sum(K(i, i)); % % calculate Σ K_ii
        
        % diagonal term 
        C(i, i) = 1-(1/tau0)+ (1/tau0) *(4 + sum_K_ij - K(i, i)) * F_prime(z_fixed_numerical,i);
        
        % non-diagonal term
        for j = 1:N
            if i ~= j
                C(i, j) = -(1/tau0) *(K(j,i)) * F_prime(z_fixed_numerical,j);
            end
        end
    end
end


