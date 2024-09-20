function dz = oneDepileptor(z, x0, K,tau0)

    N = length(z); % number of nodes
    dz = zeros(N, 1); % initialize dz vector

    for i = 1:N
        F = (1/4) * (-16/3 - sqrt(8 * z(i) - 629.6 / 27)); % define F(z) function
        G = 4 * F - z(i); % define G(z) function
        x0_i = (-4) * x0(i); % Extract the x0 value of the i-th node

        % Calculate the coupling term
        coupling_term = 0;

        % Write the function form
        for j = 1:N
            coupling_term = coupling_term + K(i, j) * ((1/4) * (-16/3 - sqrt(8 * z(j) - 629.6 / 27)) - F);
        end

        dz(i) = (1/tau0)* (G + x0_i - coupling_term); 
    end
end
