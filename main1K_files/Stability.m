function exp_sigma_sqr = Stability(covariance)

N = size(covariance,1);

exp_sigma_sqr = trace(covariance) ./ N;

end

