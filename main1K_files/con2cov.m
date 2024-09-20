function [M,err] = con2cov(C,symmetric,maxi,tol,verbose)

% Calculates the covariance matrix for the stationary multivariate
% Ornstein-Uhlenbeck process:
%
%     dX(t) = -X(t)*(I-C)dt + dW_t
%
% where W_t is a multivariate Wiener process with identity covariance
% matrix and C is the connection matrix (1st index efferent)
%
% Implements [*] Eqs. 21-23
%
% Set 'symmetric' if C is symmetric to use non-iterative formula. Note that
% the non-symmetric version will also work for symmetric C, but *not* vice-
% versa!
%
% 'maxi' gives the maximum iterations for the non-symmetric version.
% Convergence criterion is to stop when the ratios of the added terms
% approach the machine epsilon with given tolerance 'tol'. If maximum
% iterations are exceeded a warning is issued. Failure to converge
% is treated as an error. You can play off accuracy against no. of iterations.
% Try maxi = 1000, tol = 10 for starters. If it converges you should get a max.
% relative error of order 1e-15 or better.
%
% Set 'verbose' to report number of iterations and convergence results
% for non-symmetric case

%************************************************************************%
%                                                                        %
% ncomp_tools: L. Barnett, C. L. Buckley and S. Bullock (2009)           %
%                                                                        %
% [*] See http://www.secse.net/ncomp/ncomp1.pdf                          %
%                                                                        %
% IF YOU USE THIS CODE PLEASE CITE THE ABOVE AS:                         %
%                                                                        %
%   L. Barnett, C. L. Buckley and S. Bullock (2009)                      %
%   On Neural Complexity and Structural Connectivity                     %
%   Physical Review E (in review)                                        %
%                                                                        %
%************************************************************************%

err = 0;

n = size(C,1);

I = eye(n);

if symmetric
	M = inv(I-C)/2;
	return
end

M = I/2;
dM = I/2;
for i = 1:maxi
	dM = (C'*dM+dM*C)/2;
	if negligible(dM,M,tol)
	   break
	end
	M = M + dM;
end

mre = maxrelerr(dM,M);

if isnan(mre) || isinf(mre)
	fprintf(2,'ERROR (con2cov): failed to converge\n');
	err = 2;
	return;
end

if i == maxi
	mac = mean(mean(abs(M)));
	fprintf(1,'WARNING (con2cov): exceeded maximum iterations: mean abs. covariance = %d, max. relative error = %d\n',mac,mre);
	err = 1;
	return;
end

if verbose
	fprintf(1,'con2cov: iterations = %i, max. relative error = %d\n',i,mre);
end

