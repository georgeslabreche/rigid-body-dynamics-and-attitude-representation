function [ sol_lambda ] = get_principal_moments_of_inertia(JJ)
%GET_PRINCIPAL_MOMENTS_OF_INTERIA Given an internia matrix, calculate and
%return the principal momements of intertia (known  as lambdas).

syms lambda
det_fn = det(JJ-lambda*eye(size(JJ,1))) == 0;
sol_lambda = solve(det_fn, lambda);

end
