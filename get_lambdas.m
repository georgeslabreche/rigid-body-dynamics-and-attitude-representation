function [ sol_lambda ] = get_lambdas(JJ)
%GET_LAMBDAS Summary of this function goes here
%   Detailed explanation goes here

syms lambda
det_fn = det(JJ-lambda*eye(3)) == 0;
sol_lambda = solve(det_fn, lambda);

end
