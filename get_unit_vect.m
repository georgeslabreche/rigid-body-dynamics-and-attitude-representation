function [ unit_vect ] = get_unit_vect( lambda_subbed )
%GET_UNIT_VECT Summary of this function goes here
%   Detailed explanation goes here
% Subsitute these back in JJ-lambda*eye(3)

syms y z
solved_x = 1;
X = [solved_x; y; z];

fun = lambda_subbed * X == [0; 0; 0];
fun_sol = solve([fun(2), fun(3)], [y, z]);

syms a b c
unit_vect = (solved_x*a + fun_sol.y*b + fun_sol.z*c) / sqrt(solved_x^2 +  fun_sol.y^2 + fun_sol.z^2);

end
