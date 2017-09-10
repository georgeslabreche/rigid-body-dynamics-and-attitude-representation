function [ unit_vect ] = get_unit_vect(JJ, lambda)
%GET_UNIT_VECT Given an inertia matrix calculate and return the unit vector
% for a given principal moment of inertia.

syms lambda_var
lambda_subbed = subs(JJ-lambda_var*eye(3), lambda_var, lambda);

% Since the determinant of the coefficient matrix is zero, at most two of
% the three equations are independent. Thus, at most two of the three
% components of the vector can be found in terms of the third. We can
% therefore arbitrarily set x=1 and sovlve for y and z using any two of the
% independent equations.
%
% Curtis, Howard. "Moments of Interia" in Orbital Mechanics: For
% Engineering Students, 515.
syms y z
solved_x = 1;
X = [solved_x; y; z];

% Define and solve the x and y equations:
fun = lambda_subbed * X == [0; 0; 0];
fun_sol = solve([fun(2), fun(3)], [y, z]);

% Calculate and return the unit vector.
syms a b c
unit_vect = (solved_x*a + fun_sol.y*b + fun_sol.z*c) /...
    sqrt((solved_x)^2 +  (fun_sol.y)^2 + (fun_sol.z)^2);

end
