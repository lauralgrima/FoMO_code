function [com] = centerOfMass3D(X, Y, Z, M)
% Calculate the center of mass in 3 dimensions
% Inputs:
%   X, Y, Z: vectors of x, y, and z coordinates
%   M: vector of masses corresponding to each point
% Output:
%   com: [x, y, z] coordinates of the center of mass

% Check if the input vectors have the same length
if length(X) ~= length(Y) || length(X) ~= length(Z) || length(X) ~= length(M)
    error('Input vectors must have the same length.');
end

% Calculate the total mass
totalMass = sum(M);

% Calculate the center of mass coordinates
comX = sum(X .* M) / totalMass;
comY = sum(Y .* M) / totalMass;
comZ = sum(Z .* M) / totalMass;

% Return the center of mass as a vector
com = [comX, comY, comZ];

% Display the result
fprintf('Center of mass: [%.2f, %.2f, %.2f]\n', com(1), com(2), com(3));
