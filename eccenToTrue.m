function [T] = eccenToTrue(E, e)
%ECCENTOTRUE Summary of this function goes here
%   Detailed explanation goes here
T = atan2(sin(E) .* sqrt(1 - e.^2), cos(E) - e);
end

