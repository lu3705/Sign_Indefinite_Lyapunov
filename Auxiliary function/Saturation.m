% Saturation: applies input saturation with dead-zone feedback
% Inputs:
%   x  - system state
%   K1 - state feedback gain
%   K2 - dead-zone gain
% Outputs:
%   u     - control input after dead-zone correction
%   u_sat - saturated control input
%   dz    - dead-zone signal

function [u, u_sat, dz] = Saturation(x, K1, K2)

    u = K1*x;
    u_sat = min(1, max(-1, u));
    dz = u - u_sat;

    u = K1*x + K2*dz;
    u_sat = min(1, max(-1, u));
    dz = u - u_sat;

end

