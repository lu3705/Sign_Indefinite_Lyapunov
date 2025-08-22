function [u, u_sat, dz] = Saturation(x, K1, K2)

    u = K1*x;
    u_sat = min(1, max(-1, u));
    dz = u - u_sat;

    u = K1*x + K2*dz;
    u_sat = min(1, max(-1, u));
    dz = u - u_sat;

end
