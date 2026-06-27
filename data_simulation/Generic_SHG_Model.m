function I_SHG = Generic_SHG_Model( theta , theta_0 , p , q)
% Note : -1 < p < 0 ; -1 < q < 1
% theta and theta_0 is degree
% theta is Polarization angle of incident light
% I_SHG is SHG intensity of molecular dipole

I_SHG = (5 + p.^2 + q.^2) + ...
    4 * (1 + p) .* cosd(2 * (theta - theta_0)) + 2 * p .* cosd(4 * (theta - theta_0)) +...
    4 * q .* sind(2 * (theta - theta_0)) + 2 * q .* sind(4 * (theta - theta_0)) ;

end