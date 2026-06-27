function I_SHG = Single_Molecular_Model( theta , theta_0 , Pitch_Angle)
b = 2 ./ tand(Pitch_Angle) ;

I_SHG = (sind(theta - theta_0).^2 + b.* cosd(theta - theta_0).^2).^2 + ... 
    4 * sind(theta - theta_0).^2 .* cosd(theta - theta_0).^2 ;
end