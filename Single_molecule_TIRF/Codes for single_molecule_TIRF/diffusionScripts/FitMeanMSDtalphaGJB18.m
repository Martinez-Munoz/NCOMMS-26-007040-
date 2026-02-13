function F = FitMeanMSDtalphaGJB18(x, SigmaMSD, Tmsd, AvarMSD, fitlength);  
% weighted fit: minimize: weightfactor*( 4D*time^alpha + offsetAtX=0  -  MSD)
% x(1) = 4D;
% x(2) = offset;
% x(3) = alpha.

F=(1./SigmaMSD(1:fitlength)).*( (Tmsd(1:fitlength).^x(3)).*x(1) + x(2) - AvarMSD(1:fitlength) );
 