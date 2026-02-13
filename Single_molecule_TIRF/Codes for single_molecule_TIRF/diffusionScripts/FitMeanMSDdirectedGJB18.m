function F = FitMeanMSDdirectedGJB18(x, SigmaMSD, Tmsd, AvarMSD, fitlength);  
% weighted fit: minimize: weightfactor*( 4D*time + (V*time)^2 + offsetAtX=0  -  MSD)
% x(1) = 4D;
% x(2) = offset;
% x(3) = V.

F=(1./SigmaMSD(1:fitlength)).*( Tmsd(1:fitlength).*x(1) + (Tmsd(1:fitlength).*x(3)).^2 + x(2) - AvarMSD(1:fitlength) );
 