function F = FitMeanMSDlinearGJB18(x, SigmaMSD, Tmsd, AvarMSD, fitlength);  
% (similar as FitFunctionGJB10c, a linear fit through the MSD points)

% weighted fit: minimize weightfactor*( 4D*time+offsetAtX=0 - MSD) (interpret as Ax-b or f(xi)-yi from yi~f(xi))
% x(1) = 4D;
% x(2) = offset;

% fitlength=10; %%not original
 
F=(1./SigmaMSD(1:fitlength)).*(Tmsd(1:fitlength).*x(1) + x(2) -AvarMSD(1:fitlength));
 
