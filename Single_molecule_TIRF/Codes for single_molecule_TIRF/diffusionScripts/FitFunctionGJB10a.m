function F = FitFunctionGJB10a(x,MSD, MSDpointsFit);  % weighted fit: minimize weightfactor*( time*x+offsetAtX=0 - MSD) (interpret as Ax-b or f(xi)-yi from yi~f(xi))

% Only take first 1/4th of points to be more accurate, for longer traces
% only take first 3 to avoid influence due to confinement as much as
% possible. 

% Matr(:,1)=MSD(1:MSDpoints,1);
% Matr(:,2)=ones(MSDpoints,1).*MSDoffset./x;
% varia=[x;x];


F=(1./MSD(1:MSDpointsFit,3)).*(MSD(1:MSDpointsFit,1).*x(1) + x(2) - MSD(1:MSDpointsFit,2));
% F=((MSD(1:MSDpointsFit,1).*x(1) + x(2)) - MSD(1:MSDpointsFit,2));