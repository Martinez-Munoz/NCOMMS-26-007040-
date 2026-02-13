function F = FitMeanMSDconfinedGJB18(x, SigmaMSD, Tmsd, AvarMSD,offset, fitlength);  
% weighted fit: minimize: weightfactor*( <Rc^2>[1-Aexp(-4BDt/<Rc^2>)] + offsetAtX=0  -  MSD)
% x(1) = 4BD;
% x(2) = A;
% x(3) = <Rc^2>.

 F=(1./SigmaMSD(1:fitlength)).*(x(3).*(1-x(2).*exp(-(x(1)/x(3)).*Tmsd(1:fitlength)) )  + offset - AvarMSD(1:fitlength) );
 
% F=(1./SigmaMSD(5:fitlength)).*(x(3).*(1-x(2).*exp(-(x(1)/x(3)).*Tmsd(5:fitlength)) )  + offset - AvarMSD(5:fitlength) );
 