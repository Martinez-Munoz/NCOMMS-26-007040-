function [TRmobile,DiF] = msddGJB20(TR, pixelSize, Tau, Fit, Dboundary, interactive3, bigorsmall,addformerhist,genoutput, n)
% Calculates the MSD for the points (not all, see below) of a trace.
% Calculates the D out of the weighted fit for every trace.
% Plots the histogram of all D's from all the traces.
% Generates datafiles with trajectory fit parameters and histogram data.

% If interactive3 ==1, the fitted MSD will be plotted for every trace one
% after eachother. Otherwise only the average MSD and the histogram will be
% plotted.

% Input: traced points, time between frames (Tau), interactive3 and n, number of bars in histogram.

% 28-3-06: Big Error removed, D of traces after a whiped trace were counted
% multiple to the amount of wiped trajectories just before!! Result: shape
% of histogram more or less the same, but height different.

% 30-9-07: Error removed, the program didn't take into account that some
% trajectories have frames missing. They go e.g. from movieframe x to x+2 due
% to for instance blinking. The tracking program tolerates one missing 
% frame (parameter set to 1), but the building of the MSD didn't take it 
% into account! How could I've missed this for such a long time!! Well, at 
% least I've discovered and corrected it. Now I doublechecked the whole
% algorithm and it works fine with normal trajectories and trajectories
% where one timelag at the time is skipped. MORE THAN ONE FRAME MISSED ONE
% AFTER EACHOTHER IS NOT CORRECTED! And take care, if a short trace has many
% missing frames ( e.g. 1, missed, 3, missed,  5, missed, 7, 8, 9, 10, 11) 
% than the statistics of the MSD becomes poor ( # squared displacements for 
% timelag 1 upto 10: 4, 6, 3, 4, 2, 3, 2, 2, 1, 1). (VERSION 16)

% 13-1-08: error removed, that made trajectories of length 20 up to
% 23 being fit to 5 points. Trajectories with a length of 16 respectively 12 were fitted
% with 4 repectively 3 points, I made it now with 3 points, since a trajectory has length-1  
% MSD points and it is better to take 1/4th of the MSD points to get a good
% fit. Also confusing code removed, by separating column with all MSD
% points of all traces from the column with all the diffusion constants.
% The calculation of the error in the diffusion constants has been improved
% and more output is generated in text files: now the text file contains D,
% trajectory nr, Delta D, Offset, Delta offset and goodness of the fit.(VERSION 18)

% Output file of msddGJB19:
%       A  DiF2(:,1)=D
%       B  DiF2(:,2)=tracenr, the one of the list of trajectories without the movienumbers (used in this script, 1x8 matrix)!
%       C  DiF2(:,3)=error in D
%       D  DiF2(:,4)=Offset of linear fit through MSD at origin
%       E  DiF2(:,5)=DeltaMSDoffset; 
%       F  DiF2(:,6)=reduced Chi square of linear fit: Chi/degrees of freedom (0ne or two)
%       G  DiF2(:,7)=degrees of freedom: fitlength - nr. fitparameters
%       H  DiF2(:,8)=Tau;
%       I  DiF2(:,9)=number of frames of trajectory, including missing frames.
%       J  DiF2(:,10)= mean total-spot-intensity of the first 1 frames of the trace. 3,4: filtered, 5,6: from raw image, minus BG if Im right. ;
%       K  DiF2(:,11)=('r\^2')
%         if Fit==2 % in case single trajectories are fitted with Bt^alpha,
%       L      DiF(v,12)=alpha; % alpha fit parameter 
%       M      DiF(v,13)=Deltaalpha; %error in alpha fit parameter.
%         end;

if interactive3==1 % make a new frame if interactive3=1
    figure
end;

r=size(TR,1);
rr=size(TR,2);  
% NrOfTrs=TR(r,rr); %this is not the number of trajectories, but the biggest trajectory number.

MSD_D=zeros(1,4); %(r,4)
DiF=zeros(1,6);
u=1;
v=1;

mecmec=unique(TR(:,8));

for j=1:numel(mecmec) %NrOfTrs      check all trajectory numbers, if they are in the dataset
    in_inter=mecmec(j);
    tracenr=j;
    trace=zeros(1,rr);
    
    s=0;
    for i=1:r
        if TR(i,rr)==tracenr % replace zero vector 'trace' if the trajectory number 'tracenr' is present in the TR dataset.
            trace(s+1,:)=TR(i,:);
            s=s+1;
        end;
    end; 
       
    %calc MSD, D
    if trace(1,:)~=zeros(1,rr) % only calculate MSD if 'trace' is not just a zero vector.
        tracelength=size(trace,1);
        mol=zeros(tracelength,2);
        mol(:,1)=pixelSize*(trace(:,1)-min(trace(:,1)));		%X-coördinaten omreken naar µmeters en relative nul geven
        mol(:,2)=pixelSize*(trace(:,2)-min(trace(:,2)));		%Idem Y-coördinaten
        
         % ----------------Start corrected part to remove artifacts due to skipped frames -------------- 
        
        % 30-9-07: Check if the trace did skip one movieframe (e.g. due to
        % blinking). More than one skipped frame is not implemented here!
        % If a frame number is missing in the trace, it will be added in 
        % such a way that the 'missed frames' doesn't create the artifact that 
        % a missing frame and the next one will count for a single time lag. The squared displacements 
        % of the added frames will be recognized because of unrealistically big or small values andwill be removed.   
        molAA=zeros(1,2);
        molBB=zeros(1,2);        
        molAA(1,:)=mol(1,:); % iniciate the two matrices that can be substracted from eachother with a certain timelag.
        molBB(1,:)=mol(1,:); 
        s=0;
        for i=1:tracelength-1 % Check the whole trace
            if trace(i+1,end-1) ~= trace(i,end-1)+1  %Did the trace skip a movieframe?
                molAA(i+1+s,:)=[1000 1000]; % add the missing frame, with unrealistic x and y
                molBB(i+1+s,:)=[500 500]; % if lag is such that a squared displacement is an added - added frame, the value will still be unrealistic: 2000 - 1000 = 1000
                s=s+1; % the new trace index becomes bigger by 's' due to the adding of frames. 
                molAA(i+1+s,:)=mol(i+1,:); % give the frame after the added frame a higher index number and go to the next i
                molBB(i+1+s,:)=mol(i+1,:);
            else
                molAA(i+1+s,:)=mol(i+1,:);
                molBB(i+1+s,:)=mol(i+1,:);
                s=s; % the new trace becomes bigger by 's' due to the adding of s extra frames before this one. 
            end;
        end;           
        
        newtracelength=size(molAA,1);
        MSD=zeros(1,3);
        
        for i=1:(newtracelength - 1)                             %MSD punt 1 (X2-X1, ... Xn-Xn-1) tot en met n-1 (Xn-X1) van de trace met lengte n.
            DISP=zeros(1,2);
            molBB=[0 0											%Posities op i-1
                molBB];
            t=size(molAA);										%bepalen laatste meting;
            molAA=molAA(2:t(1,1),:);                            %telkens een MSD waarde minder bij een grotere stap door met 2 te beginnen;
            molBB=molBB(2:t(1,1),:);	    					%colommen even lang maken door B inkorten door lengte A te nemen;
            DISP=(molAA-molBB).^2;							    %squared displacement per coördinaat
            
            s=1;                 % remove the squared displacements which are unrealistically big, those are the ones based on a skipped frame.
            disp=zeros(1,2);
            for h=1:t(1,1)-1
                if DISP(h,1)<90000
                    disp(s,:)=DISP(h,:);
                    s=s+1;
                end;
            end;
           
            MSD(i,2)=mean(disp(:,1)+disp(:,2));				    %MSD 
            MSD(i,1)=i*Tau;                                     %Time
            MSD(i,3)=std(disp(:,1)+disp(:,2))/sqrt(t(1,1)-1);   %Error = std/sqrt(number of MSDs used for the average MSD(t) of one trace)
        end;
        
        MSDpoints=newtracelength-1; 
        
%         if MSDpoints<16  % e.g. (15)/4=3.75 >> fix >> =3
%              MSDpointsFit=fix((MSDpoints)/4); %(number of MSD points in trajectory)/4, rounded to zero.
%         else
%             MSDpointsFit=4; % this makes most sense for Tau = 0.05. With 6 points fitted MSD, 
%         end;                % the error in D becomes errorD=MSDoffset/4*Timelag=0.01/4*0.05*6=0.008 µm2/s
%                             % and this is for me an acceptable error. If Tau=0.1, only 3 points are needed.
%                             % >>> think about this. The error is probably
%                         % That's why I made it 4 instead of 6.
%                             % 4 turned out to be good, to cancel out
%                             % 'oscillations in the MSD plot as well. Less
%                             % gives incorrect fits for low D.
         MSDpointsFit=4;                   
        % ---------------- End corrected part to remove artifacts due to skipped frames --------------              

        % MSDoffset is one of the two fitparameters. Other option: fixed
        % offset calculated from spot intensity. See below, doesnt work well
        % for (fast or slow) moving spots on cells. Due to Ifluctuations caused by rotation?       
        % weighted fit: minimize weightfactor*( time*x-offsetAtX=0 - MSD) (interpret as Ax-b or f(xi)-yi from yi~f(xi))
        % starts at the point x0 and finds a minimum of the sum of squares of the functions 
        % described in fun. fun should return a vector of values and not the sum of squares 
        % of the values. (The algorithm implicitly sums and squares fun(x).)
        
        if Fit==1
            fitlength=4; %6
            x0=[0.1,0.005]; % x(1) = 4D; x(2) = offset;
            [X,resnorm,residual,eflag,out,lambda,jaco] = lsqnonlin(@FitFunctionGJB10a,x0,[0,0.0001],[40,0.05],optimset('display','off'),MSD, MSDpointsFit);%(@FitMeanMSDlinearGJB18,x0,[0,0],[40,0.05],[], SigmaMSD, Tmsd, AvarMSD, fitlength);
            TmsdFit=[0;MSD(1:round(fitlength),1)];
            FitMSD=TmsdFit.*X(1) + X(2);
            XiQuadr=resnorm; %Xi=sum[(fit-obs)/std]^2 where std is integrated into the fitfunction, as the weight factor.
            dgrfreedom=fitlength-2; % 2 fitparameters
            goodness=XiQuadr/dgrfreedom;
            errorX=reshape(full(sqrt(goodness*diag(inv(jaco'*jaco))))',size(X) );   
            x=X(1);
            Deltax=errorX(1);
            MSDoffset=X(2);
            DeltaMSDoffset=errorX(2);
        elseif Fit==2
            fitlength=fix((MSDpoints)/4); %(number of MSD points in trajectory)/4, rounded to zero.
            x0=[0.1,0.005,1]; % x(1) = B = 4D at timelag = 1 sec; x(2) = offset; x(3) = alpha.
            [X,resnorm,residual,eflag,out,lambda,jaco] = lsqnonlin(@FitMeanMSDtalphaGJB18,x0,[0,0,0],[40,0.05,5],optimset('display','off'), MSD(:,3), MSD(:,1), MSD(:,2), fitlength);
            TmsdFit=[0:(0.05*Tau):((fitlength+0.5)*Tau)]; 
            FitMSD=(TmsdFit.^X(3)).*X(1) + X(2);
            XiQuadr=resnorm; %Xi=sum[(fit-obs)/std]^2 
            dgrfreedom=fitlength-3; % 3 fitparameters
            goodness=XiQuadr/dgrfreedom; % reduced Chi^2, 'XiQuadr'.
            errorX=reshape(full(sqrt(goodness*diag(inv(jaco'*jaco))))',size(X) );
            x=X(1);
            Deltax=errorX(1);
            MSDoffset=X(2);
            DeltaMSDoffset=errorX(2);
            alpha=X(3);
            Deltaalpha=errorX(3);
        end;
        
        %built column T, MSD, deltaMSD, Tracenr per trace. Each trace will be added underneath to
        %the same colums.
        MSD_D(u:u+MSDpoints-1,1)=MSD(1:MSDpoints,1); %time
        MSD_D(u:u+MSDpoints-1,2)=MSD(1:MSDpoints,2); %MSD
        MSD_D(u:u+MSDpoints-1,3)=MSD(1:MSDpoints,3); %deltaMSD
        MSD_D(u:u+MSDpoints-1,4)=tracenr;%tracenr
        u=u+MSDpoints;
        
        %buit column D, trajectory number, Delta D, offset, Delta offset,
        %goodness of the fit, degrees of freedom, Tau and trajectory length:
        %Sorry all parameters are mixed. It is due to the history of the
        %program.
        DiF(v,1)=x/4; %D
        DiF(v,2)=in_inter; % trajectory nr in TR, not in relation to movienr.
        DiF(v,3)=Deltax/4; %error in D
        DiF(v,4)=MSDoffset; % Offset of linear fit through MSD at origin
        DiF(v,5)=DeltaMSDoffset; 
        DiF(v,6)=goodness; %reduced Chi square of linear fit: Chi/degrees of freedom (0ne or two)
        DiF(v,7)=dgrfreedom; %degrees of freedom: fitlength - nr. fitparameters
        DiF(v,8)=Tau;
        DiF(v,9)=newtracelength; %number of frames of trajectory, including missing frames.
        DiF(v,10)=mean(trace(1,5));    % mean total-spot-intensity of the first 1 frames of the trace. 3,4: filtered, 5,6: from raw image, minus BG if Im right. ;
        DiF(v,11)=mean(trace(1,6));    % mean r^2 of the first 1 frames of the trace. 3,4: filtered, 5,6: from raw image, minus BG if Im right. ;
        if Fit==2 % in case single trajectories are fitted with Bt^alpha,
            DiF(v,12)=alpha; % alpha fit parameter 
            DiF(v,13)=Deltaalpha; %error in alpha fit parameter.
        end;
        v=v+1;
        
        %option to plot MSD with fit, per trajectory.
        if interactive3==1
            errorbar(MSD(1:MSDpoints,1),MSD(1:MSDpoints,2),MSD(1:MSDpoints,3),'o');
            hold on
            plot(MSD(1:MSDpoints,1),MSD(1:MSDpoints,2),'*',TmsdFit,FitMSD,'-','LineWidth',2);
            ymax=2*max(MSD(1:MSDpoints,2));
            axis([0 1.5*MSD(MSDpoints,1) 0 ymax]);
            text(0.2*1.5*MSD(MSDpoints,1),0.8*ymax, [int2str(tracenr)] ,'HorizontalAlignment','left','FontSize',16, 'color',[0.7 0.7 0.7], 'FontWeight', 'bold');   %tracenr
            text(0.2*1.5*MSD(MSDpoints,1),0.7*ymax, ['D = ',num2str(x/4), ' \pm ',num2str(Deltax/4)] ,'HorizontalAlignment','left','FontSize',16, 'color',[0.7 0.7 0.7], 'FontWeight', 'bold');   %D
            text(0.2*1.5*MSD(MSDpoints,1),0.6*ymax, ['Offset = ',num2str(MSDoffset), ' \pm ',num2str(DeltaMSDoffset/4)], 'HorizontalAlignment','left','FontSize',16, 'color',[0.7 0.7 0.7], 'FontWeight', 'bold');   %D
            hold off
            pause
        end;
    end; % end only calculate MSD when trace is not just a zero vector.
end; % end check all trajectory numbers, if they are in the dataset

% Save linear fit data of the individual trajectories. A file will be made
% (or over writtten if present) in the work folder. One file is in .m, the
% other in ascii. If necessary, one can also save the MSD points of the
% individual trajectories (column MSD_D). The data in the files is
% separated by tabs and columns are seperated by a square symbol.
% The columns are: D, trajectory nr, error D, offset, error in offset,
% goodness (reduced Chi squared of the lin fit), degrees of freedom, Tau,
% number of frames in trajectory (missing frames included).


if genoutput==1
    if Fit==1 % save all fit parameters and trajectory information. Sorry all parameters are mixed.
        DataFile=[DiF(:,1)'; DiF(:,2)'; DiF(:,3)'; DiF(:,4)'; DiF(:,5)'; DiF(:,6)'; DiF(:,7)'; DiF(:,8)'; DiF(:,9)'; DiF(:,10)'; DiF(:,11)']';
    elseif Fit==2 % save also alpha and alplha error.
        DataFile=[DiF(:,1)'; DiF(:,2)'; DiF(:,3)'; DiF(:,4)'; DiF(:,5)'; DiF(:,6)'; DiF(:,7)'; DiF(:,8)'; DiF(:,9)'; DiF(:,10)'; DiF(:,11)'; DiF(:,12)'; DiF(:,13)']';
    end;
    DiF2=DiF';
    save('msddGJBout','DiF2'); 
    save('msddGJBoutASCII.txt','DataFile', '-ASCII' , '-tabs');   
end;

% Make the plot:
nrtr=size(DiF,1); % number of diffusion constants, traces

figure
if min( abs(DiF(:,1)) )<0.001
    La=-3;
    A=0.001;
    LA=log10(A);
else
    a=min( abs(DiF(:,1)) );
    La=floor(log10(a));
    LA=log10(1*10^La);
end;
if max( abs(DiF(:,1)) )>1
    b=max( abs(DiF(:,1)) );
    Lb=ceil(log10(b));
    LB=log10(1*10^Lb);
else
    Lb = 0;
    B=1;
    LB=log10(B);
end;

y = logspace(La,Lb,n); %generates n points between decades 10^a and 10^b.
yy = logspace(LA,LB,(LB-LA+1));


hold on
if addformerhist==1   % load histogram file from workspace and plot
    load('datalasthistogram','xout2','N2');
    hist(DiF(:,1),y);
    [N,xout]=hist(DiF(:,1),y);
    stairs(xout2/10^(1.3/n),(N2/N2(1))*N(1),'-r'); 
else % save histogram file when the former one  won't be used in this histogram
    DiF2=DiF;
    y2=y;
    [N2,xout2]=hist(DiF(:,1),y);
    save('datalasthistogram','N2','xout2','nrtr'); 
%     aa=['Xbins','Ycounts'];
%     DataText=cellstr(aa)';
    DataFile=[xout2; N2]';
    save('datalasthistogramASCII.txt','DataFile', '-ASCII' , '-tabs');     %,'DataText'
    hist(DiF(:,1),y);
end;
%XScale      {linear} | log
set(gca,'XScale','log');
% set(gca,'XTick',yy)  % Ive to improve this!!!!!!!!!!!!!!!!!!!!!
% set(gca,'XTickLabel',yy)
set(gca,'XLim',[0.8*10^(LA) 10^(LB)])
Y=get(gca,'YLim');
set(gca,'YLim',[Y(1) Y(2)*1.1])
text(2*10^(LA),0.7*Y(2), [int2str(nrtr), ' trajectories'] ,'HorizontalAlignment','left','FontSize',16, 'color',[0.7 0.7 0.7], 'FontWeight', 'bold');   % # of traces
Largestbar=max(N2);
LineXdboundary=ones(2)*Dboundary;
LineYdboundary=[0 Largestbar*0.5];
line(LineXdboundary,LineYdboundary);
hold off
title('Diffusion coefficients');
xlabel('D (µm^2/sec)');
ylabel('Number of trajectories');

if Fit==2
    figure
    yyy=linspace(0,2,n);
    hist(DiF(:,12),yyy) %histogram of alpha
end;


%%% --------------- select traces according to Dboundary ------------------

s=1;
immobile=0;
sizeDiF=size(DiF,1);

if bigorsmall==1 % only trajectories bigger or equal  
    for i=1:sizeDiF
        if DiF(i,1)>=Dboundary
            selected(s)=DiF(i,2);
            s=s+1;
        else
            immobile=immobile+1;
        end;
    end;

    s=1;   
    sizew=size(selected,2);
    for j=1:sizew;                     % for every selected trace,
        for i=1:r                      % pick out the TR rows of the paticular trace
            if TR(i,rr)==selected(j)
                TRmobile(s,:)=TR(i,:);
                s=s+1;
            end;
        end;
    end;
    
elseif bigorsmall==0 %only trajectories smaller than
    for i=1:sizeDiF
        if DiF(i,1)<Dboundary
            selected(s)=DiF(i,2);
            s=s+1;
        else
            immobile=immobile+1;
        end;
    end;

    s=1;   
    sizew=size(selected,2);
    for j=1:sizew;                     % for every selected trace,
        for i=1:r                      % pick out the TR rows of the paticular trace
            if TR(i,rr)==selected(j)
                TRmobile(s,:)=TR(i,:);
                s=s+1;
            end;
        end;
    end;
else  % all trajectories if bigorsmall =~ 1 or 0. 
    TRmobile=TR;
end
    

% % 
% % % Some discussion:
% % % The size of MSDpoints is size(trace,1)-4, because there are size(trace,1)-1 MSD
% %         % points, the last MSD value has no calculated error and can thus not be weighted
% %         % in the average and the error in the other two values seems to be to much random.
% %         % A better error (weight factor) should be found, which takes the randomness of
% %         % STD(Square Displacements) into account!! Or is it not possible to
% %         % filter out statistically irrelevant datapoints automatically? In
% %         % that case, I should just take only the first 1/4th of the points
% %         % for a fit of a single trajectory. This is what I do now, see
% %         % FitfunctionGJB5. 
% %         % I still think that for the weighted average of the MSDs of all!! 
% %         % the individual trjectories together, a better weightfactor should 
% %         % be found. Now some points are more important due to a by coincidence 
% %         % smaller STD. Another way to solve this problem, is to connect all 
% %         % trajectories together and regard it as one long trajectory, and then 
% %         % use only 1/4 th of the total nr of points. In both cases one
% %         % assumes simple diffusion in all trajectories. But, in case of the
% %         % weighted average, additional information is still present,
% %         % because there are more short trajectories than long ones and
% %         % these also have other characteristics which can be found back in
        % the average MSD plot by means of a bent at the place where most 
        % short trajectories end.
        % At the moment the weighted fit (FitFunction2GJB5(x)) in the weighted 
        % average MSD plot of all trajectories, uses 1/4th of the points of
        % the total MSD points of all traces together, OR the miniumum tracelength. 
        % The smallest value will be chosen. Keeping in mind
        % that we can make one single trajectory of all trajectories and
        % than taking 1/4th of the MSD points, this number is ok as long as we dont have 
        % many short trajectories. Then we better chose the minimum
        % trajectory length from trackpeaksGJB5. This will be a good number
        % (maybe too small) as long as there are more than 4 trajectories.
        
        % here I make the MSD fit offset per trajectory dependent on the mean spot 
        % intensity of that particular trace. Offset ~ STDx^2 ~ (1/snr)^2 ~ (Ibg + Isignal)/Isignal^2   
%        Intensity=mean(trace(:,3));   
%        MSDoffset=0.80*((Intensity+500)/Intensity^2); %  ~(N/S)^2 = (STDmolecule^2 + STDbg^2) / Imolecule^2 
%        Intensityraw=mean(trace(:,5));
%        MSDoffset=1.50*((Intensityraw)/(Intensityraw-5000)^2);
        
%          Fixed MSDoffset, determined from
%          intensity, see above. Unfortunately this doesn't work well for
%          moving spots on cells. It is fine for fixed spots. Why? Not
%          due to smear-out effect and not due to setup accuracy. I don't
%          know.
%          x0=0.1;             
%          [x,resnorm] = lsqnonlin(@FitFunctionGJB10b,x0,[],[],[],MSD, MSDpointsFit, MSDoffset); 



% ------------------------------------ make histogram of Ds of immobile
% trajectories, where the fit is not restricted to D=0. There will be
% positive and negative Ds around the mean value of zero. -------------

% [X,resnorm] = lsqnonlin(@FitFunctionGJB10a,x0,[-40,0.001],[40,0.05],[],MSD, MSDpointsFit);
% Here the -40 is the direction coefficient A in Ax + B, which is negative now.


% shift immobile Ds to see how a distribution will be of almost not moving
% trajectories (part of understanding D distribution immobile spots):
% for i=1:nrtr
%     DiF(i,1)=DiF(i,1)+0.004;
% end;


% figure
% 
% [N3,xout3]=hist(DiF(:,1),15);
% save('datalasthistogramASCII2.txt','N3','xout3', '-ASCII' , '-tabs');
% hist(DiF(:,1),15);
% 
% % resultfitD = fit_ML_normal(DiF(:,1)); % find fit parameters
% % stdfit=sqrt(resultfitD.sig2);    % standard deviation, sigma of fit [µm^2/s]
% %nn=max(hist(N3));
% %set(gca,'YLim',[0 1.1*nn])
% %hAx=gca;%[ 0 1.1*nn];
% %plot_normal2( XX,resultfitD,hAx,1 )
% 
% X0=[1,1];
% coeff=fmins('gafit',X0,[],[],xout3,N3);
% a=coeff(1); 
% b=coeff(2);
% XX= min(xout3):(max(xout3)-min(xout3))/100:max(xout3); 
% y=a*exp(-b*XX.^2);
% sig=1/sqrt(b);
% hold on
% hist(DiF(:,1),15);
% plot(XX,y);
% text(0.8*min(xout3),0.8*max(N3), ['Sigma = ', num2str(sig), ' µm\^2 / sec'] ,'HorizontalAlignment','left','FontSize',12, 'color',[0 0 0], 'FontWeight', 'bold');  
% hold off



% -----------------------------------------------end improvisationary part