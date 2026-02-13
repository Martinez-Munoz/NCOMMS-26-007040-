function diffusionCoefficients = calculateDiffusion(TRtotal,pixelSize,Dboundary,fitMode,varargin) 
% Description: Script to analyze the TRtotal matrix.
% Output: Dout, input/output files of msddGJB18 (histogram),  MeanMSDanalysisGJB18 (mean MSD and fit 
% data) and an output file containing all parameters of the individual trajectories (parameters from 
% linear fit and parameters like spot number, intensity, movie number, and more).
%                                                 
% Input: TRtotal.
%        pixelSize in microns
%        Dboundary in microns^2/s: maximum Diffusion coefficient for
%            immobile particles
%        fitMode = 'nofit', 'free', 'alpha', 'directed', 'confined'
% 
% The script contains a number of functions, to select and analyze subsets of the dataset and to 
% visualize the behavior of these.  These functions can be turned on or off by commenting and they can 
% be moved or repeated to achieve the desired perspective on the data. Selections can be made on 
% (minimum or maximum or in a window) trajectory lengths and trajectories can be truncated 
% (‘LengthAdjustAndSelect12’). The initial diffusion constants (Dini fit through the first 3 or 4 time 
% lags in the individual trajectory MSD plot) are calculated and put in a histogram (‘msddGJB17’). A 
% subset of this histogram or all trajectories can be selected for further analysis. Output: the script 
% generates a column of Dini of the selected trajectories (using DoftracesGJB1) which can be used to 
% combine or fit histograms in other functions (MultihistGJB1, FitDhistogramGJB1).
% 
% The following options are available for the analysis:
% msddGJB18: histogram of Dini of all individual selected trajectories. One can run the script 
% ‘AnalyseAllFilesGJB18_MAIN’once, a log file will be made of the histogram data. If one runs the script 
% again, one can set an option to plot the histogram of the first run as a line plot through the 
% histogram of the second run. Binning can be adjusted. 
% A log file with the data is generated in the work folder of Matlab.
%
% Dout=OutputGJB18: generates an output file containing all parameters of the individual trajectories 
% (parameters from linear fit and parameters like spot number, intensity, movie number, and more).
%
% watchdeltapeakGJB16: generate two histograms of all the one-time-lag-displacements (?x and ?y) of all 
% the selected trajectories together. The histogram is fitted with a normal function and the STD and 
% mean are given. Also an indication of the diffusion constant is given, but this value is biased by the 
% STD due to the accuracy.
%   
% MeanMSDanalysisGJB18: generates the mean square displacement plot of all the selected trajectories 
% together, with the number of used steps per time lag on the right y axis. The MSD plot can be fitted 
% with functions of models for simple diffusion, anomalous diffusion, directed diffusion and corralled / 
% confined diffusion. 
% A log file with the data is generated in the work folder of Matlab.
% 
% ScatteranalysisGJB18: generates scatter plots of trajectory length, spot intensity or spot radius of 
% gyration as a function of Dini. Also a scatter plot of ‘single trajectory mean square displacement 
% offset at the origin’ as a function of spot intensity is generated. The spot intensity and radius of 
% gyration are taken from the raw image stack. For the spot intensity, a camera background has been 
% subtracted before integrating all the intensities of all the pixels within the chosen spot radius. 
% 
% CumProbAnalysisGJB17: generates the cumulative probability distribution of all selected trajectories 
% together and fits this one with a function which represents two components of Brownian diffusion. From 
% the fit at different time lags, three plots are generated: one representing the fraction (alpha) of 
% diffusion component 1, two representing the square displacement of diffusion component 1 and three 
% representing the square displacement component 2. The two square displacement curves are fitted with a 
% linear function to estimate the diffusion constants of the different components.

if nargin<=4
    outputSuffix='';
else
    outputSuffix=varargin{1};
end

close all

r=size(TRtotal,1);
rr=size(TRtotal,2);
Tau=TRtotal(1,(rr-1)); % framerate, time between two frames

TRcont=zeros(r,(rr-2));         % put all trajectories of all movies in one matrix, 
TRcont(1,:)=TRtotal(1,1:(rr-2));% with increasing number.
TRcont(1,(rr-2))=1;
s=1;
for i=2:r
    if TRtotal(i,rr)==TRtotal(i-1,rr) && TRtotal(i,rr-2)==TRtotal(i-1,rr-2)
        TRcont(i,1:(rr-3))=TRtotal(i,1:(rr-3));
        TRcont(i,rr-2)=s;
    else
        s=s+1;
        TRcont(i,1:(rr-3))=TRtotal(i,1:(rr-3));
        TRcont(i,rr-2)=s;
    end;
end;

%%% -------------- analysis: trace selection part -----------------------------------------

% --- Discriminate on trajectory length, truncate trajectories if wanted. -----------------
selectlength=8; %21 17 13 8 select the traces which are 'selectlength' or longer, out of the dataset. 
% if 'selectlength' is equal to or smaller than tracelength, all trajectories will be selected.
truncate=2;
interactive5=0; % 0=all selected trajectories, 1=all but truncated, 2=only longer than, 
% 3=only shorter or equal to truncate
TRlength = LengthAdjustAndSelect12(TRcont, selectlength, truncate, interactive5);

% TRlength=TRtotal;
% ---- The individual trajectory MSD's are calculated and Ds are put in a histogram -------
% --- Split off immobile fraction ---



Fit=1; % Choose fit function, 1 is linear over first 3 or 4 points, 2 is A*t^alpha over first 1/4th of MSD points.  

       % Si eliges Fit=2, el software te da el alpha value por
       % trayectoria.Aparece en la columna 12. La columna 13 es el error
       % del alphavalue



interactive3=0;  % if interactive3==1, the individual fitted MSD curves are displayed.
bigorsmall=1; % if bigorsmall==1, all trajectories bigger or equal to Dboundary will be taken.
              % if bigorsmall== 0, all trajectories smaller to Dboundary will be taken.
              % if bigorsmall is not 1 or 0, the boundary will not be considered.
addformerhist=0; % if = 1, add the last histogram made, into this histogram in the form of a line. 
% The data of the former histogram is saved in a text and matlab file in the workspace.
% For ease of use the histogram data file in the matlab workspace will not be replaced if addformerhist = 1!
genoutput=0; % generate output files with trajectory fit data of ALL TRAJECTORIES IN INPUT msddGJB18(TRlength)!!
             % (not histogram data), 0 = no, 1 = yes.
bars=30;
[TRmobile,DiF] = msddGJB20(TRlength, pixelSize, Tau, Fit, Dboundary, interactive3 ,bigorsmall, addformerhist,genoutput, bars);

% COSS
diffusionCoefficients = DiF(:,1);
trajectoryIdx = unique(TRtotal(:,8));
fnMask=fullfile('mask.tif');
fnOut=fullfile('results','TrackingPackage','tracks',['diffusionCoefficients' outputSuffix '.txt']);
fnOutMobile=fullfile('results','TrackingPackage','tracks',['diffusionCoefficientsMobile' outputSuffix '.txt']);
if exist(fnMask,'file')
    mask=imread(fnMask);
    maskLabelVote = zeros(length(diffusionCoefficients),max(mask(:))+1);
    for j=1:size(TRlength,1)
        xj=round(TRlength(j,1));
        yj=round(TRlength(j,2));
        idx=round(TRlength(j,8));
        maskj=mask(yj,xj)+1;
        maskLabelVote(idx,maskj)=maskLabelVote(idx,maskj)+1;
    end
    maskLabel = zeros(length(diffusionCoefficients),1);
    fhOut = fopen(fnOut,'w');
    fhOutMobile = fopen(fnOutMobile,'w');
    for j=1:length(diffusionCoefficients)
        [dummy, maskLabel(j)] = max(maskLabelVote(j,:));
        fprintf(fhOut,'%d %f %d\n',trajectoryIdx(j),diffusionCoefficients(j),maskLabel(j)-1);
        if diffusionCoefficients(j)>Dboundary
            fprintf(fhOutMobile,'%d %f %d\n',trajectoryIdx(j),diffusionCoefficients(j),maskLabel(j)-1);
        end
    end
    fclose(fhOut);
    fclose(fhOutMobile);
else
    fhOut = fopen(fnOut,'w');
    fhOutMobile = fopen(fnOutMobile,'w');
    for j=1:length(diffusionCoefficients)
        fprintf(fhOut,'%d %f\n',trajectoryIdx(j),diffusionCoefficients(j));
        if diffusionCoefficients(j)>Dboundary
            fprintf(fhOutMobile,'%d %f\n',trajectoryIdx(j),diffusionCoefficients(j));
        end
    end
    fclose(fhOut);
    fclose(fhOutMobile);
end

disp(['Mean of the diffusion coefficient = ' num2str(mean(diffusionCoefficients))])
disp(['95% Percentile of the diffusion coefficient = ' num2str(quantile(diffusionCoefficients,.95))])

% TRmobile1 = msddGJB19(TRmobile, Tau, 1, 0.01, 0 ,0, 0,0, bars);
% TRmobile=TRmobile1;
% TRbar=MSDhistGJB20(TRlength, Tau, Fit, interactive3, addformerhist,genoutput, 10);

%  Fit=1; %1=straight line (fit 4 (3) first points), 2=t^alpha, 3=directed, 4=confined (1/4 of MSD points used for fits).
%  show=0; % show every trajectory one by one after eachother.
%  genoutput=0; % Generate output files if = 1, with all trajectory information, also real frame nr and movienumber.
% %  DiF2=OutputGJB20(TRmobile,TRtotal,Fit,genoutput,show); % must be one to generate output file used by ScatterAnalysis
%Dout=OutputGJB20xy(TRmobile,TRtotal,Fit,genoutput,show); % if fit is 1,
% this fct will give the x and y coordinates of first frame of a trajectory. Handy for
% colocalization of two-color data.

% %%% -------------- analysis part ----------------------------------------------------------
% % --- Watch the Displacement histogram of the followed peaks from one frame to the other --
% DT = watchdeltapeakGJB16(TRmobile);
% Ddeltax=DT/Tau
% Output0=DT/Tau; % The average D is calculated from the SD of a Gaussian fit to the histogram
% 


%% Average MSD





%%------ The Mean MSD is calculated out of all squared displacements, output D & plot ------
truncateMSD=8;

interactiveMSD=0;% 0=all selected trajectories, 1=all but truncated, 2=only longer than, 
% 3=only shorter or equal to truncate


plotlength=13; %% Puedes modificarlo! Timelag maximo que quieres

if strcmp(fitMode,'nofit')==1
    plotfit=0;
elseif strcmp(fitMode,'free')==1
    plotfit=1;
elseif strcmp(fitMode,'alpha')==1
    plotfit=2;
elseif strcmp(fitMode,'directed')==1
    plotfit=3;
elseif strcmp(fitMode,'confined')==1
    plotfit=4;
end
% plotfit=3; % 0=no fit through MSD, 1=straight line, 2=t^alpha, 3=directed, 4=confined.  %%Puedes modificarlo!


AvarD_MSD=MeanMSDanalysisGJB18(TRmobile,pixelSize,Tau,selectlength,truncateMSD,interactiveMSD,plotlength,plotfit);





%% Cumulative Probability Distribution   CPD!





% components=2; %%Puedes modificarlo!
% 
% fit1=2; % 0=no fit through MSD, 1=straight line, 2=t^alpha, 3=directed, 4=confined.  %%Puedes modificarlo!
% 
% fit2=2; % 0=no fit through MSD, 1=straight line, 2=t^alpha, 3=directed, 4=confined.  %%Puedes modificarlo!
% 
% plotlength=13;
% 
% [alfa1,alfa1error,D1,D1error,D2,D2error] =CumProbAnalysisGJB20(TRmobile,Tau,components,fit1,fit2,plotlength)
% % % % % % % 
% % % % % % % % Output1=Dout; %TRmobile; %
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

% Output file from OutputGJB20:

% Filename:
%     if Fit==1 % save linear fit parameters and trajectory information. Sorry all parameters are mixed.
%         Fname='OutLinearGJBASCII.txt'
%     elseif Fit==2 % Dini is replaced by D(1 sec), save also alpha and alplha error 
%         Fname='OutTalphaGJBASCII.txt'
%     elseif Fit==3 % Dini is replaced by D from directed fit save also V and V error 
%         Fname='OutDirectedGJBASCII.txt'
%     elseif Fit==4 % save linear fit (Dini) parameters and confined fit parameters as well. Goodness and degrees of freedom are only saved for the confined fit!
%         Fname='OutConfinedGJBASCII.txt'
%     end;

% Saved parameters:
%         DiF2(:,1)=D
%         DiF2(:,2)=tracenr, the one related to a movie!
%         DiF2(:,3)=error in D
%         DiF2(:,4)=Offset of linear fit through MSD at origin
%         DiF2(:,5)=DeltaMSDoffset; 
%         DiF2(:,6)=reduced Chi square of linear fit: Chi/degrees of freedom (0ne or two)
%         DiF2(:,7)=degrees of freedom: fitlength - nr. fitparameters
%         DiF2(:,8)=Tau;
%         DiF2(:,9)=number of frames of trajectory, including missing frames.
%         DiF2(:,10)= mean total-spot-intensity of the first 1 frames of the trace. 3,4: filtered, 5,6: from raw image, minus BG if Im right. ;
%         DiF2(:,11)=('r\^2')
%         DiF2(:,12)=the particular movie of the dataset TRtotal, where the trajectory belongs to.
%         if Fit==2 % in case single trajectories are fitted with Bt^alpha,
%             DiF(v,13)=alpha; % alpha fit parameter 
%             DiF(v,14)=Deltaalpha; %error in alpha fit parameter.
%         elseif Fit==3
%             DiF(v,13)=V; % V (speed, µm/s) fit parameter 
%             DiF(v,14)=DeltaV; %error in V fit parameter.
%         elseif Fit==4 % x(1) = 4BD; x(2) = A; x(3) = <Rc^2>.
%             DiF(v,13)=X(1)/4;    % fit parameter B*D, where D ~ Dini, DiF(v,1)
%             DiF(v,14)=errorX(1); % error in B*D
%             DiF(v,15)=X(2);      % fit parameter A
%             DiF(v,16)=errorX(2); % error in fit parameter A
%             DiF(v,17)=X(3);      % fit parameter <Rc^2>, related to size confinement zone
%             DiF(v,18)=errorX(3); % error in fit parameter <Rc^2>
%         end;



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

% Output file of MeanMSDanalysisGJB18:
% Colums 'tlag','MeanMSD','ErrorMSD','TmsdFit' and 'FitMSD'.



% % ----------- Remove traces ------------------------------ --------------------------------
% wiped='22 26 30 72 77 102 131 132 151 154 179 188 200 204 218 248'; % convert cell output into string
% TRwiped=removetracesGJB16(TRmobile,wiped,1,1,1);
% Output1=TRwiped;
% dummy = msddGJB18(TRwiped, Tau, Dboundary, 0,2,0,0, bars);

