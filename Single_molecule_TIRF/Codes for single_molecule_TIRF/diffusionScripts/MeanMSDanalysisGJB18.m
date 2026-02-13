function output = MeanMSDanalysisGJB18(TR, pixelSize, Tau, selectlength, truncate, interactive5, plotlength, plotfit)
% Plots the weighted mean of all mean square displacements from all the traces. 
% This mean MSD is ...

% Input: traced points, time between frames (Tau)
% Output: D average out of mean MSD fit.

r=size(TR,1);
rr=size(TR,2); 
NrOfTrs=TR(r,rr); %this is not the number of trajectories, but the biggest trajectory number.
%truncate trajectories to make them all the same length, so all
%trajectories will have the same weight in the mean MSD curve. If longer
%trajectories have different diffusion behaviour than shorter ones, the
%longer won't get dominant, which would normally be the case.


% take the trajectories which have a length 'selectlength' or longer, to work with.
u=1;
w=1;
for j=1:NrOfTrs                    % for every trace, 
         tracenr=j;           
         for i=1:r                      % pick out the TR rows of the paticular trace
             if TR(i,rr)==tracenr
                 if  (u)==selectlength
                     selected(w)=tracenr;  % a list of trajectories which have selectlength, for next option
                     w=w+1;
                 end;
                 u=u+1;
             else 
                 u=1;
             end;
         end;   
 end;
 
 s=1;   
 sizew=size(selected,2);
 for j=1:sizew;                      % for every trace longer / equal than selectlength,   
     for i=1:r                      % pick out the TR rows of the paticular trace
         if TR(i,rr)==selected(j)
             TRselect(s,:)=TR(i,:);
             s=s+1;
         end;
     end;   
 end;
TR=TRselect;   % work further with the selected trajectories, without looking at the other ones. 
r=size(TR,1);  % Now there are no shorter trajectories in the dataset and a good analysis can be 
rr=size(TR,2); % made byond the tracelength of the tracking algorithm itself. There will be no 
% bias due to the smaller amount of traces with longer than 'tracelength' lengths. Now we can 
% safely put the truncate higher than the 'tracelength' of the algorithm. We put the truncate 
% simply as high as the 'selectlength', but not further.

% second option: truncate every trajectory at timepoint 'truncate'.
u=1;
s=1;
w=1;
for j=1:NrOfTrs                    % for every trace, 
         tracenr=j;           
         for i=1:r                      % pick out the TR rows of the paticular trace
             if TR(i,rr)==tracenr
                 if u<=truncate
                     TRtrunc(s,:)=TR(i,:);
                     s=s+1;
                 elseif (u-1)==truncate
                     l(w)=tracenr;  % a list of trajectories which are truncated, for next option
                     w=w+1;
                 end;
                 u=u+1;
             else 
                 u=1;
             end;
         end;   
 end;
 
 % third option: take only trajectories longer than 'truncate'. 
 if interactive5>=2 
     s=1;   
     sizel=size(l,2); % l is the small letter L!!
     TRsubstract=zeros(r,rr);
     for j=1:sizel                    % for every trace longer than truncate,   
         for i=1:r                      % pick out the TR rows of the paticular trace
             if TR(i,rr)==l(j)
                 TRlong(s,:)=TR(i,:);
                 TRsubstract(i,:)=TR(i,:);
                 s=s+1;
             end;
         end;   
     end;
 end;

 % fourth option: take only trajectories shorter than or equal to 'truncate'. 
 if interactive5==3 
     s=1;
     TRshortzeros=TR-TRsubstract; % get only traces shorter than..
     for i=1:size(TRshortzeros,1)  % remove zeros for cosmetics
         if TRshortzeros(i,rr)~=0
             TRshort(s,:)=TRshortzeros(i,:);
             s=s+1;
         end;
     end;
 end;

 if interactive5==0 % no truncation, all trajectories fully analysed.
     TRtobedone=TR;
 elseif interactive5==1 % truncate trajectories to 'truncate', preferably to min. tracelength
     TRtobedone=TRtrunc;
 elseif interactive5==2 % only take trajectories longer than 'truncate'. 
     TRtobedone=TRlong;
 elseif interactive5==3 % only take trajectories shorter than 'truncate'.
     TRtobedone=TRshort;
 end;
 

for i=1:plotlength  %  msd plot length
    lag=i;
    tlag(i)=lag*Tau;
    dtotal=0;
    dtotal = MeanMSD_GJB16(TRtobedone,pixelSize,lag); %derive all the squared displacements for a certain timelag, of all given trajectories in TR.
    MeanMSD(i)=mean(dtotal);
    SizeMSD(i)=size(dtotal,1);
    ErrorMSD(i)=std(dtotal)/sqrt(SizeMSD(i));
end

%make this fit without global, input parameters via lsqnonlin!
SigmaMSD=ErrorMSD;
Tmsd=tlag; 
AvarMSD=MeanMSD; 


if plotfit<=1 %>>>>FitFunctionGJB10c
    fitlength=4; %6
    x0=[0.1,0.005]; % x(1) = 4D; x(2) = offset;
    [X,resnorm,residual,eflag,out,lambda,jaco] = lsqnonlin(@FitMeanMSDlinearGJB18,x0,[-40,0],[40,0.05],optimset('display','off'), SigmaMSD, Tmsd, AvarMSD, fitlength);
    TmsdFit=[0,Tmsd];   % add a point at 0 for the fit line 
    FitMSD=TmsdFit.*X(1) + X(2);
    XiQuadr=resnorm; %Xi=sum[(fit-obs)/std]^2 where std is integrated into the fitfunction, as the weight factor.
    dgrfreedom=fitlength-2; % 2 fitparameters
    goodness=XiQuadr/dgrfreedom;
    errorX=reshape(full(sqrt(goodness*diag(inv(jaco'*jaco))))',size(X) );
    output= [X(1)/4 errorX(1)/4];
elseif plotfit==2
    fitlength=plotlength; 
    x0=[0.1,0.005,1]; % x(1) = 4D; x(2) = offset; x(3) = alpha.
    [X,resnorm,residual,eflag,out,lambda,jaco] = lsqnonlin(@FitMeanMSDtalphaGJB18,x0,[-40,0,0],[40,0.05,5],optimset('display','off'), SigmaMSD, Tmsd, AvarMSD, fitlength);
    TmsdFit=[0,Tmsd];   % add a point at 0 for the fit line
    FitMSD=(TmsdFit.^X(3)).*X(1) + X(2);
    XiQuadr=resnorm; %Xi=sum[(fit-obs)/std]^2 
    dgrfreedom=fitlength-3; % 3 fitparameters
    goodness=XiQuadr/dgrfreedom; % reduced Chi^2, 'XiQuadr'.
    errorX=reshape(full(sqrt(goodness*diag(inv(jaco'*jaco))))',size(X) );
    output= [X(1)/4 errorX(1)/4];    
elseif plotfit==3
    fitlength=plotlength;
    x0=[0.1,0.005,0]; % x(1) = 4D; x(2) = offset; x(3) = V.
    [X,resnorm,residual,eflag,out,lambda,jaco] = lsqnonlin(@FitMeanMSDdirectedGJB18,x0,[-40,0,0],[40,0.05,1],optimset('display','off'), SigmaMSD, Tmsd, AvarMSD, fitlength);
    TmsdFit=[0,Tmsd];   % add a point at 0 for the fit line
    FitMSD=TmsdFit.*X(1) + (TmsdFit.*X(3)).^2 + X(2); 
    XiQuadr=resnorm;
    dgrfreedom=fitlength-3; % 3 fitparameters
    goodness=XiQuadr/dgrfreedom;
    errorX=reshape(full(sqrt(goodness*diag(inv(jaco'*jaco))))',size(X) );
    output= [X(1)/4 errorX(1)/4];
elseif plotfit==4
    fitlength=4; %6
    x0=[0.1,0.005]; % x(1) = 4D; x(2) = offset;
    [X2,resnorm2,residual2,eflag2,out2,lambda2,jaco2] = lsqnonlin(@FitMeanMSDlinearGJB18,x0,[-40,0],[40,0.0005],optimset('display','off'), SigmaMSD, Tmsd, AvarMSD, fitlength);
    XiQuadr2=resnorm2; %Xi=sum[(fit-obs)/std]^2 where std is integrated into the fitfunction, as the weight factor.
    dgrfreedom2=fitlength-2; % 2 fitparameters
    goodness2=XiQuadr2/dgrfreedom2;
    errorX2=reshape(full(sqrt(goodness2*diag(inv(jaco2'*jaco2))))',size(X2) );
    output= [X2(1)/4 errorX2(1)/4];
    
    fitlength=plotlength;
    x0=[0.1,1,0.02]; % x(1) = 4BD; x(2) = A; x(3) = <Rc^2>.
    [X,resnorm,residual,eflag,out,lambda,jaco] = lsqnonlin(@FitMeanMSDconfinedGJB18,x0,[-40,0,0.001],[40,2,100],optimset('display','off'), SigmaMSD, Tmsd, AvarMSD,X2(2), fitlength);
    TmsdFit=[0,Tmsd];   % add a point at 0 for the fit line
    FitMSD=X(3).*(1-X(2).*exp(-(X(1)/X(3)).*TmsdFit) ) + X2(2);
    XiQuadr=resnorm;
    dgrfreedom=fitlength-3; % 3 fitparameters
    goodness=XiQuadr/dgrfreedom;
    errorX=reshape(full(sqrt(goodness*diag(inv(jaco'*jaco))))',size(X) ); % [1,3]);
end;

%plot mean MSD with error and fit
figure
errorbar(tlag,MeanMSD,ErrorMSD,'o'); 
xlabel('Timelag (sec)');
if interactive5==0 % no truncation, all trajectories fully analysed.
     title('Average MSD, all trajectories fully analyzed');
 elseif interactive5==1 % truncate trajectories to 'truncate', preferably to min. tracelength
     title(['Average MSD, all trajectories truncated to '  int2str(truncate) ' .']);
 elseif interactive5==2 % only take trajectories longer than 'truncate'. 
     title(['Average MSD, trajectories longer than '  int2str(truncate) ' .']);
 elseif interactive5==3 % only take trajectories shorter than 'truncate'.
     title(['Average MSD, trajectories shorter or equal than ' int2str(truncate) ' .']);
end;


hold on
if plotfit==1 %plot linear fit in figure.
    [AX,H1,H2] = plotyy(TmsdFit,FitMSD,tlag,SizeMSD);  
    set(H1,'LineStyle','-', 'color', 'r')
    set(H2,'LineStyle',':')
    text(0.1,1*max(AvarMSD), ['D = ',num2str(X(1)/4,'%11.2g'), ' \pm ',num2str(errorX(1)/4,'%11.1g'),' µm^2/s (first 4 points)'] ,'HorizontalAlignment','left','FontSize',12);   %
    text(0.1,0.9*max(AvarMSD), ['Offset = ',num2str(X(2),'%11.2g'), ' \pm ',num2str(errorX(2),'%11.1g'),' µm^2'] ,'HorizontalAlignment','left','FontSize',12);   %
    text(0.1,0.8*max(AvarMSD), ['Goodness fit = ',num2str(goodness)] ,'HorizontalAlignment','left','FontSize',12);
elseif plotfit==2 %plot t^alpha fit in figure.
    [AX,H1,H2] = plotyy(TmsdFit,FitMSD,tlag,SizeMSD); 
    set(H1,'LineStyle','-', 'color', 'r')
    set(H2,'LineStyle',':')
    text(0.1,1*max(AvarMSD), ['D = ',num2str(X(1)/4,'%11.2g'), ' \pm ',num2str(errorX(1)/4,'%11.1g'),' µm^2/s'] ,'HorizontalAlignment','left','FontSize',12);   %
    text(0.1,0.9*max(AvarMSD), ['Offset = ',num2str(X(2),'%11.2g'), ' \pm ',num2str(errorX(2),'%11.1g'),' µm^2'] ,'HorizontalAlignment','left','FontSize',12);   %
    text(0.1,0.8*max(AvarMSD), ['Alpha = ',num2str(X(3),'%11.2g'), ' \pm ',num2str(errorX(3),'%11.1g')] ,'HorizontalAlignment','left','FontSize',12);
    text(0.1,0.7*max(AvarMSD), ['Goodness fit = ',num2str(goodness)] ,'HorizontalAlignment','left','FontSize',12);
elseif plotfit==3 %plot directed motion fit in figure.
    [AX,H1,H2] = plotyy(TmsdFit,FitMSD,tlag,SizeMSD); 
    set(H1,'LineStyle','-', 'color', 'r')
    set(H2,'LineStyle',':')
    text(0.1,1*max(AvarMSD), ['D = ',num2str(X(1)/4,'%11.2g'), ' \pm ',num2str(errorX(1)/4,'%11.1g'),' µm^2/s'] ,'HorizontalAlignment','left','FontSize',12);   %
    text(0.1,0.9*max(AvarMSD), ['Offset = ',num2str(X(2),'%11.2g'), ' \pm ',num2str(errorX(2),'%11.1g'),' µm^2'] ,'HorizontalAlignment','left','FontSize',12);   %
    text(0.1,0.8*max(AvarMSD), ['V = ',num2str(X(3),'%11.2g'), ' \pm ',num2str(errorX(3),'%11.1g'),' µm/s.'] ,'HorizontalAlignment','left','FontSize',12);
    text(0.1,0.7*max(AvarMSD), ['Goodness fit = ',num2str(goodness)] ,'HorizontalAlignment','left','FontSize',12);
elseif plotfit==4 %plot corraled motion fit in figure.
    [AX,H1,H2] = plotyy(TmsdFit,FitMSD,tlag,SizeMSD); 
    set(H1,'LineStyle','-', 'color', 'r')
    set(H2,'LineStyle',':')
    text(0.1,1*max(AvarMSD), ['B*D = ',num2str(X(1)/4,'%11.2g'), ' \pm ',num2str(errorX(1)/4,'%11.1g'),' µm^2/s'] ,'HorizontalAlignment','left','FontSize',12);   %
    text(0.1,0.9*max(AvarMSD), ['A = ',num2str(X(2)),' \pm ',num2str(errorX(2),'%11.1g')] ,'HorizontalAlignment','left','FontSize',12);
    text(0.1,0.8*max(AvarMSD), ['Offset = ',num2str(X2(2),'%11.2g'), ' \pm ',num2str(errorX2(2),'%11.1g'),' µm^2'] ,'HorizontalAlignment','left','FontSize',12);   %
    text(0.1,0.7*max(AvarMSD), ['<r^2> = ',num2str(X(3),'%11.2g'), ' \pm ',num2str(errorX(3),'%11.1g'),' µm^2'] ,'HorizontalAlignment','left','FontSize',12);
    text(0.1,0.6*max(AvarMSD), ['D_{initial} = ',num2str(X2(1)/4,'%11.2g'), ' \pm ',num2str(errorX2(1)/4,'%11.1g'),' µm^2/s (first 4 points)'] ,'HorizontalAlignment','left','FontSize',12);
    text(0.1,0.5*max(AvarMSD), ['Goodness fit = ',num2str(goodness)] ,'HorizontalAlignment','left','FontSize',12);
elseif plotfit==0 %plot corraled motion fit in figure.
    [AX,H1,H2] = plotyy(TmsdFit,FitMSD,tlag,SizeMSD); 
    set(H1,'LineStyle','none')
    set(H2,'LineStyle',':')
end;

% set two times Y axis and two linestyles
set(get(AX(1),'Ylabel'),'String','MSD (µm^2)')
set(get(AX(2),'Ylabel'),'String','Number of steps')
hold off    

% put some data in file, to be able to present plot in Origin. 
tlag=[0,tlag];
MeanMSD=[0,MeanMSD];
ErrorMSD=[0,ErrorMSD];
DataFile=[tlag; MeanMSD; ErrorMSD; TmsdFit; FitMSD]';
save('datalastMSD_ASCII.txt','DataFile', '-ASCII' , '-tabs'); 
       
%x = lsqcurvefit(fun,x0,xdata,ydata)
%x = lsqnonlin(fun,x0,lb,ub,options,P1,P2,...) passes the problem-dependent 
%parameters P1, P2, etc. directly to the function fun. Pass an empty matrix 
%for options to use the default values for options.