
function [S_MSS_def,Scaling_parameters_def,MSD_def]=momentScalingSpectrum(TRtotal1,pixelsize,tau)


wzxr=1;

Scaling_parameters_def=[];

MSD_def=[];

uni_TRtotal1=unique(TRtotal1(:,8));


for wwz=1:numel(unique(TRtotal1(:,8)))
    
    xcv=find(TRtotal1(:,8)==uni_TRtotal1(wwz));
    
    TRtotal=TRtotal1(xcv,:);
    
jjw=1;

for v=0:6 % The moments
    
for ll=1:round(size(TRtotal,1)/4) %Take only into account 25% of the trajectory
    
    disp_x=[];
    disp_y=[];
    
    for mm=1:(size(TRtotal,1)-ll)
        
        disp_x_inter=TRtotal(mm+ll,1)-TRtotal(mm,1);
        disp_x=[disp_x,disp_x_inter];
        disp_y_inter=TRtotal(mm+ll,2)-TRtotal(mm,2);
        disp_y=[disp_y,disp_y_inter];
        
    end
    
    disp_x=disp_x.*pixelsize;
    
    disp_y=disp_y.*pixelsize;
    
    disp_z=disp_x.^2+disp_y.^2; 
    
    disp_z=sqrt(disp_z);
    
    disp_z=disp_z.^v;%% The moments v
    
    MSD(ll,jjw+1)=mean(disp_z);
    MSD(ll,jjw)=ll*tau;
    MSD(ll,jjw+2)=std(disp_z)/sqrt(numel(disp_z)); %SEM of the MSD
    
end

jjw=jjw+3;


end

kks=1;

for wxz=1:3:size(MSD,2)
    
    MSD_inter=MSD(1:10,wxz:wxz+1); %10 points
    
    MSD_to_fit=log10(MSD_inter); %logarithm
    
    cfun=fittype('rr*dd+ee','coefficients',{'rr','ee'},'independent','dd'); 
     
    %%% JR %%%options=fitoptions('Method','LinearLeastSquares');
    options=fitoptions('Method','NonlinearLeastSquares','StartPoint',[0 0]);
    
    [cfun,dfe,Jacobian]=fit(MSD_to_fit(:,1),MSD_to_fit(:,2),cfun,options);
   
    slope=cfun.rr;offset=cfun.ee; %%Extracting fittingsparameter
    
    nor_resi=norm(Jacobian.residuals); goodness=nor_resi/dfe.dfe;errorX=full(sqrt(goodness*diag(inv(Jacobian.Jacobian'*Jacobian.Jacobian)))); %%Calculating errors as Gert-Jan.
   
    Dif_coeff=.25*exp(offset); 
    
    Scaling_parameters(kks,1)=slope; %Scaling coefficient 
    Scaling_parameters(kks,2)=offset; % offset, related to D 
    Scaling_parameters(kks,3)=errorX(1); %error in tilt
    Scaling_parameters(kks,4)=errorX(2); % error in Offset
    Scaling_parameters(kks,5)=Dif_coeff; %Diffusion coefficient in um2/s
    Scaling_parameters(kks,6)=errorX(2)*Dif_coeff; %The error of the diffusion coefficient: The error propagation of an exponential (a*exp(b*x)) is (b*error(x)*a*exp(b*x))). In this case, b=1.
    Scaling_parameters(kks,7)=uni_TRtotal1(wwz); %Trajectory number
    Scaling_parameters(kks,8)=kks-1;%The moment
    kks=kks+1;
end


xx=0:6;

xx=xx';

yy=Scaling_parameters(:,1);

cfun=fittype('rr*dd+ee','coefficients',{'rr','ee'},'independent','dd'); 
     
%%% JR %%% options=fitoptions('Method','LinearLeastSquares');
options=fitoptions('Method','NonlinearLeastSquares','StartPoint',[0,0]);
    
[cfun,dfe,Jacobian]=fit(xx,yy,cfun,options);
   
slope_S_MSS=cfun.rr;offset_S_MSS=cfun.ee; %%Extracting fittingsparameter
    
nor_resi=norm(Jacobian.residuals); goodness=nor_resi/dfe.dfe;errorX_S_MSS=full(sqrt(goodness*diag(inv(Jacobian.Jacobian'*Jacobian.Jacobian)))); %%Calculating errors as Gert-Jan.
   
S_MSS_def(wzxr,1)=slope_S_MSS; %slope of the Moment Scaling Spectrum
S_MSS_def(wzxr,2)=offset_S_MSS; %Offset of the MSS
S_MSS_def(wzxr,3)=errorX_S_MSS(1); %Error of the MSS
S_MSS_def(wzxr,4)=errorX_S_MSS(2); %Error of the offset
S_MSS_def(wzxr,5)=uni_TRtotal1(wwz); %trajectory number



Scaling_parameters_def=[Scaling_parameters_def;Scaling_parameters];

MSD(:,22)=uni_TRtotal1(wwz);%trajectory number

MSD_def=[MSD_def;MSD];



clear Scaling_parameters MSD

wzxr=wzxr+1;



end

