
function [Classification_trajectories]=classifyLongTrajectories(TRtotal,pixelsize,Dboundary,outputSuffix)

unique_trajectories=unique(TRtotal(:,8));
rr=size(TRtotal,2);
tau=TRtotal(1,(rr-1)); % framerate, time between two frames

jj=1;

Nmax = numel(unique_trajectories);
for ll=1:Nmax

    disp(['Analyzing trajectory=' num2str(ll) '/' num2str(Nmax) ' ' datestr(now)])

    xx=find(TRtotal(:,8)==unique_trajectories(ll));
    
    TRtotal_trial=TRtotal(xx,:);
    
    [DiF] = AnalyseAllFilesGJB20_MAIN_MSS(TRtotal_trial,pixelsize,Dboundary);
    
    DiF1=DiF;
    
    Din=DiF1(1,1);
    
    length_trajectory=numel(xx);
    
    num_traj=10;
        
    [TRtotal_simulated]=simulateDirectedTrajectories(Din,pixelsize,tau,num_traj,length_trajectory);

    [S_MSS_def,Scaling_parameters_def,MSD_def]=momentScalingSpectrum(TRtotal_simulated,pixelsize,tau);
    
%     [mean,std]=normfit(S_MSS_def(:,1)); 
%     lower_threshold=mean-2*std; %   95% of the distribution
%     upper_threshold=mean+2*std;
    
    lower_threshold = quantile(S_MSS_def(:,1),0.025);
    upper_threshold = quantile(S_MSS_def(:,1),0.975);
    
    [S_MSS_def,Scaling_parameters_def,MSD_def]=momentScalingSpectrum(TRtotal_trial,pixelsize,tau);
    
    classified=2; %brownian. Default value.
    
    if S_MSS_def(1,1)<lower_threshold
        
        classified=1; %confined trajectory
        
    elseif S_MSS_def(1,1)>upper_threshold
        
        classified=3; %directed transport
        
    end
        
          
        
      Classification_trajectories(jj,1)=unique_trajectories(ll);   
  
      Classification_trajectories(jj,2)=classified;
    
      Classification_trajectories(jj,3)=S_MSS_def(1,1);
    
    
    
    jj=jj+1;
    
end

fnOut=fullfile('results','TrackingPackage','tracks',['trajectoryClassification' outputSuffix '.txt']);
fhOut = fopen(fnOut,'w');
for j=1:size(Classification_trajectories,1)
    fprintf(fhOut,'%d %d %f\n',Classification_trajectories(j,1),Classification_trajectories(j,2),Classification_trajectories(j,3));
end
fclose(fhOut);


    