function TRtotal_all=readTrajectories(timeLapse,varargin)
    %% trajectories=readTrajectories(0.098); % Read all trajectories, TimeLapse (time between two frames)=0.098s
    %% trajectories=readTrajectories(0.098,[4,5,28]); % Read all trajectories except 4, 5 and 28, TimeLapse (time between two frames)=0.098s
    
load(fullfile('results','TrackingPackage','tracks','Channel_1_tracking_result.mat'));
ww=1;

for ll=1:size(tracksFinal,1)
    inter_data=tracksFinal(ll).tracksCoordAmpCG;
    if size(inter_data,1)==1
        rrww=isnan(inter_data(:));
        if sum(rrww)>0
            rrww_inter=find(rrww);
            inter_data=inter_data(1:(rrww_inter(1)-1));% if there is a blinking, take the trajectory till the first blinking event
        end
        ii=1;
        for mm=1:8:max(size(inter_data))
            tracks(ii,:)=inter_data(mm:mm+7);
            ii=ii+1;
        end
        frame_number=tracksFinal(ll).seqOfEvents;
        frame_number_def=frame_number(1,1):(((numel(inter_data)/8)+frame_number(1,1))-1);%frame_number(2,1);
        TRtotal1=tracks;
        TRtotal1(:,3)=10000;
        TRtotal1(:,5)=10000;
        TRtotal1(:,4)=1.5;
        TRtotal1(:,6)=1.5;
        TRtotal1(:,7)=frame_number_def;
        TRtotal1(:,8)=ll; %ll is the "spot number"
        TRtotal1(:,9)=timeLapse;
        TRtotal1(:,10)=1;
        si_TRtotal1=size(TRtotal1,1);
        TRtotal_all(ww:ww+si_TRtotal1-1,:)=TRtotal1;
        ww=ww+si_TRtotal1;
        clear tracks TRtotal1
    end
end 

if nargin==2
    b=varargin{1};
    d=size(b);
    for t=1:max(d)
        z=b(t);
        [x] = find (TRtotal_all(:,8)==z);
        TRtotal_all(x,8)=0;
    end
    [x] = find (TRtotal_all(:,8)>0);
    TRtotal_all=TRtotal_all(x,:);
end

   
   
   
   
   