function TRlength = LengthAdjustAndSelect12(TR, selectlength, truncate, interactive5)
% function to truncate trajectories to make them all the same length, so all
% trajectories will have the same weight in the mean MSD curve. If longer
% trajectories have different diffusion behaviour than shorter ones, the
% longer won't get dominant, which would normally be the case. One can also
% choose to select the part of the trajectories which is bigger than or
% smaller than 'truncate'. If one wants to see if shorter trajectories have
% a different contribution than bigger ones, one can let the short ones out
% by selecting all trajectories which are bigger as or equal to selectlength. 

% Input: 
% > traced points, TR 
% > selectlength, select the traces which are 'selectlength' or longer, out of the dataset. 
% if 'selectlength' is equal to or smaller than tracelength, all trajectories will be selected.
% > truncate, to make all traces the same length
% > interactive5, 0=all selected trajectories, 1=all but truncated, 2=only longer than, 
% 3=only shorter or equal to truncate

% Output: TRlength, the trajectories, selected and/or adjusted (truncated)

r=size(TR,1);
rr=size(TR,2); 
NrOfTrs=TR(r,rr); %this is not the number of trajectories, but the biggest trajectory number.

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
     TRlength=TR;
 elseif interactive5==1 % truncate trajectories to 'truncate', preferably to min. tracelength
     TRlength=TRtrunc;
 elseif interactive5==2 % only take trajectories longer than 'truncate'. 
     TRlength=TRlong;
 elseif interactive5==3 % only take trajectories shorter than 'truncate'.
     TRlength=TRshort;
 end;