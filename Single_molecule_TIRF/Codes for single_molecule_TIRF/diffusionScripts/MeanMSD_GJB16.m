% Gebaseerd op Cumprob_GJB6   

function dtotal = MeanMSD_GJB16(TR,pixelSize,lag);

r=size(TR,1);
rr=size(TR,2); 
NrOfTrs=TR(r,rr); %this is not the number of trajectories, but the biggest trajectory number.

u=0;
for j=1:NrOfTrs % check all trajectory numbers, if they are in the dataset
    tracenr=j;
    trace=zeros(1,rr);
    
    s=0;
    for i=1:r
        if TR(i,rr)==tracenr % replace zero vector 'trace' if the trajectory number 'tracenr' is present in the TR dataset.
            trace(s+1,:)=TR(i,:);
            s=s+1;
        end;
    end; 

    % calc Mean Squares, part adapted from msddGJB16, with correction for missing frames in trajectories!
    if trace(1,:)~=zeros(1,rr) % only calculate if 'trace' is not just a zero vector.
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
        
        %Squared Displacements berekenen door kolommen molAA ten opzichte van molBB te verschuiven en van elkaar af te trekken.
        if lag<newtracelength
            DISP=zeros(1,2);
            molAA=molAA(lag+1:end,:);                          % AA later laten beginnen met timelag
            molBB=molBB(1:end-lag,:);	    					% BB eerder laten ophouden met timelag
            DISP=(molAA-molBB).^2;							    %squared displacement per coördinaat
            
            s=1;                 % remove the squared displacements which are unrealistically big, those are the ones based on a skipped frame.
            disp=zeros(1,2);
            for h=1:size(DISP,1)
                if DISP(h,1)<90000
                    disp(s,:)=DISP(h,:);
                    s=s+1;
                end;
            end;
            
            d = sum((disp)')'; % differencex^2 + differencey^2     (trace(lag+1:end,1:2)-trace(1:end-lag,1:2)).^2 
            
            % built column: 
            dlength=size(d);
            dtotal(u+1:u+dlength)=d; 
            u=u+dlength;
        end;
    end; % end only calculate when trace is not just a zero vector.
    
end; % end check all trajectory numbers, if they are in the dataset

dtotal=dtotal';
