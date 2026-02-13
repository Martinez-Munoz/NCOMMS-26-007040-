



%% Esta funcion es para separar las trayectorias por su longitud. Solo tienes que decidir el numero de puntos por el que quieres separarlas. 

% Este valor "length_threshold" lo tendras que poner tu. La funcion te
% genera dos TRtotal. shortTrajectories son todas las trayectorias
% con una duracion menor del valor que tu le has puesto.
% longTrajectories son todas las trayectorias con una duracion mayor
% del valor que tu has elegido.


function [shortTrajectories, longTrajectories] = separateTrajectoriesByLength(TRtotal,length_threshold)

whole_trajectories=unique(TRtotal(:,8));


shortTrajectories=[];


longTrajectories=[];

for ll=1:numel(whole_trajectories)
    
    xx=find(TRtotal(:,8)==whole_trajectories(ll));
    
    TRtotal_inter=TRtotal(xx,:);
    
    if numel(xx)>(length_threshold-1)
        
        longTrajectories=[longTrajectories; TRtotal_inter];
        
    else
        
        shortTrajectories=[shortTrajectories;TRtotal_inter];
        
    end
    
end


