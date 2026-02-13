function [TRtotal_simulated]=simulateDirectedTrajectories(Din,pixelsize,delta_t_in,num_traj,length_trajectory)

D=Din/pixelsize^2; %% Din has to be in um2/s

delta_t=delta_t_in/100; %delta_t has to be in seconds.
                        %The time delta_t has to be 100 times less because the software generates
                        %10000 points and choose every 100 frames to end up
                        %having 100 points
                        
Velocity=1; %it has to be in um/s

Velocity_in_pixel=Velocity/pixelsize; %Velocity in pixel

Velocity_in_pixel_x_component=sqrt(Velocity_in_pixel^2/2);% V^2=Vx^2+Vy^2 and Vx=Vy

axis_contribution=Velocity_in_pixel_x_component*delta_t; %Displacement in 1D which will be added to the displacement coming from the brownian motion









%num_traj has to be the power of a given number. f.e. num_traj=4 since
%2^2=4.

% The localization accuracy (local_a) has to be also in pixel.

%noise=(local_a)^2;

r=1;v=1;w=1;

x_y_1=round(random('unif',100,256,round((sqrt(num_traj)))));%generate all x/y starting positions

x_y_1=x_y_1(:);

x_y_2=round(random('unif',100,256,round((sqrt(num_traj))))); x_y_2=x_y_2(:);

si_x_y2=size(x_y_2); si_x_y1=size(x_y_1);

x_y=(x_y_1); x_y((si_x_y1(1)+1):(2*si_x_y1(1)))=x_y_2;


for z=1:num_traj


Variance=2*D*delta_t; % It is "2D" because we generate 1 trajectory in x/y.


Std_d=sqrt(Variance);

%var_x=randmat(0,Variance,1,9999); var_y=randmat(0,Variance,1,9999);

var_x=random('norm',0,Std_d,500); var_y=random('norm',0,Std_d,500);

var_x=var_x(:); var_y=var_y(:); 

% var_x=var_x+axis_contribution;var_y=var_y+axis_contribution; %Making the simulation directed transport. If you want brownian movement, simply comment this line.

x=x_y(w,1); 

y=x_y((w+1),1);


for i=1:(length_trajectory*100-1)
    j=i+1;
    x(j,1)=x(i,1)+var_x(i,1);
    y(j,1)=y(i,1)+var_y(i,1);
    
end


x_new=x(1:100:length_trajectory*100);
 
y_new=y(1:100:length_trajectory*100);


TRtotal(r:(r+(length_trajectory-1)),1)=x_new;TRtotal(r:(r+(length_trajectory-1)),2)=y_new;

TRtotal(r:(r+(length_trajectory-1)),3)=20000;TRtotal(r:(r+(length_trajectory-1)),4)=1.36;TRtotal(r:(r+(length_trajectory-1)),5)=20000;

TRtotal(r:(r+(length_trajectory-1)),6)=1.36; TRtotal(r:(r+(length_trajectory-1)),7)=1:length_trajectory; TRtotal(r:(r+(length_trajectory-1)),8)=z;TRtotal(r:(r+(length_trajectory-1)),9)=delta_t_in; TRtotal(r:(r+(length_trajectory-1)),10)=1;


si_TR=size(TRtotal); 

r=si_TR(1)+1;

w=z+2;

clear x_new y_new





end

TRtotal_simulated=TRtotal;


















