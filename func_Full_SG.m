function SGstresses=func_Full_SG(accel,X_vec,time,density)
% This function is written to compute stress gauge stresses using the
    % stress gauge equations

%Author: Andrew Matejunas

%completed: 2022-05-31

%Change log
    %2022-09-22- Added output of surface averaged acceleratations

%Function input arguments (note; Unless otherwise noted all units are in
%their associated standard SI units i.e (kg,m,s, Pa...)
    %accel- structure of acceleration arrays with fields
        %x- accelerations in the x direction
        %y- in the y direction
            %array dimensions [# of Ypoints, # X points, # time increments]
    %X_vec- Vector of X coordinates
    %time- Structure of time data
    %density- material density in kg/m^3

%Function output arguments
    %SGstresses- structure of calculated stress-gage stresses with
        %dimensions [# X points, time steps] and fields
            %x- normal stresses in X direction
            %s- in-plane shear stresses
            %Ay_surfx- Surfaced average y accelerations
            %Ax_surfx- Surface averaged x accelerations
     

%% convert into forms used by the code
time_vec=time.vec;
Ay=accel.y;
Ax=accel.x;

%% compute average accelerations
Ay_avx=mean(Ay);
Ax_avx=mean(Ax);

%% compute surface averaged acceleration in the x direction from top edge
    %to y0
    %allocate memory
Ay_surfx=zeros(length(X_vec),length(time_vec));
Ax_surfx=zeros(length(X_vec),length(time_vec));
SGstresses.s=zeros(length(X_vec),length(time_vec));
SGstresses.x=zeros(length(X_vec),length(time_vec));

for m=1:length(X_vec)
    for n=1:length(time_vec)
        Ay_surfx(m,n)=mean(Ay_avx(1,1:m,n));
        Ax_surfx(m,n)=mean(Ax_avx(1,1:m,n));
%% compute stress gage
SGstresses.s(m,n)=density*X_vec(m)*Ay_surfx(m,n);
SGstresses.x(m,n)=density*X_vec(m)*Ax_surfx(m,n);
SGstresses.AY_surfx=Ay_surfx;
SGstresses.Ax_surfx=Ax_surfx;
        
    end
end


end