function Shear_SG=func_ShearSG(accel,X_vec,time,density)
% This function is written to compute shear stress using the shear stress
% equation

%% convert into forms used by the code
time_vec=time.vec;
Ay=accel.y;

%% compute average accelerations
Ay_avx=mean(Ay);

%% compute surface averaged acceleration in the x direction from top edge
    %to y0
    %allocate memory
Ay_surfx=zeros(length(X_vec),length(time_vec));
Shear_SG=zeros(length(X_vec),length(time_vec));
for m=1:length(X_vec)
    for n=1:length(time_vec)
        Ay_surfx(m,n)=mean(Ay_avx(1,1:m,n));
%% compute shear stress gage
Shear_SG(m,n)=density*X_vec(m)*Ay_surfx(m,n);
        
    end
end


end