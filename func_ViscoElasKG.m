function C_Visco=func_ViscoElasKG(K,G)
%Author: Andrew Matejunas
%Date Completed: 2021-09-27

%current version
%Change Log:

%This function is written to compute the stiffness tensor from Bulk
    %Modulus, K, and shear modulus, G, inputs
    
%Function input arguments
    %K-  Vector of bulk modulus in Pa
    %G- Vector of shear modulus in Pa\

%Function outoput arguments
    %C_elas(i,j,k,l,m)- elasticity tensors
        %(i,j,k,l,:)- 4th order tensor indexes
        %(:,:,:,:,m)- number of elements
        
%% convert bulk and shear modulus to Lame' constant
Lame=K-2/3*G;

%% define Kroeniker delta
Kdelta=eye(3); %3x3 identity matrix

%Calculate elasticity matrix
for i=1:3;
    for j=1:3;
        for k=1:3;
            for l=1:3;          
                for m=1:length(G)
                    C_elas(i,j,k,l,m)=Lame(m)*Kdelta(i,j)*Kdelta(k,l)...
                        +G(m)*(Kdelta(i,k)*Kdelta(j,l)+Kdelta(i,l)*Kdelta(j,k));
                end
            end
        end
    end
end
    
C_Visco=C_elas;
end