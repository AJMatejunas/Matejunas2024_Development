function [SG,accel,strain,X_vec]=func_KinFieldsDownsampling(SG,accel,...
    strain,X_vec,time,Xds,Yds)

% This script is written to downsample data for VFM minimizations
% Author: Andrew Matejunas

%Version History/Changelog
    %2022-06-02- Converted original script to a function for added
        %flexibility, effiecency, and ease of integration with existing
        %codes. Also converted notation to more closely match other codes.
    %2022-07-14: Fixed Y downsampling


%% Choose downsampling parameters
if nargin<6 %allows downsampling factor to be defined if it is not already
            % defined in fuction inputs
            prompt='Downsampling factor';
            title='define downsample options';
            default_factor={'5'};
            
            % define whether downsampling occurs in x coordinates, y coordinates, or
            % both
            
            quest='Downsample in X, Y, or both';
            
            DownOpts=questdlg(quest,'downsampling options',...
                'X only', 'Y_only', 'X and Y', 'X and Y')
            
            %downsampling factor
            SampFact=str2num(cell2mat(inputdlg(prompt,title,[1],default_factor)));
            
            clear prompt title
            
            % convert downsampling factor into separate factors for X and Y
                switch DownOpts
                    case 'X only'
                        Xds=SampFact;
                        Yds=1;
                    case 'Y only'
                        Yds=SampFact;
                        Xds=1;
                    case 'X and Y'
                        Xds=SampFact;
                        Yds=SampFact;
                end
       
end



 X_vec=downsample(X_vec,Xds);
        
        
         for n=1:length(time.vec)
            %Note downsampling occurs on each column of a vector so to
            %downsample the x coordinates the matrices must be transposed
            %prior to downsampling then transposed back. 
             
            %% accelerations downsampled along X
            accel.xDSX(:,:,n)=downsample(squeeze(accel.x(:,:,n))',Xds)';
            accel.yDSX(:,:,n)=downsample(squeeze(accel.y(:,:,n))',Xds)';
            
            %% strains downsampled along X
            strain.xDSX(:,:,n)=downsample(squeeze(strain.x(:,:,n))',Xds)';
            strain.yDSX(:,:,n)=downsample(squeeze(strain.y(:,:,n))',Xds)';
            strain.sDSX(:,:,n)=downsample(squeeze(strain.s(:,:,n))',Xds)';
                                  
                       
        end
   

         SG.x=downsample(SG.x,Xds);
         SG.s=downsample(SG.s,Xds);
            for n=1:length(time.vec)
                %% accelerations downsampled along Y
                accel.xDSXY(:,:,n)=downsample(squeeze(accel.xDSX(:,:,...
                    n)),Yds);
                accel.yDSXY(:,:,n)=downsample(squeeze(accel.yDSX(:,:,...
                    n)),Yds);
            
                %% strains downsampled along Y
                strain.xDSXY(:,:,n)=downsample(squeeze(strain.xDSX(:,:,...
                    n)),Yds);
                strain.yDSXY(:,:,n)=downsample(squeeze(strain.yDSX(:,:,...
                    n)),Yds);
                strain.sDSXY(:,:,n)=downsample(squeeze(strain.sDSX(:,:,...
                    n)),Yds);
            end
            
         accel.x=accel.xDSXY;
         accel.y=accel.yDSXY;
         strain.x=strain.xDSXY;
         strain.y=strain.yDSXY;
         strain.s=strain.sDSXY;
         
accel=rmfield(accel,{'xDSX','yDSX','xDSXY','yDSXY'});
strain=rmfield(strain,{'xDSX','yDSX','sDSX','xDSXY','yDSXY','sDSXY'});
end