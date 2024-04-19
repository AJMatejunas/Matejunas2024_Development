% This script is written to downsample data for VFM minimizations
%% initialize
clear all
%% Load SG data and define test designation
SGname=uigetfile('Choose Stress gauge data to downsample');
DSdesig=inputdlg('enter test designation');

load(SGname)
testdeg=char(DSdesig);
clear desig DSdesig

%% Choose downsampling parameters
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
   %% Record raw values
        Raw.accel=accel;
        Raw.strain=strain;
        Raw.SG=SG;
        Raw.Shear_SG=Shear_SG;
        Raw.X=X_vec;
        %Raw.Y=Y_vec;
       
%% Downsample the data
        
   X_vec=downsample(Raw.X,SampFact)';
        
        
   clear accel strain SG Shear_SG
        for n=1:length(time.vec)
            
            %% accelerations downsampled along X
            accel.xDSX(:,:,n)=downsample(squeeze(Raw.accel.x(:,:,n))',SampFact)';
            accel.yDSX(:,:,n)=downsample(squeeze(Raw.accel.y(:,:,n))',SampFact)';
            
            %% strains downsampled along Y
            strain.xDSX(:,:,n)=downsample(squeeze(Raw.strain.x(:,:,n))',SampFact)';
            strain.yDSX(:,:,n)=downsample(squeeze(Raw.strain.y(:,:,n))',SampFact)';
            strain.sDSX(:,:,n)=downsample(squeeze(Raw.strain.s(:,:,n))',SampFact)';
                                  
                       
        end
    %%    
    switch DownOpts
        case 'X only'
            accel.x=accel.xDSX;
            accel.y=accel.yDSX;
            strain.x=strain.xDSX;
            strain.y=strain.yDSX;
            strain.s=strain.sDSX;
        
        SG=downsample(Raw.SG,SampFact);
        Shear_SG=downsample(Raw.Shear_SG,SampFact);
        case 'Y only'
           
           for n=1:length(time.vec)
            %% accelerations downsampled along X
            accel.x(:,:,n)=downsample(squeeze(Raw.accel.x(:,:,n)),SampFact);
            accel.y(:,:,n)=downsample(squeeze(Raw.accel.y(:,:,n)),SampFact);
            
            %% strains downsampled along Y
            strain.x(:,:,n)=downsample(squeeze(Raw.strain.x(:,:,n)),SampFact);
            strain.y(:,:,n)=downsample(squeeze(Raw.strain.y(:,:,n)),SampFact);
            strain.s(:,:,n)=downsample(squeeze(Raw.strain.s(:,:,n)),SampFact);
           end
           
        case 'X and Y'
         SG=downsample(Raw.SG,SampFact);
         Shear_SG=downsample(Raw.Shear_SG,SampFact);
            for n=1:length(time.vec)
                %% accelerations downsampled along X
                accel.x(:,:,n)=downsample(squeeze(accel.xDSX(:,:,n)),SampFact);
                accel.y(:,:,n)=downsample(squeeze(accel.yDSX(:,:,n)),SampFact);
            
                %% strains downsampled along Y
                strain.x(:,:,n)=downsample(squeeze(strain.xDSX(:,:,n)),SampFact);
                strain.y(:,:,n)=downsample(squeeze(strain.yDSX(:,:,n)),SampFact);
                strain.s(:,:,n)=downsample(squeeze(strain.sDSX(:,:,n)),SampFact);
           end               
    end    
    
%% save downsampled data
save(strcat(testdeg,'_DSdata.mat'));