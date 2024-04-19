%This script is written to compare displacement, strain and acceleration
    %fields between pure finite element data, Grid method data with the
    %edges replaced, and corrected grid method data
    
    
 %Author: Andrew Matejunas
 
 %Date Completed: 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %% initialize
 clear all, close all, clc
 
 %% Choose Data files
 
 %Pure FE file
 [FEfile,FEpath]=uigetfile('*.dat',...
     'Choose pure finite element data file');
 
 %Corrected GM file
 [GMcorrfile,GMcorrpath]=uigetfile('*.dat',...
     'Choose GM data with corrected edges');
 
 %GM file with edge replaced
 [GMedgefile,GMedgepath]=uigetfile('*.dat',...
     'Choose GM data with finite element edges');
 
 %%  Load Data files
 fprintf('loading Pure Finite Element file \n')
 FE=load(strcat(FEpath,'/',FEfile));
 FE.path=FEpath;
 FE.file=FEfile;
 
 fprintf('loading grid method file with corrected edges \n')
 GMcorr=load(strcat(GMcorrpath,'/',GMcorrfile));
 GMcorr.file=GMcorrfile;
 GMcorr.path=GMcorrpath;
 GMcorr.X_vec=GMcorr.pos.x;
 GMcorr.Y_vec=GMcorr.pos.y;
 
 fprintf('Loading GM file with FE data substitutions along the edges \n')
 GMedge=load(strcat(GMedgepath,'/',GMedgefile));
 GMedge.file=GMedgefile;
 GMedge.path=GMedgepath;
 GMedge.X_vec=GMedge.pos.x;
 GMedge.Y_vec=GMedge.pos.y;
 
 fprintf('data loaded')
 
 %% Set up titles for the kinematic fileds plots
 
 Source1='Finite Element';
 Source2='GM (edge corrections)';
 Source3='GM with FE edge substitutions';
 
 %% Generate Plots
 answer=func_Compare3kinFields(FE,GMcorr,GMedge,Source1,Source2,Source3);
    
 