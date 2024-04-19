% This script is written to combine two separate spatial smoothing
    % intervals into a single plot for a viscoelastic material
    % characterized with a single element Generalized Maxwell model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
close all; clear variables; clc

%% Choose Folder to save combined files in

%% Load G errors for smaller kernels
[GLSname,GLSpath]=uigetfile('#.')
