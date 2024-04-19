%% Initialize
clear all, close all, clc

%% Choose Grid Method Data File
[GMfile,GMpath]=uigetfile('*.mat',...
    'Choose data file containing experimental grid method data');

%% Load Grid Method Data File
fprintf('Loading Grid Method Data \n')
load(strcat(GMpath,'/',GMfile));
fprintf('Loading Complete \n')

%% Define Test Designation
TestDeg=char(cell2mat(inputdlg('Test Designation')));

%% Evaluate Stress Gage Equations and Generate PLots
SG = func_IBIIExperimentalSG(accel,strain,pos,time,TestDeg,...
    material.rho);

%% Save
fprintf('Saving Data \n')
savename=strcat(TestDeg,'_SGdata.mat');
save(savename)
fprintf('Done \n')
