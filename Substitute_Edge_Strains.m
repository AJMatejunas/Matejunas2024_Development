%This code is intended to investigate the effects of errors in the grid
%method displacement results near the edges of the specimens by replacing
%them with interpolated FE data.

%% Initialize 
clear all; close all; clc

%% Select FE and GM data Files
[FEfile,FEpath]=uigetfile('*.mat','Choose file containing pure FE data');
[GMfile,GMpath]=uigetfile('.mat','Choose GM data file');

%% load FE Data
fprintf('Loading Finite Element Data \n')
FEdata=load(strcat(FEpath,'/',FEfile));
FE.strain=FEdata.strain;
FE.disp=FEdata.disp;
FE.accel=FEdata.accel;
FE.X_vec=FEdata.X_vec;
FE.Y_vec=FEdata.Y_vec;

%% Free up RAM by deleting unnecessary variables
% clear FEdata

%% load GM Data
fprintf('Loading GM data \n')
GM=load(strcat(GMpath,'/',GMfile));
%NOTE: This was changed for effeciency of using the output of this 
    %substitution with my existing codes
    
% GM.strain=GM.data.strain;
% GM.disp=GM.data.disp;
% GM.accel=GM.data.accel;
 GM.X_vec=GM.pos.x;
 GM.Y_vec=GM.pos.y;
% GM.grid=GM.data.grid;
% GM.time=GM.data.time;

testdeg=GM.testdeg;

%% Substitutionn options
fprintf('Defining substitution options')

%define parameters to substitute
quest='Which Parameters to perform substitutions on';
SubParam.fields=questdlg(quest,'Substitution Parameters',...
    'Strain Only','Strain and Accel','Disp, Strain, Accel',...
    'Disp, Strain, Accel')
SubParam.num=str2num(cell2mat(inputdlg('How Many Pixels to Substitute')));

switch SubParam.fields
    case 'Strain Only'
        SubParam.StrainSub=true;
        SubParam.DispSub=false;
        SubParam.AccelSub=false;
    case 'Strain and Accel' %only differentiated parameters are replaced
        SubParam.StrainSub=true;
        SubParam.DispSub=false;
        SubParam.AccelSub=true;
    case 'Disp, Strain, Accel'
        SubParam.StrainSub=true;
        SubParam.DispSub=true;
        SubParam.AccelSub=true;
end
        

%% Perform Substitution
fprintf('performing substitution \n')
[SubStrain,SubDisp,SubAccel,SubRecord]=func_SubEdgeFieldsV1(SubParam,...
    FE,GM);

%% replace raw data with substituted data
GM.strain=SubStrain;
GM.disp=SubDisp;
GM.accel=SubAccel;
GM.SubRecord=SubRecord;

%% Save data
savename=strcat(testdeg,'_GM_subData.mat');
save(savename,'-Struct','GM')