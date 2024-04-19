% This code is written to extract the inputs for the cost function
    %evalation code
    
%% load data file
name=uigetfile('.mat');
load(name);

%% define test designation
prompt='input test designation';
testdeg=char(inputdlg(prompt));

%% save appropriate workspace variables
save(strcat(testdeg,'_SGdata.mat'),'SG','time','X_vec','strain','accel');