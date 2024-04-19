%This code is written to censor data (remove it from consideration) in an
    %IBII test. Use it to remove data points that may be the sources of
    %error in parameter identification. Such as if there is a shear stress
    %concentration on a half height impact.
    
    
%Change log:
    %2022-02-09: Fixed bug where testdeg was improperly saved
    %          : Removed option to censor data in Y. This breaks the stress
    %          gage equation. 
 %% initialize
 clear all, close all, clc
 
%% Choose raw data to load
SGname=uigetfile('choose Stress Gauge data to Censor');
Censdesig=inputdlg('enter test designation');

%% Load raw data
load(SGname);
testdeg=char(Censdesig);
clear desig Censdesig

%% choose censoring parameters
% count the number of data points for an educated guess on how many data
% points to remove
%% 
Xcount=num2str(length(X_vec));
Ycount=num2str(length(Y_vec));

Xquest=strcat('Censoring options in X ',Xcount,' datapoints');
Xopt=questdlg(Xquest,'X censoring options',...
    'Impact only','Free only','Impact and Free','Impact and Free');
clear Xquest

dlgtitle='X censoring options';
Ydlgtitle='Y censoring options';
switch Xopt
    case 'Impact only'
    prompt='How many datapoints to remeove from the impact edge?';
    default_factor={'10'}; %how many datapoints to remove from the edge
    ImpCens=str2num(cell2mat(inputdlg(prompt,...
        dlgtitle,...
        1,... %size of the input window
        default_factor)));
    
    case 'Free only'
        prompt='How many datapoints to remove from the free edge?';
        default_factor='10';
        FreeCens=str2num(cell2mat(inputdlg(prompt,...
            dlgtitle,...
            1,... %size of the input window
            default_factor)));
            
    case 'Impact and Free'
       prompt={'how many data points to remove from the free edge?',...
               'how many data points to remove form the impact edge?'};
       default_factor={'10','10'};
       CensOpt=str2num(cell2mat(inputdlg(prompt,...
         dlgtitle,...
         [1,35],... %size of the input window
         default_factor)));
     FreeCens=CensOpt(1);
     ImpCens=CensOpt(2);
end
clear prompt

%Y censoring option is removed. It may be worth looking into when using
    %parametric virtual fields
Yopt='No';
%Yquest=strcat('Censor top and bottom edges? ',Ycount,' datapoints in Y');

%Yopt=questdlg(Yquest,'Y censoring options',...
    %'Yes','No','Yes');
clear Yquest

switch Yopt
    case 'Yes'
    prompt='How many datapoints to remeove from the top and bottom';
    default_factor={'10'};
    YCens=str2num(cell2mat(inputdlg(prompt,...
            Ydlgtitle,...
            [1,35],... %size of the input window
            default_factor)));
    clear prompt
end

if exist('accel','var')==0
    accel.x=zeros(size(strain));
    accel.y=zeros(size(strain));
end

if exist('avg_Ystrain','var')==0
    avgY_strain=zeros(size(avgX_strain));
end

% 
% %% create Raw Data Structure
% %Kinematic fields
% Raw.strain=strain;
% 
% Raw.SG=SG;
% Raw.Shear_SG=Shear_SG;
% Raw.accel=accel;
% %geometric information
% Raw.X_vec=X_vec;
% Raw.Y_vec=Y_vec;
% 
% %averages
% Raw.avg_accel=squeeze(avg_accel);
% Raw.avgX_strain=squeeze(avgX_strain);
% Raw.avgXY_strain=squeeze(avgXY_strain);
% Raw.avgY_strain=squeeze(avgY_strain);

avg_accel=squeeze(avg_accel);
avgX_strain=squeeze(avgX_strain);
avgXY_strain=squeeze(avgXY_strain);
avgY_strain=squeeze(avgY_strain);

%% Censor data in X
switch Xopt
    case 'Impact only'
%X_vector
X_vec((end-ImpCens):end)=[];

%censor strain
strain.x(:,(end-ImpCens):end,:)=[];
strain.y(:,(end-ImpCens):end,:)=[];
strain.s(:,(end-ImpCens):end,:)=[];

%censor accel
accel.x(:,(end-ImpCens):end,:)=[];
accel.y(:,(end-ImpCens):end,:)=[];

%censor SG
SG((end-ImpCens):end,:)=[];
Shear_SG((end-ImpCens):end,:)=[];

%averages
avg_accel((end-ImpCens):end,:)=[];
avgX_strain((end-ImpCens):end,:)=[];
avgXY_strain((end-ImpCens):end,:)=[];
avgY_strain((end-ImpCens):end,:)=[];

    case 'Free only'
%X_vector
X_vec(1:FreeCens)=[];

%censor strain
strain.x(:,1:FreeCens,:)=[];
strain.y(:,1:FreeCens,:)=[];
strain.s(:,1:FreeCens,:)=[];

%censor accel
accel.x(:,1:FreeCens,:)=[];
accel.y(:,1:FreeCens,:)=[];

%censor SG
SG(1:FreeCens,:)=[];
Shear_SG(1:FreeCens,:)=[];

%averages
avg_accel(1:FreeCens,:)=[];
avgX_strain(1:FreeCens,:)=[];
avgXY_strain(1:FreeCens,:)=[];
avgY_strain(1:FreeCens,:)=[];    

    case 'Impact and Free'
%% Impact end
%X_vector
X_vec((end-ImpCens):end)=[];

%censor strain
strain.x(:,(end-ImpCens):end,:)=[];
strain.y(:,(end-ImpCens):end,:)=[];
strain.s(:,(end-ImpCens):end,:)=[];

%censor accel
accel.x(:,(end-ImpCens):end,:)=[];
accel.y(:,(end-ImpCens):end,:)=[];

%censor SG
SG((end-ImpCens):end,:)=[];
Shear_SG((end-ImpCens):end,:)=[];

%averages
avg_accel((end-ImpCens):end,:)=[];
avgX_strain((end-ImpCens):end,:)=[];
avgXY_strain((end-ImpCens):end,:)=[];
avgY_strain((end-ImpCens):end,:)=[];

%% Free end
%X_vector
X_vec(1:FreeCens)=[];

%censor strain
strain.x(:,1:FreeCens,:)=[];
strain.y(:,1:FreeCens,:)=[];
strain.s(:,1:FreeCens,:)=[];

%censor accel
accel.x(:,1:FreeCens,:)=[];
accel.y(:,1:FreeCens,:)=[];

%censor SG
SG(1:FreeCens,:)=[];
Shear_SG(1:FreeCens,:)=[];

%averages
avg_accel(1:FreeCens,:)=[];
avgX_strain(1:FreeCens,:)=[];
avgXY_strain(1:FreeCens,:)=[];
avgY_strain(1:FreeCens,:)=[];    


end
%% Censor data in Y
switch Yopt
    case 'Yes' 
end
%% Save Censored Data
fprintf('Saving Censored Data \n')
save(strcat(testdeg,'_CensData.mat'));

