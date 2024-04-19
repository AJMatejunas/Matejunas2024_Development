% This script is written to generate an input mat file for viscoelastic
% constitutive property verification from a matlab code

%Author: Andrew Matejunas
clear all, close all, clc

%% Define test designation
testdeg=char(inputdlg('input test designation','test designation'));

%% Define whether calculations will be E, or KG formulation

quest='E or K&G formulation';
formtype=questdlg(quest,'formulation type','E and nu','K and G','K and G');

%% define density
dlgtitle='density';
prompt='input density (kg/m^3)';
dims=[1,50];
MatProps.rho=str2double(cell2mat(inputdlg(prompt,dlgtitle,dims)));
    

%% Generate MatProps data structure
switch formtype
    case 'K and G'  
%% Set all E and nu to 0 
MatProps.E0=0;
MatProps.Einf=0;
MatProps.Ei=0;
MatProps.nu=0;


%% Determine whether K and G will be directly input or calculated from E 
% and nu
    
quest='directly input or calculate K and G';
KGdef=questdlg(quest,'K G definition method',...
    'Direct input','calculate from E and nu',...
    'calculate from E and nu');

switch KGdef
    case 'calculate from E and nu'
        
        %% input E and nu
        dlgtitle='Define E and nu';
        prompt={'long term elastic modulus','vector of element moduli',...
                'Long Term Poissons ratio','vector of Element nu'};
        dims=[1,50;1,50;1,50;1,50];
        Enu=inputdlg(prompt,dlgtitle,dims);
        
        Einf=str2num(cell2mat(Enu(1)));
        Ei=str2num(cell2mat(Enu(2)));
        E0=sum(Ei)+Einf;
        nuinf=str2num(cell2mat(Enu(3)));
        nui=str2num(cell2mat(Enu(4)));
        
        %% Calculate K from E and nu
        MatProps.Kinf=Einf/(3*(1-2*nuinf));
        MatProps.Ki=Ei./(3*(1-2*nui));
        MatProps.K0=MatProps.Kinf+sum(MatProps.Ki);
        
        %% Calculate G from E and nu
        MatProps.Ginf=Einf/(2*(1+nuinf));
        MatProps.Gi=Ei./(2*(1+nui));
        MatProps.G0=MatProps.Ginf+sum(MatProps.Gi);
        
    case 'Direct input'
        %% input K and G
       dlgtitle='define K and G';
       prompt={'Long term K', 'vector of elemental K',...
           'Long term G','vector of elemental G'};
       dims=[1,50;1,50;1,50;1,50];
       KG=inputdlg(prompt,dlgtitle,dims);
       
       %% define K
       MatProps.Kinf=str2double(cell2mat(KG(1)));
       MatProps.Ki=str2num(cell2mat(KG(2)));
       MatProps.K0=MatProps.Kinf+sum(MatProps.Ki);
       
       %% define G
       MatProps.Ginf=str2double(cell2mat(KG(3)));
       MatProps.Gi=str2num(cell2mat(KG(4)));
       MatProps.G0=MatProps.Ginf+sum(MatProps.Gi);
       
end

       
    case 'E and nu'
      %% input E and nu
        dlgtitle='Define E and nu';
        prompt={'long term elastic modulus','vector of element moduli',...
                'Long Term Poissons ratio','vector of Element nu'};
        dims=[1,50;1,50;1,50;1,50];
        Enu=inputdlg(prompt,dlgtitle,dims); 
        
        MatProps.Einf=str2double(cell2mat(Enu(1)));
        MatProps.Ei=str2num(cell2mat(Enu(2)));
        MatProps.E0=MatProps.Einf+sum(MatProps.Ei);
        MatProps.nu=str2double(cell2mat(Enu(3)));
        MatProps.nui=str2double(cell2mat(Enu(4)));
end 
        
%% define time constants
tau=str2num(cell2mat(inputdlg('vector of time constants (s)',...
    'time constant',[1,50])));
MatProps.tau=tau;


%% Save MatProps stucture
switch formtype
    case 'K and G'
        save(strcat(testdeg,'_KGprops.mat'),'MatProps' );
    case 'E and nu'
        save(strcat(testdeg,'_EnuProps.mat'),'MatProps');
end

        