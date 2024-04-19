function [DesigVec,SaveDirVec,FileNameVec]=...
    func_BatchProcessVarsV2(SweepDesig,MainDir)
% This script is written to produce Vectors with the necessary variables to
%process parametric sweeps with up to 3 variables in a batch instead of
%individually

%Author: Andrew Matejunas
%Date: 2022-08-25

%Change Log:

%% Choose source of data structures
quest='Data location source?';
DataMethod=questdlg(quest,'dataMethod','Manual','From File','Manual');

prompt='input number of sweep variables (1-3)';
NumVar=str2num(char(inputdlg(prompt)));

switch DataMethod
    case 'From File'
        [SweepFile,SweepPath]=uigetfile('*.mat','Select sweep data file');
        load(strcat(SweepPath,'/',SweepFile));
        
    case 'Manual'
        prompt={' Function for Sweep Variable 1',...
            'Function for Sweep Variable 2',...
            'Function for Sweep Variable 3'};
        SweepVars=inputdlg(prompt);
        
        SweepVec1=eval(char(SweepVars{1}));
        SweepChar1=num2str(SweepVec1');
        
        if NumVar>=2
            SweepVec2=eval(char(SweepVars{2}));
            SweepChar2=num2str(SweepVec2');
        end
        
        if isempty(SweepVars{3})==3
            SweepVec3=eval(char(SweepVars{3}));
            SweepChar3=num2str(SweepVec3');
        end
        
        if NumVar==3
            DesigArray=cell(length(SweepVec1),length(SweepVec2),length(SweepVec3));
            SaveArray=cell(length(SweepVec1),length(SweepVec2),length(SweepVec3));
            FileNameArray=cell(length(SweepVec1),length(SweepVec2),...
                length(SweepVec3));
            prompt='Input characters of filename after the end of the test designation';
            EndChars=char(inputdlg(prompt));
            for m=1:length(SweepVec1)
                for n=1:length(SweepVec2)
                    for k=1:length(SweepVec3)
                        DesigArray{m,n,k}=strcat(SweepDesig,'_',SweepChar1(m),'_',...
                            SweepChar2(n),'_',SweepChar3(k));
                        SaveArray{m,n,k}=strcat(MainDir,'\',DesigArray{m,n,k},...
                            '\Matlab\');
                        FileNameArray{m,n,k}=strcat(SaveArray{m,n,k},...
                            DesigArray{m,n,k},EndChars);
                    end
                end
            end
        end
        
        if NumVar==2
            DesigArray=cell(length(SweepVec1),length(SweepVec2));
            SaveArray=cell(length(SweepVec1),length(SweepVec2));
            FileNameArray=cell(length(SweepVec1),length(SweepVec2));
            prompt='Input characters of filename after the end of the test designation';
            EndChars=char(inputdlg(prompt));
            for m=1:length(SweepVec1)
                for n=1:length(SweepVec2)
                    DesigArray{m,n}=strcat(SweepDesig,'_',SweepChar1(m),'_',...
                        SweepChar2(n));
                    SaveArray{m,n}=strcat(MainDir,'\',DesigArray{m,n},...
                        '\Matlab\');
                    FileNameArray{m,n}=strcat(SaveArray{m,n},...
                        DesigArray{m,n},EndChars);
                end
            end
        end
        
        if NumVar==1
            DesigArray=cell(length(SweepVec1));
            SaveArray=cell(length(SweepVec1));
            FileNameArray=cell(length(SweepVec1));
            prompt='Input characters of filename after the end of the test designation';
            EndChars=char(inputdlg(prompt));
            for m=1:length(SweepVec1)
                DesigArray{m}=strcat(SweepDesig,'_',SweepChar1(m));
                SaveArray{m}=strcat(MainDir,'\',DesigArray{m},...
                    '\Matlab\');
                FileNameArray{m}=strcat(SaveArray{m},...
                    DesigArray{m},EndChars);
                
            end
        end
        
        DesigVec=reshape(DesigArray,[numel(DesigArray),1]);
        SaveDirVec=reshape(SaveArray,[numel(SaveArray),1]);
        FileNameVec=reshape(FileNameArray,[numel(FileNameArray),1]);
        
        save(strcat(MainDir,'\',SweepDesig,'_SweepVecs.mat'))
end
end