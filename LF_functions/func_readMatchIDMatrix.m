function [DataMatrix] = func_readMatchIDMatrix(DataPath,FirstFile)
%This function is intendend to read match ID matrices into Matlab for use
%in IBII experimental analysis

%Author: Andrew Matejunas
%Date Created: 2023/08/02

%Version history/change log:





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find extensions of the file
perloc=find(FirstFile=='.'); %finds periods in the file name
FileExt=FirstFile(perloc(1):end); %Determines the filed extension

%% Find Repeating part of file string
usLoc=find(FirstFile(1:perloc(1))=='_'); %Finds underscores before frame numbers
startStr=FirstFile(1:usLoc(end));

%% Sort data Files
fileStruct=dir([DataPath,'\',startStr,'*',FileExt]);
fileCell={fileStruct.name}';

sortedFiles=cell([length(fileCell)+1,1]);
sortedFiles{2} = FirstFile;  % places the first image file as first entry in sorted image cell array
for i = 2:length(sortedFiles)
    checkStr1 = [startStr,num2str(i),FileExt];   
    nn = length(num2str(i));
    if nn==1
        num = ['00' num2str(i)];
    elseif nn==2
        num = ['0' num2str(i)];
    else
        num = num2str(i);
    end
    checkStr2 = [startStr,num,FileExt];

    for j = 2:length(fileCell)       
        if strcmp(checkStr1,fileCell{j}) || strcmp(checkStr2,fileCell{j})
            sortedFiles{i} = fileCell{j};
            break
        end
    end
end
fileCell = sortedFiles;

%% Import Data From First File
TempMatrix=readmatrix(cell2mat(strcat(DataPath,'/',fileCell(2))));


DataCell=cell([length(fileCell),1]);
%DataMatrix(:,:,1)=TempMatrix;

%% Import Remaining files Into Cell Array
for k=2:length(fileCell)
    %% 
    TempFile=cell2mat(strcat(DataPath,'/',fileCell(k)));
    TempMatrix=readmatrix(TempFile);
    DataCell{k}=TempMatrix;
end

%% Record Cell Sizes
    %match ID exports change size based on the movement of the subsets
    %within the ROI (Lagrangian Description) 
CellSizes=cellfun(@size,DataCell,'UniformOutput',false);
CellSizes{1}=CellSizes{2};

DataSizes=cell2mat(CellSizes);
DataRows=min(DataSizes(:,1));
DataCol=min(DataSizes(:,2));
DataMatrix=zeros([DataRows,DataCol,length(fileCell)]);

%% Record Data Matrix
for k=2:length(fileCell)
    TempMatrix=cell2mat(DataCell(k));
    DataMatrix(:,:,k)=TempMatrix(1:DataRows,1:DataCol);
end



end