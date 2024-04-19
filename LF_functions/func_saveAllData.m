function func_saveAllData(saveFiles,test,specimen,material,time,pos,disp,accel,strain,strainRate,stress,...
        identStiffSG,identStiffVFMan,identStiffVFOpt,identStrength,fracture)
% Author: Jared Van Blitterswyk
% PhotoDyn group, University of Southampton
% Last updated: 21/03/2017
% Edited by: Lloyd Fletcher

% Get the names of all the variables passed to the function and assign them
% to the config or the processed data file
configVarStrs = {};
processedVarStrs = {};
for i = 2:nargin
    % Check if this variable should be put in the config file
    for c = 1:length(saveFiles.configVars)
        if saveFiles.configVars(c) == i
            configVarStrs{c} = inputname(i);
            break
        end
    end
    % Check if this variable should be put in the processed data file
    for p = 1:length(saveFiles.processedVars)
        if saveFiles.processedVars(p) == i
            processedVarStrs{p} = inputname(i);
            break
        end
    end   
end

% If the variable string cell containers are empty then something is wrong
if isempty(configVarStrs) || isempty(processedVarStrs)
    error('Specified variables to save were not passed to the function.')
end

% Check if files exist - if so, delete so they can be written over
if exist([saveFiles.path,saveFiles.processedData],'file')==2
  delete([saveFiles.path,saveFiles.processedData]);
end
if exist([saveFiles.path,saveFiles.configFile],'file')==2
  delete([saveFiles.path,saveFiles.configFile]);
end
if exist([saveFiles.path,saveFiles.README],'file')==2
  delete([saveFiles.path,saveFiles.README]);
end

%--------------------------------------------------------------------------
% Save the config .mat file
% Loop over all config variables and append them to the config file
for v = 1:length(configVarStrs)
    if v == 1
        save([saveFiles.path,saveFiles.configFile],configVarStrs{v});
    else
        save([saveFiles.path,saveFiles.configFile],configVarStrs{v},'-append');
    end
end

%--------------------------------------------------------------------------
% Save processed data to .mat file
for v = 1:length(processedVarStrs)
    if v == 1
        save([saveFiles.path,saveFiles.processedData],processedVarStrs{v});
    else
        save([saveFiles.path,saveFiles.processedData],processedVarStrs{v},'-append');
    end
end

%--------------------------------------------------------------------------
% Write the README file

% Open the README file in text write mode
fID = fopen(char([saveFiles.path,saveFiles.README]), 'wt');
% Write a header at the top of the file
today = date;
fprintf(fID, strcat(today, '\n'));
fprintf(fID, '-------------------------------------------------------------\n');
fprintf(fID, 'README: Grid Method and Post Processing\n');
fprintf(fID, '-------------------------------------------------------------\n');
fprintf(fID, 'Test Description: %s \n',test.description);
fprintf(fID, '\n');
fprintf(fID, '-------------------------------------------------------------\n');
fprintf(fID, 'TEST DETAILS:\n');
fprintf(fID, '-------------------------------------------------------------\n');
fprintf(fID, 'Test Material: %s \n',material.name);
fprintf(fID, 'Projectila and Waveguide Material: %s \n',test.projMaterial);
fprintf(fID, 'Projectile Length: %.1f mm\n',test.projLength);
fprintf(fID, 'Waveguide Length: %.1f mm\n',test.wgLength);
fprintf(fID, 'Nominal Projectile Velocity: %.1f m/s\n',test.impactVel);
fprintf(fID, 'Camera Frame Rate: %.1f Mfps \n',time.frameRate*10^-6);
fprintf(fID, '\n');
% Write out the contents of the configuration file
fprintf(fID, '-------------------------------------------------------------\n');
fprintf(fID, 'Configuration File: %s\n',saveFiles.configFile);
fprintf(fID, '-------------------------------------------------------------\n');
fprintf(fID, 'Data structures contained in the file:\n');
for v = 1:length(configVarStrs)
    fprintf(fID, '\t%s\n',configVarStrs{v});
end
fprintf(fID, '\n');
% Write out the contents of the prcocessed data file
fprintf(fID, '-------------------------------------------------------------\n');
fprintf(fID, 'Processed Data File: %s\n',saveFiles.processedData);
fprintf(fID, '-------------------------------------------------------------\n');
fprintf(fID, 'Data structures contained in the file:\n');
for v = 1:length(processedVarStrs)
    fprintf(fID, '\t%s\n',processedVarStrs{v});
end
fprintf(fID, '\n');

% Close the file
fclose(fID); 
clear fID fileID;


end

