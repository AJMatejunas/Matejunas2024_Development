function [ output_args ] = func_extrapShearAtEdgesToBCs( input_args )
% Authors: Jared Van-Blitterswyk
% PhotoDyn Group, University of Southampton
% Date Created: 21/11/2018
% Date Edited: 24/1/2018

[sy,sx,st] = size(strain.xy);
ind1 = round(smoothingOpts.spatialKernal(1)/2); % position inwards from left edge for data reconstruction
ind2 = sx - ind1; % position inwards from right edge for data reconstruction
strain.sRE = strain.s;
% reconstruct edge shear strains
for i = 1:time.numFrames % loop through number of frames
    for j = 1:sy % loop through each position on the specimen
        temp1 = strain.s(j,ind1,i);
        temp2 = strain.s(j,ind2,i);
        slope1 = temp1/pos.x(ind1);
        slope2 = -temp2/(pos.x(end) - pos.x(ind2));
        for k = 1:ind1
            strain.sRE(j,k,i) = 0 + pos.x(k)*slope1; 
        end
        for m = ind2:sx
            strain.sRE(j,m,i) = temp2 + (pos.x(m)-pos.x(ind2))*slope2; 
        end
    end
end



end

