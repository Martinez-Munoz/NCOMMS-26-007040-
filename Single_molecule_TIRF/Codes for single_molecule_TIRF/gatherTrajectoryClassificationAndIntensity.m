function gatherTrajectoryClassificationAndIntensity(trajectorySuffix)
    fnTrackDir = fullfile('results','TrackingPackage','tracks');
    fnTrajectories=fullfile(fnTrackDir,['trajectoryClassification' trajectorySuffix '.txt']);
    fnDiffusions=fullfile(fnTrackDir,'diffusionCoefficients.txt');
    fnIntensities=fullfile(fnTrackDir,'meanSpotIntensities.txt');
    fnOut=fullfile(fnTrackDir,['trajectoryClassificationAndIntensities' trajectorySuffix '.txt']);
    fnLog=fullfile(fnTrackDir,['log_trajectoryClassificationAndIntensities' trajectorySuffix '.txt']);
    
    trajectories=textread(fnTrajectories);
    intensities=textread(fnIntensities);
    diffusionCoefficients=textread(fnDiffusions);
    fhOut = fopen(fnOut,'w');
    
    diary(fnLog)
    for i=1:size(trajectories,1)
        spotIdx = trajectories(i,1);
        movementType = trajectories(i,2);
        firstMoment = trajectories(i,3);
        
        idx = find(intensities(:,1)==spotIdx);
        if ~isempty(idx)
            intensity = intensities(idx,2);
            maskRegion = intensities(idx,3);
            idxC = find(diffusionCoefficients(:,1)==spotIdx);
            if ~isempty(idxC)
                diffusion=diffusionCoefficients(idxC,2);
            else
                diffusion=-1;
            end
            disp(['spot= ' num2str(spotIdx) ...
                  ' movementType= ' num2str(movementType) ...
                  ' MSSFirstMoment= ' num2str(firstMoment) ...
                  ' spotCorrectedIntensity= ' num2str(intensity) ...
                  ' diffusionCoefficient= ' num2str(diffusion) ...
                  ' maskRegion=' num2str(maskRegion)])
            fprintf(fhOut,'%d %d %f %f %f %d\n',spotIdx,movementType,firstMoment,...
                intensity,diffusion,maskRegion);
        end
    end
    diary off
    fclose(fhOut);
end