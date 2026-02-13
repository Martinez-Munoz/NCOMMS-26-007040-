function gatherDiffusionAndIntensity(varargin)
    if size(varargin)>0
        diffusionSuffix=varargin{1};
    else
        diffusionSuffix='';
    end

    fnTrackDir = fullfile('results','TrackingPackage','tracks');
    fnDiffusion=fullfile(fnTrackDir,['diffusionCoefficients' diffusionSuffix '.txt']);
    fnIntensities=fullfile(fnTrackDir,'meanSpotIntensities.txt');
    fnOut=fullfile(fnTrackDir,['diffusionCoefficientsAndIntensities' diffusionSuffix '.txt']);
    fnLog=fullfile(fnTrackDir,['log_diffusionCoefficientsAndIntensities' diffusionSuffix '.txt']);
    
    diffusions=textread(fnDiffusion);
    intensities=textread(fnIntensities);
    fhOut = fopen(fnOut,'w');
    
    diary(fnLog)
    for i=1:size(diffusions,1)
        spotIdx = diffusions(i,1);
        diffusionCoefficient = diffusions(i,2);
        if size(diffusions,2)==3
           maskRegion = diffusions(i,3);
        else
            maskRegion = 1;
        end
        
        idx = find(intensities(:,1)==spotIdx);
        if ~isempty(idx)
            intensity = intensities(idx,2);
            disp(['spot= ' num2str(spotIdx) ...
                  ' diffusion= ' num2str(diffusionCoefficient) ...
                  ' spotCorrectedIntensity= ' num2str(intensity) ...
                  ' maskRegion=' num2str(maskRegion)])
            fprintf(fhOut,'%d %f %f %d\n',spotIdx,diffusionCoefficient,...
                intensity,maskRegion);
        end
    end
    diary off
    fclose(fhOut);
end
