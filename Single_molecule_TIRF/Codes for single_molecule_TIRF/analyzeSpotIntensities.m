function analyzeSpotIntensities(varargin)

    % Configure behavior
    spotRadius = 1;
    maxBackground = 0;
    subtractBackground = 1;
    onlyInitialTrajectories = 0;
    Nbleach = 10;
    trackTrajectory = 1;
    extendTrajectory = 0;
    showImages = 0;
    showIntensityProfiles = 0;
    meanLength=-1;
    backgroundMethod = 0; % 3 = automatic through segmentation, 
                          % 0 = manual per cell
                          % 1 = manual per frame
                          % 2 = automatic around spot
                          % 4 = convex hull of spots over all frames
    k0Percentile = 0.5; % For method = 0
    k0R = 4*spotRadius; % For method = 3 or 4
    excludeTrajectories = [];
    
    if mod(size(varargin),2)~=0
        error('The number of variable arguments must be even')
    end
    for i=1:2:length(varargin)
        if strcmp(varargin{i},'maxBackground')==1
            maxBackground=varargin{i+1};
        elseif strcmp(varargin{i},'spotRadius')==1
            spotRadius=varargin{i+1};
        elseif strcmp(varargin{i},'subtractBackground')==1
            subtractBackground=varargin{i+1};
        elseif strcmp(varargin{i},'onlyInitialTrajectories')==1
            onlyInitialTrajectories=varargin{i+1};
        elseif strcmp(varargin{i},'trackTrajectory')==1
            trackTrajectory=varargin{i+1};
        elseif strcmp(varargin{i},'extendTrajectory')==1
            extendTrajectory=varargin{i+1};
        elseif strcmp(varargin{i},'Nbleach')==1
            Nbleach=varargin{i+1};
        elseif strcmp(varargin{i},'showImages')==1
            showImages=varargin{i+1};
        elseif strcmp(varargin{i},'showIntensityProfiles')==1
            showIntensityProfiles=varargin{i+1};
        elseif strcmp(varargin{i},'backgroundMethod')==1
            backgroundMethod=varargin{i+1};
        elseif strcmp(varargin{i},'backgroundPercentile')==1
            k0Percentile=varargin{i+1};
        elseif strcmp(varargin{i},'backgroundRadius')==1
            k0R=varargin{i+1};
        elseif strcmp(varargin{i},'meanLength')==1
            meanLength=varargin{i+1};
        elseif strcmp(varargin{i},'excludeTrajectories')==1
            excludeTrajectories=varargin{i+1};
        else
            error(['Parameter ' varargin{i} ' is unrecognized'])
        end
    end
    disp(['spotRadius =              ' num2str(spotRadius)])
    disp(['maxBackground =           ' num2str(maxBackground)])
    disp(['subtractBackground =      ' num2str(subtractBackground)])
    disp(['onlyInitialTrajectories = ' num2str(onlyInitialTrajectories)])
    disp(['trackTrajectory =         ' num2str(trackTrajectory)])
    disp(['extendTrajectory =        ' num2str(extendTrajectory)])
    disp(['showImages =              ' num2str(showImages)])
    disp(['showIntensityProfiles =   ' num2str(showIntensityProfiles)])
    disp(['backgroundMethod =        ' num2str(backgroundMethod)])
    disp(['k0Percentile =            ' num2str(k0Percentile)])
    disp(['k0R =                     ' num2str(k0R)])
    disp(['excludeTrajectories =     ' num2str(excludeTrajectories)])
    
    warning('off','all')
    fnTracksDir=fullfile('results','TrackingPackage','tracks');
    fnRootImage=fullfile('videoSeq');
    fnTrack=fullfile(fnTracksDir,'Channel_1_tracking_result.mat');
    load(fnTrack);

    fnMask=fullfile('mask.tif');
    mask=[];
    numberOfRegionsInMask=1;
    if exist(fnMask,'file')
        mask=imread(fnMask);
        numberOfRegionsInMask=max(mask(:));
    end

    fhMean=fopen(fullfile(fnTracksDir,'meanSpotIntensities.txt'),'w');
    fhFrame=fopen(fullfile(fnTracksDir,'spotIntensitiesByFrame.txt'),'w');
    fnDiary=fullfile(fnTracksDir,'log.txt');
    diary(fnDiary);
    
    % Check first and last frames with events
    veryFirstImage = 1e38;
    veryLastImage = -1e38;
    for i=1:length(tracksFinal)
        idx=tracksFinal(i).seqOfEvents;
        firstImage=idx(1,1)-1;
        lastImage=idx(2,1)-1;
        if firstImage < veryFirstImage
            veryFirstImage = firstImage;
        end
        if lastImage > veryLastImage
            veryLastImage = lastImage;
        end
    end
    fnFrames = dir(fullfile('videoseq','video*tif'));
    lastImageInDisk=size(fnFrames,1)-1;
    
    % Construct a mask for this cell
    if backgroundMethod==3
        cellMask=[];
        disp('Constructing cell mask')
        for img=veryFirstImage:veryLastImage        
            fnImg=fullfile(fnRootImage,sprintf('video%04d.tif',img));
            I=imread(fnImg);

            if maxBackground>0
                I(I>maxBackground) = maxBackground/10;
            end

            % Segment the cell
            level = graythresh(I);
            BW = im2bw(I,level);
            if size(cellMask,1)==0
                cellMask=BW;
            else
                cellMask=cellMask | BW;
            end
        end
        cellMask=imopen(cellMask,strel('disk',3));
        cellMask=bwlabel(cellMask);
        M=max(cellMask(:));
        count=zeros(M,1);
        for m=1:M
            count(m)=sum(cellMask(:)==m);
        end
        [M,im]=max(count);
        cellMask=cellMask==im;
    end

    bgX = k0R*cosd(0:45:315);
    bgY = k0R*sind(0:45:315);
    
    if backgroundMethod==4
        fnImg=fullfile(fnRootImage,sprintf('video%04d.tif',veryFirstImage));
        I=imread(fnImg);
        cellMask=uint8(I*0);
        X=[];
        Y=[];
        for i=1:length(tracksFinal)
            if any(excludeTrajectories==i)
                continue
            end
            x=tracksFinal(i).tracksCoordAmpCG(1:8:end); % x is in nm
            y=tracksFinal(i).tracksCoordAmpCG(2:8:end); % y is in nm
            X=[X x];
            Y=[Y y];
        end
        
        k=convhull(X,Y);
        % plot(X,Y,'.'); hold on; plot(X(k),Y(k))
        for y=1:size(cellMask,1)
            for x=1:size(cellMask,2)
                if inpolygon(x,y,X(k),Y(k))
                    cellMask(y,x)=1;
                end
            end
        end
        cellMask = cellMask==1;
    end
    
    if backgroundMethod==0 % Manual, per cell
        fnImg=fullfile(fnRootImage,sprintf('video%04d.tif',veryFirstImage));
        I=imread(fnImg);
        figure(1)
        subplot(111);
        imagesc(I); colormap gray;
        disp('Select 8 background points to track along the series');
        [x,y]=ginput(8);

        allBrightness=[];
        for j=1:length(x)
            xx=round(x(j));
            yy=round(y(j));

            for img=veryFirstImage:veryLastImage
                fnImg=fullfile(fnRootImage,sprintf('video%04d.tif',img));
                I=imread(fnImg);
                patch=double(I((yy-spotRadius:yy+spotRadius),(xx-spotRadius:xx+spotRadius)));
                allBrightness=[allBrightness; mean(patch(:))];
            end
        end
        k0 = quantile(allBrightness,0.95);
    end
    
    % Process all spots
    disp('Processing trajectories')
    for i=1:length(tracksFinal)
        if any(excludeTrajectories==i)
            continue
        end
        x=tracksFinal(i).tracksCoordAmpCG(1:8:end); % x is in nm
        y=tracksFinal(i).tracksCoordAmpCG(2:8:end); % y is in nm
        I=tracksFinal(i).tracksCoordAmpCG(4:8:end);
        idx=tracksFinal(i).seqOfEvents;
        firstImage=idx(1,1)-1;
        lastImage=idx(2,1)-1;
        if firstImage ~= veryFirstImage && onlyInitialTrajectories
            continue
        end
        
        jj=1;
        intensitiesAlongTrajectory=[];
        meanIntensityAlongTrajectory=[];
        rawMeanIntensityAlongTrajectory=[];
        k0AlongTrajectory=[];
        validCoordinates=[];
        maskRegionAlongTrajectory=zeros(numberOfRegionsInMask+1,1);
        lastImageForTrack=lastImage;
        if extendTrajectory
            lastImageForTrack=lastImageInDisk;
        end
        for img=firstImage:lastImageForTrack
            fnImg=fullfile(fnRootImage,sprintf('video%04d.tif',img));
            I=double(imread(fnImg));
            
            if ~isnan(x(jj)) && ~isnan(y(jj))
                % Estimate background for this frame
                if backgroundMethod==3 || backgroundMethod==4 % Automatic through segmentation
                    k0=quantile(I(cellMask(:)),k0Percentile);
                elseif backgroundMethod==1 % Manual, per frame
                    figure(1)
                    subplot(111);
                    imagesc(I); colormap gray;
                    disp('Select 8 background points in this frame');
                    [xi,yi]=ginput(8);

                    allBrightness=[];
                    for j=1:length(xi)
                        xx=round(xi(j));
                        yy=round(yi(j));

                        patch=double(I((yy-spotRadius:yy+spotRadius),(xx-spotRadius:xx+spotRadius)));
                        allBrightness=[allBrightness; mean(patch(:))];
                    end
                    k0 = quantile(allBrightness,0.95);
                elseif backgroundMethod==2 % Automatic through radius
                    xi=round(x(jj)+bgX);
                    yi=round(y(jj)+bgY);
                    allBrightness=[];
                    for j=1:length(xi)
                        xx=round(xi(j));
                        yy=round(yi(j));

                        if xx-spotRadius>=1 && xx+spotRadius<=size(I,2) && ...
                           yy-spotRadius>=1 && yy+spotRadius<=size(I,1)
                            patch=double(I((yy-spotRadius:yy+spotRadius),(xx-spotRadius:xx+spotRadius)));
                            allBrightness=[allBrightness; mean(patch(:))];
                        end
                    end
                    k0 = quantile(allBrightness,0.95);
                end
                
                % Measure the intensities around the spot
                xx=round(x(jj));
                yy=round(y(jj));
                if xx-spotRadius>=1 && xx+spotRadius<=size(I,2) && yy-spotRadius>=1 && yy+spotRadius<=size(I,1)
                    patch=double(I((yy-spotRadius:yy+spotRadius),(xx-spotRadius:xx+spotRadius)));
                    rawMeanIntensityAlongTrajectory=[rawMeanIntensityAlongTrajectory; mean(patch(:))];
                    if subtractBackground
                        patch=patch-k0;
                    end
                    k0AlongTrajectory=[k0AlongTrajectory; k0];
                    intensitiesAlongTrajectory=[intensitiesAlongTrajectory; patch(:)];
                    meanIntensityAlongTrajectory=[meanIntensityAlongTrajectory; mean(patch(:))];
                    validCoordinates=[validCoordinates img];
                    maxSpot = max(patch(:));

                    if size(mask,1)>0
                        maskRegion=mask(yy,xx);
                        maskRegionAlongTrajectory(maskRegion+1)=maskRegionAlongTrajectory(maskRegion+1)+1;
                    else
                        maskRegion=0;
                    end
        
                    disp(['    spot=' num2str(i) ' frame=' num2str(img) ...
                          ' x=' num2str(xx) ' y=' num2str(yy) ...
                          ' k0=' num2str(k0) ...
                          ' spotRawIntensity=' num2str(rawMeanIntensityAlongTrajectory(end)) ...
                          ' spotCorrectedIntensity=' num2str(meanIntensityAlongTrajectory(end)) ...
                          ' maxCorrectedSpotIntensity=' num2str(maxSpot) ...
                          ' maskRegion=' num2str(maskRegion)])
                    fprintf(fhFrame,'%d %d %f %f %f %d %f %f %f\n',i,img,...
                        mean(rawMeanIntensityAlongTrajectory(end)),meanIntensityAlongTrajectory(end),...
                        maxSpot,maskRegion,xx,yy,k0);
                    
                    if showImages && backgroundMethod==3
                      figure(1);
                      subplot(131); imagesc(I); colormap gray
                      subplot(132); imagesc(cellMask);
                      hold on
                      plot(x(jj),y(jj),'o');
                      subplot(133);
                      I(I<k0)=max(I(:));
                      imagesc(I);
                      text(round(x(jj)),round(y(jj)),'X','EdgeColor',[.7 .9 .7]);
                      disp('Press any key');
                      pause
                    end
                end
            end
            
            if trackTrajectory && img<lastImage
                jj=jj+1;
            end
        end
        
        % Check photobleaching
        if length(meanIntensityAlongTrajectory)>2*Nbleach
            x1=meanIntensityAlongTrajectory(1:Nbleach);
            x2=meanIntensityAlongTrajectory(end-Nbleach:end);
            [H,pval]=ttest2(x1,x2,'tail','right');
            if pval<0.05
                disp('   ******* Possible photobleaching')
            end
        end
        
        % Write summary for this trajectory
        [M,iM]=max(maskRegionAlongTrajectory);
        if meanLength==-1
            meanIntensity = mean(intensitiesAlongTrajectory(:));
        else
            lengthToUse = min(meanLength,length(intensitiesAlongTrajectory(:)));
            meanIntensity = mean(intensitiesAlongTrajectory(1:lengthToUse));
        end
        disp([' spot=' num2str(i) ' meanCorrectedSpotIntensity along frames=' num2str(meanIntensity) ' majoritarian region=' num2str(iM)])
        fprintf(fhMean,'%d %f %d\n',i,meanIntensity,iM);

        if showIntensityProfiles
          figure(1);
          subplot(111);
          if trackTrajectory
            plot(validCoordinates,rawMeanIntensityAlongTrajectory,validCoordinates,k0AlongTrajectory);
            xlabel('Frame number')
          else
            plot(1:length(rawMeanIntensityAlongTrajectory),rawMeanIntensityAlongTrajectory,...
                 1:length(k0AlongTrajectory),k0AlongTrajectory);
          end
          legend('Spot mean raw intensity','Background intensity')
          disp('Press any key');
          pause
        end
    end
    fclose(fhMean);
    fclose(fhFrame);
    diary off


    