filelistdata2 = dir('Z:\Maria Frechou\20210228Exp41DREADDsRetroAAV\2*\*\DG\cleaned_traces_data.npy')

%filelistdata2 = dir('\\data.einsteinmed.org\users\Maria Frechou\Ex9\*\CA1\cleaned_traces_data.npy'); dir('\\data.einsteinmed.org\users\Maria Frechou\Ex9\*\mossy\cleaned_traces_data.npy');

MaxIterations=10000; %set number of iterations

%TreadmillVersion=1; % Choose 0 for old treadmill (pre 01Jan2021) datasets A and B , choose 1 for new treadmill, dataset C
FrameoutPresent=0; % If == 0 there is not frameout data, if == 1 framout data is present



for files=1:size(filelistdata2,1)
    
    TreadmillVersion=1;
    
    levels = strfind(filelistdata2(files).folder, '\');
    
    
    %set where each of the concatenated movies ends, so the transitions can be
    %correctly evaluated

    mov1 = dir([filelistdata2(files).folder(1:levels(end)) 'DG\Chan*_001*.tif']);
    mov2 = dir([filelistdata2(files).folder(1:levels(end)) 'DG_000\Chan*_001*.tif']);
    mov3 = dir([filelistdata2(files).folder(1:levels(end)) 'DG_001\Chan*_001*.tif']);
    
    movie1End=numel(mov1);
    movie2End=numel(mov2);
    movie3End=numel(mov3);
    movieEnds=[movie1End; movie2End; movie3End];
    
    rpmdata=readNPY([filelistdata2(files).folder '\mset.npy']); %read RPM texture data
    rfiddata=readNPY([filelistdata2(files).folder '\tset.npy']); %read RFID texture data
    traces=readNPY([filelistdata2(files).folder '\cleaned_traces_data.npy']); %read deconvolved traces
    
    rfiddata=rfiddata'; %may be needed if loaded from .npy file
    rpmdata=rpmdata';   %may be needed if loaded from .npy file
    
    if TreadmillVersion==0
    
        detectionThreshold = 300;
    
    end

    if TreadmillVersion==1
    
        detectionThreshold = 700;
        rpmdata = rpmdata-1.372462080290498;
    end
    
    %imagingmeta=xml2struct([filelistdata2(files).folder '\Experiment.xml']);
    
    % Analog input data from RPM and RFID sensors is slighly longer than
    % Ca2+ movies (recording only stops after movie is saved)
    % This section removes final segment of these recordings
    
%     frate=str2num(imagingmeta.ThorImageExperiment.LSM.Attributes.frameRate); %determine Ca2+ imaging frame rate
%     mlength=size(traces,2)/frate; % caculate Ca2+ movie length in seconds
%     lastrecording=ceil(mlength*30000); %last datapoint in RPM and RFID data @30kHz that matches Ca2+ movie
%     
%     rfiddata(lastrecording:end)=[]; %trim rfid data
%     rpmdata(lastrecording:end)=[]; %trim RPM data
    
    
    %First find transitions and resample RFID and RPM data to size of movies (traces)
    
    cd(filelistdata2(files).folder);
    
    %Find transitions in original 30kHz data as this is more robust
    
    transitions = findchangepts(rfiddata,'Statistic', 'mean', 'MinThreshold',detectionThreshold); %May need to tweak threshold, 300 works for 30kHz data
    
    %Plot transitions before trimming data
    
    hold off
    yvaltran=rfiddata(transitions);
    plot(rfiddata);
    hold on
    plot(transitions,yvaltran,'+r');
    savefig('transitions-pretrim.fig');
    saveas(gcf,'transitions-pretrim.jpg');
    
    %Then convert these points to timepoints in movie frame units
    
    movieTransitions=transitions.*(size(traces,2)/size(rfiddata,2));
    movieTransitions=ceil(movieTransitions);
    
    %Conversely, convert movieEnds to digitizer rpm/rfid time unids
    
    digitizerEnds = movieEnds./(size(traces,2)/size(rfiddata,2));
    digitizerEnds =round(digitizerEnds);
    if digitizerEnds(end)>size(rfiddata,2)
        digitizerEnds(end)=size(rfiddata,2);
    end
    
    
    
%     %Find index of zone transitions that occur between concatenated movies 
%      
%     transBetweenMovies=find(abs(movieTransitions-movieEnds)<=1);
%     [~,transBetweenMovies]=ind2sub([size(movieEnds,1), size(transitions,2)],transBetweenMovies);
    

    
    %Remove initial part of each movie (except for first) until a zone transition is reached (as
    %cannot be sure of the position of the animal until that)
    
        
       for movieBreak=(size(movieEnds,1)-1):-1:1
           transitionAfterNewMovie = movieTransitions-movieEnds(movieBreak);
           transitionAfterNewMovie(transitionAfterNewMovie < 0)=NaN; %set threshold as 0 frames but use ceil above because of rounding errors
           [~,transitionAfterNewMovie]=min(transitionAfterNewMovie);
           traces(:,movieEnds(movieBreak):movieTransitions(transitionAfterNewMovie))=[]; 
           rfiddata(digitizerEnds(movieBreak):transitions(transitionAfterNewMovie))=[];
           rpmdata(digitizerEnds(movieBreak):transitions(transitionAfterNewMovie))=[];       
       end
       
     %Also remove first segment of first movie (until first transition is
     %reached)
     
     traces(:,1:movieTransitions(1))=[];
     rfiddata(1:transitions(1))=[];
     rpmdata(1:transitions(1))=[];
     
     %Now, make a new list of transitions in the trimmed traces (use 30KHz version) 
    
    transitions = findchangepts(rfiddata,'Statistic', 'mean', 'MinThreshold',detectionThreshold); %May need to tweak threshold, 300 works for 30kHz data
    
    %Plot transitions with 30 kHz data before trimmming movie 
    
    hold off
    yvaltran=rfiddata(transitions);
    plot(rfiddata);
    hold on
    plot(transitions,yvaltran,'+r');
    savefig('transitions-pretrim30k.fig');
    saveas(gcf,'transitions-pretrim30k.jpg');
       
    
    %Then convert these points to timepoints in traces file
    
    transitions=transitions.*(size(traces,2)/size(rfiddata,2));
    transitions=round(transitions);
    
    rfiddata=imresize(rfiddata(1,:),[1,size(traces,2)]);
    rpmdata=imresize(rpmdata(1,:),[1,size(traces,2)]);
    
    
    %Plot transitions after trimming data
    
    hold off
    yvaltran=rfiddata(transitions);
    plot(rfiddata);
    hold on
    plot(transitions,yvaltran,'+r');
    savefig('transitions-posttrim.fig');
    saveas(gcf,'transitions-posttrim.jpg');
    
    %Then make a list of when the animal is in each zone...
    nzones = size(transitions,2)+1;
  
    zoneMeanVolt=zeros(nzones,1);
    beltZone=zeros(nzones,1);
    
    %...for this find mean voltage of each zone...
    
    %...and attribute each mean voltage to a texture
%     if TreadmillVersion==0
%         for zone=1:1:nzones
%             if (round(zoneMeanVolt(zone)-0.60)==0)
%                 beltZone(zone)=1;
%             elseif (round(zoneMeanVolt(zone)-1.8)==0)
%                 beltZone(zone)=2;
%             elseif (round(zoneMeanVolt(zone)-2.3)==0)
%                 beltZone(zone)=3;
%             elseif (round(zoneMeanVolt(zone)-1.2)==0)
%                 beltZone(zone)=4;
%             end
%         end
%     end
%     if TreadmillVersion==1
%           for zone=1:1:nzones
%             if (round(zoneMeanVolt(zone)-0.60)==0)
%                 beltZone(zone)=1;
%             elseif (round(zoneMeanVolt(zone)-1.2)==0)
%                 beltZone(zone)=2;
%             elseif (round(zoneMeanVolt(zone)-2.4)==0)
%                 beltZone(zone)=3;
%             elseif (round(zoneMeanVolt(zone)-1.8)==0)
%                 beltZone(zone)=4;
%             end
%           end
%     end
%     transitions(2,:) = transitions; % time points of transitions
%     transitions(1,:) = beltZone(2:end);
    
    for zone=1:1:nzones
        
        if (1<zone)&&(zone<nzones)
            zoneMeanVolt(zone)=mean(rfiddata(transitions(zone-1):transitions(zone)));
        elseif zone==1
            zoneMeanVolt(zone)=mean(rfiddata(1:transitions(zone)));
        elseif zone==nzones
            zoneMeanVolt(zone)=mean(rfiddata(transitions(zone-1):end));
        end    
    end
    
    %...and attribute each mean voltage to a texture
    % -0.6, -1.2, -2.4, -1.8
    
    for zone=1:1:nzones
        if (round(zoneMeanVolt(zone)-0.60)==0)
            beltZone(zone)=1;
        elseif (round(zoneMeanVolt(zone)-1.8)==0)
            beltZone(zone)=2;
        elseif (round(zoneMeanVolt(zone)-2.3)==0)
            beltZone(zone)=3;
        elseif (round(zoneMeanVolt(zone)-1.2)==0)
            beltZone(zone)=4;
        
        end
    end
    
    %Create variable with start and end of each zone
    
    zoneStartEnd = zeros(nzones,2);
    
    for zone=1:1:nzones
        
        if (1<zone)&&(zone<nzones)
            zoneStartEnd(zone,:)=[transitions(zone-1),(transitions(zone)-1)]; % subtract 1 so that border point isn't accounted twice
        elseif zone==1
            zoneStartEnd(zone,:)=[1,(transitions(zone)-1)];
        elseif zone==nzones
            zoneStartEnd(zone,:)=[transitions(zone-1),size(rfiddata,2)];
        end    
    end
    
          
    
    %Count number of times that mouse travels through each of 4 textures
    
        textureone=zeros(sum(beltZone==1),1);
        texturetwo=zeros(sum(beltZone==2),1);
        texturethree=zeros(sum(beltZone==3),1);
        texturefour=zeros(sum(beltZone==4),1);
    
    %Find distance for each trip through each texture
    for zone=1:nzones
        
        distmeasure(zone)=sum(rpmdata(zoneStartEnd(zone,1):zoneStartEnd(zone,2))); %Find distance for each trip through each texture
    end
    
    textureone = distmeasure(beltZone==1); %distance traveled in texture one
    texturetwo = distmeasure(beltZone==2); %distance traveled in texture two
    texturethree = distmeasure(beltZone==3); %distance traveled in texture three
    texturefour = distmeasure(beltZone==4);%distance traveled in texture four
    
    %Find number of data points for each zone
    
    zonesize=zeros(nzones,1); %initialize variable
    
    for zone=1:nzones
        
        zonesize(zone,1)=zoneStartEnd(zone,2)-zoneStartEnd(zone,1); % number of data points for each trip through each texture
    end
    
    zonesizeone= zonesize(beltZone==1); %number of data points traveled in texture one
    zonesizetwo = zonesize(beltZone==2); %number of data points traveled in texture two
    zonesizethree = zonesize(beltZone==3); %number of data points traveled in texture three
    zonesizefour = zonesize(beltZone==4);%number of data points traveled in texture four
    
    %Find cumulative distance along each zone, this needs to be a cell array
    %as each zone is a different size
    
    distancex{nzones,1}=[]; %initialize cell array
    
    %Create x (distance) for each zone crossing
    
    for zone=1:nzones
        
        distancex{zone}=cumsum(rpmdata(zoneStartEnd(zone,1):zoneStartEnd(zone,2))); %integrate distance for each trip through each texture
    end
    
    %and organize this data by zones and laps
    
    xone=cell(size(textureone));
    newzone=1;
    for zone=1:1:nzones
        if beltZone(zone)==1
            xone{newzone}= distancex{zone}; % x positions of data points in texture one
            newzone=newzone+1;
        end
    end
    
    xtwo=cell(size(texturetwo));
    newzone=1;
    for zone=1:1:nzones
        if beltZone(zone)==2
            xtwo{newzone}= distancex{zone}; % x positions of data points in texture two
            newzone=newzone+1;
        end
    end
    
    xthree=cell(size(texturethree));
    newzone=1;
    for zone=1:1:nzones
        if beltZone(zone)==3
            xthree{newzone}= distancex{zone}; % x positions of data points in texture three
            newzone=newzone+1;
        end
    end
    
    xfour=cell(size(texturefour));
    newzone=1;
    for zone=1:1:nzones
        if beltZone(zone)==4
            xfour{newzone}= distancex{zone}; % x positions of data points in texture four
            newzone=newzone+1;
        end
    end
    
    %Then we build histograms of the time spent at each position
    %first define bin edges
    
%     edgesone=(0:0.15:max(textureone));
%     edgestwo=(0:0.15:max(texturetwo));
%     edgesthree=(0:0.15:max(texturethree));
%     edgesfour=(0:0.15:max(texturefour));
%     
%     edgesone(end)=max(textureone);
%     edgestwo(end)=max(texturetwo);
%     edgesthree(end)=max(texturethree);
%     edgesfour(end)=max(texturefour);
    
    edgesone=20;
    edgestwo=20;
    edgesthree=20;
    edgesfour=20;
    
    %then build histograms of time spent
    
    xone=xone'; %transposing may be needed due to inconsistencies in data structure of .npy vs .h5 files
    xtwo=xtwo';
    xthree=xthree';
    xfour=xfour';
    
    [histtimeone,~]=histcounts(cell2mat(cellfun(@transpose,xone,'UniformOutput',false)),edgesone);
    [histtimetwo,~]=histcounts(cell2mat(cellfun(@transpose,xtwo,'UniformOutput',false)),edgestwo);
    [histtimethree,~]=histcounts(cell2mat(cellfun(@transpose,xthree,'UniformOutput',false)),edgesthree);
    [histtimefour,~]=histcounts(cell2mat(cellfun(@transpose,xfour,'UniformOutput',false)),edgesfour);
    
    %--------------------------------------------------------
    %From this point on analysis is cell by cell
    %--------------------------------------------------------
    
    logTraces=logical(traces); %make traces into logical values


    %initialize output structure

    CellTuning{size(traces,1),7}=[];
    
    for neuron = 1:1:size(traces,1)
        
        [firePosOne, firePosTwo, firePosThree, firePosFour] = placeFiring2(nzones, beltZone, zoneStartEnd, distancex, logTraces(neuron,:));
        
        
        %create histogram of firing positions
        [histone,~]=histcounts(cell2mat(cellfun(@transpose,firePosOne,'UniformOutput',false)),edgesone);
        [histtwo,~]=histcounts(cell2mat(cellfun(@transpose,firePosTwo,'UniformOutput',false)),edgestwo);
        [histthree,~]=histcounts(cell2mat(cellfun(@transpose,firePosThree,'UniformOutput',false)),edgesthree);
        [histfour,~]=histcounts(cell2mat(cellfun(@transpose,firePosFour,'UniformOutput',false)),edgesfour);
        
        %and normalize by the time spent in each position from the histogram
        %calculated above
        
        norhistone=histone./histtimeone;
        norhisttwo=histtwo./histtimetwo;
        norhistthree=histthree./histtimethree;
        norhistfour=histfour./histtimefour;
        
        norhistall=[norhistone, norhisttwo, norhistthree, norhistfour];
        norhistall=norhistall/sum(norhistall);
        
        %Create lookup table for the center of histogram bins from 0 to 2pi (i.e.
        %commonly referred at theta, whereas norhistall is rho), this could be
        %moved to place above since independent of cell
        
        polangle=linspace(0,2*pi, size(norhistall,2));
        zoneborders=[size(norhistone,2), size(norhistone,2)+size(norhisttwo,2), size(norhistone,2)+size(norhisttwo,2)+size(norhistthree,2)];
        
        %optional plot of tuning for each cell
        
        %polarplot(polangle, norhistall)
        
        %add all components
        
        tuningphasor=0;
        for i = 1:size(norhistall,2)
            
            tuningphasor=tuningphasor + (norhistall(i)*exp(1i*polangle(i)));
            
        end
        
        tuningidx = abs(tuningphasor);
        
        %now find if this is significant by bootstraping to create null
        %distribution
        
        randTuningDist = tuningShuffle2(nzones, beltZone, zoneStartEnd, distancex, logTraces(neuron,:), MaxIterations, edgesone, histtimeone, edgestwo, histtimetwo, edgesthree, histtimethree, edgesfour, histtimefour,polangle);
        
        SigThreshold=prctile(randTuningDist,95);
        SortedItrTuning=sort(randTuningDist);
        IndicesHigher=find(SortedItrTuning > tuningidx);
        if numel(IndicesHigher)~=0
    
            PVal=(MaxIterations-IndicesHigher(1))/MaxIterations;
        else
            PVal=1/MaxIterations;
        end
        
        %Save information of each cell: the firing histogram (rho),the
        %angles (theta), the tuning index and the alpha value for 95th percentile
        
        CellTuning{neuron,1}=norhistall;
        CellTuning{neuron,2}=polangle;
        CellTuning{neuron,3}=tuningidx;
        CellTuning{neuron,4}=SigThreshold; %threshold for alpha=0.05
        CellTuning{neuron,5}=zoneborders; %borders between zones (angles)
        CellTuning{neuron,6}=PVal; %p-value of tuningidx
        CellTuning{neuron,7}=randTuningDist; %save null distribution
        
        
    end
    
    save('CellTuning_NP.mat', 'CellTuning');
    
    
    clearvars -except filelistdata2 files movieEnds MaxIterations %clears all vars except those needed to loop through different mice
    
    
end




