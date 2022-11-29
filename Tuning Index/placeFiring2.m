function [xposfireone, xposfiretwo, xposfirethree, xposfirefour] = placeFiring2(nZones, texture, startAndEnd, distanceX, activeindices)
    
    activepts=zeros(nZones,1);
    activeindices=logical(activeindices); %indices of active points for all zones
    
    for zone=1:nZones
        
        activepts(zone,1)=sum(activeindices(startAndEnd(zone,1):startAndEnd(zone,2))); % number of active (firing) data points for each trip through each texture
    end
    
    %number of active points in each zone
    
    activeptsone= activepts(texture==1,:); %number of active data points (i.e. time) in texture one
    activeptstwo = activepts(texture==2,:); %number of active data points (i.e. time) in texture two
    activeptsthree = activepts(texture==3,:); %number of active data points (i.e. time) in texture three
    activeptsfour = activepts(texture==4,:); %number of active data points (i.e. time) in texture four
    
    % Now get x-position of each active point
    
    xposfire{nZones,1}=[]; %initialize cell array
    
    activeindicesz{nZones,1}=[];
    
    for zone=1:nZones
        
        activeindicesz{zone}=activeindices(startAndEnd(zone,1):startAndEnd(zone,2)); %integrate distance for each trip through each texture
    end
    
    
    
    for zone=1:nZones
        
        xposfire{zone}=distanceX{zone}(activeindicesz{zone}); % x positions of active (firing) data points for each trip through each texture
    end
    
    %and organize this firing by zones and laps
    
    xposfireone=cell(size(activeptsone));
    newzone=1;
    for zone=1:1:nZones
        if texture(zone)==1
            xposfireone{newzone}= xposfire{zone}; % x positions of active data points in texture one
            newzone=newzone+1;
        end
    end
    
    xposfiretwo=cell(size(activeptstwo));
    newzone=1;
    for zone=1:1:nZones
        if texture(zone)==2
            xposfiretwo{newzone}= xposfire{zone}; % x positions of active data points in texture two
            newzone=newzone+1;
        end
    end
    
    xposfirethree=cell(size(activeptsthree));
    newzone=1;
    for zone=1:1:nZones
        if texture(zone)==3
            xposfirethree{newzone}= xposfire{zone}; % x positions of active data points in texture three
            newzone=newzone+1;
        end
    end
    
    xposfirefour=cell(size(activeptsfour));
    newzone=1;
    for zone=1:1:nZones
        if texture(zone)==4   
            xposfirefour{newzone}= xposfire{zone}; % x positions of active data points in texture four
            newzone=newzone+1;
        end
    end
    
end

    

    
    
   


