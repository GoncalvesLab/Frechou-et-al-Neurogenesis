function  ItrTuningIdx = tuningShuffle2(nz, textre, zse, distx, trace, MaxItr, edg1, htime1, edg2, htime2, edg3, htime3, edg4, htime4, pangle)

%initialize surrogate tuning distribution

rng('shuffle'); %random number generator will be based on time, i.e. always different

ItrTuningIdx=zeros(1,MaxItr);


parfor itr = 1:1:MaxItr
    
    ItrReSamp = circshift(trace,randi(size(trace,2)),2);
    
    [ItrXPosFireOne, ItrXPosFireTwo, ItrXPosFireThree, ItrXPosFireFour] = placeFiring2(nz, textre, zse, distx, ItrReSamp);
    
    %then create histogram
    [ItrHistOne,~]=histcounts(cell2mat(cellfun(@transpose,ItrXPosFireOne,'UniformOutput',false)),edg1);
    [ItrHistTwo,~]=histcounts(cell2mat(cellfun(@transpose,ItrXPosFireTwo,'UniformOutput',false)),edg2);
    [ItrHistThree,~]=histcounts(cell2mat(cellfun(@transpose,ItrXPosFireThree,'UniformOutput',false)),edg3);
    [ItrHistFour,~]=histcounts(cell2mat(cellfun(@transpose,ItrXPosFireFour,'UniformOutput',false)),edg4);
    
    ItrNorHistOne=ItrHistOne./htime1;
    ItrNorHistTwo=ItrHistTwo./htime2;
    ItrNorHistThree=ItrHistThree./htime3;
    ItrNorHistFour=ItrHistFour./htime4;
    
    ItrNorHistAll=[ItrNorHistOne, ItrNorHistTwo, ItrNorHistThree, ItrNorHistFour];
    ItrNorHistAll=ItrNorHistAll/sum(ItrNorHistAll);
    
    %add all components to calculate putative tuning phasor
    
    ItrTuningPhasor=0;
    for i = 1:size(ItrNorHistAll,2)
        
        ItrTuningPhasor=ItrTuningPhasor + (ItrNorHistAll(i)*exp(1i*pangle(i)));
        
    end
    
    ItrTuningIdx(itr) = abs(ItrTuningPhasor);
    
end







