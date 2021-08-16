function map = convertData2Grid(resultsP)
    xyP = resultsP.xy;
    xyMin = min(xyP,[],1);
    xyMax = max(xyP,[],1);    
    xP0 = unique(xyP(:,1));    
    pas = xP0(2) - xP0(1);
   
    [X,Y] = meshgrid(xyMin(1):pas:xyMax(1),xyMin(2):pas:xyMax(2));
    
    ij = 1 +  floor((xyP - repmat(xyMin,[size(xyP,1),1]))./pas);
         
    
    map = [];
    map.X = X;
    map.Y = Y;
    map.ij = ij;
    map.DT =  resultsP.DT;
    map.VAL = nan([size(X),numel(map.DT),]);
    

    
    for i = 1:numel(map.DT)
       mi = map.VAL(:,:,i);
       mi(ij(:,2)+(ij(:,1)-1)*(size(X,1))) =  resultsP.dFCOVER(:,i);
       map.VAL(:,:,i) = mi; 
    end
           
end