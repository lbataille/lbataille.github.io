function lGO = convertReg2Grids(resultsP,filt)
    if nargin ==1
        filt =false;
    end
    lGO = [];
    fFN =  {"CCx";...
            "CGC";...
            "CC0";...
            "CDC";...
            "UPVT0";...
            "UPVTSEN";...
            "TIC";...
            "R2";...
            "RMSE"}; 
    for i = 1:numel(resultsP)
        mapSFi = generateGridFields(resultsP(i).xy,resultsP(i).fFCOVER,resultsP(i).fFCOVER.Properties.VariableNames);
        
        if ~any(strcmp(resultsP(i).fFCOVER.Properties.VariableNames,'Yield'))
            mapSYi = generateGridFields(resultsP(i).xy,resultsP(i).ACRS,{'Yield'});
        else
            disp("coucou")
            mapSYi = generateGridFields(resultsP(i).xy,resultsP(i).fFCOVER,{'Yield'});
            
        end
        lGOi = [];
        R2i = processStaticMaps(mapSFi,"R2");
        for j=1:numel(fFN)
            
            if filt
                goFij = processStaticMaps(mapSFi,fFN{j});
                mapSFij = mapSFi.(fFN{j});
                mapSFij(R2i.Z<0.8) = nan;
                mapSFi.(fFN{j}) = mapSFij;
                goFij0 = processStaticMaps1(mapSFi,fFN{j});                 
                goFij.Z(R2i.Z <0.8) = goFij0.Z(R2i.Z <0.8); 
            else
                goFij = processStaticMaps(mapSFi,fFN{j});
            end
            lGOi.(fFN{j}) = goFij;

        end
        if filt
            goYij = processStaticMaps(mapSYi,"Yield");
            mapSYij = mapSFi.(fFN{j});
            mapSYij(R2i.Z<0.8) = nan;
            mapSYi.(fFN{j}) = mapSYij;
            goYij0 = processStaticMaps1(mapSYi,"Yield");                 
            goYij.Z(R2i.Z <0.8) = goYij0.Z(R2i.Z <0.8); 
        else
            goYij = processStaticMaps(mapSYi,"Yield");
        end        


        lGOi.Yield = goYij;
        lGO = [lGO;lGOi];
    end
end

function map = generateGridFields(xyP,data,fields)
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
    for i = 1:length(fields)
        if (class(data.(fields{i}))=="datetime")
            map.(fields{i}) = NaT(size(X));
        else
            map.(fields{i}) = nan(size(X));
        end
    end
    
    for i = 1:length(fields)
       map.(fields{i})(ij(:,2)+(ij(:,1)-1)*(size(X,1))) = data.(fields{i}); 
    end
           
end


function goG = processStaticMaps(mapG,field)
    goG = GRIDobj(mapG.X,mapG.Y,mapG.(field));
    goGZm = nanmedfilt2(goG.Z,[5,5]);
    dbw = bwdist(isnan(goG.Z),'quasi-euclidean');
    perci = prctile(goG.Z-goGZm,99,'all');
    goG.Z(((goG.Z-goGZm > perci))) = goGZm(((goG.Z-goGZm > perci)));
    goG.Z(dbw*goG.cellsize<20) = nan;
end


function goG = processStaticMaps1(mapG,field)
    goG = GRIDobj(mapG.X,mapG.Y,mapG.(field));
    goGZm = nanmedfilt2(goG.Z,[7,7]);
    dbw = bwdist(isnan(goG.Z),'quasi-euclidean');
    perci = prctile(goG.Z-goGZm,99,'all');
    goG.Z = goGZm;
%     goG.Z(dbw<2) = nan;
end
