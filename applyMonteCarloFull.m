function MC = applyMonteCarloFull(resultsP,GF,pheno,nRep)
    pathOuputAC = 'C:\Users\Laurent\Documents\MATLAB\AquaCropOS_v60a\Output1\'; 
    dT = resultsP.DT;
    dFCOVER = resultsP.dFCOVER;
    idxK =  knnsearch(resultsP.xy,resultsP.xy,'K',25,'Distance','euclidean');
    perc = 0:5:100;
    paramPERC = table('Size',[numel(perc),7],...
            'VariableTypes',repmat("double",[7,1]),...
            'VariableNames',{'CCx';'CGC';'CC0';'CDC';'UPVT0';'UPVTSEN';'Yield'});   
    
    dn1 = linspace(datenum(dT(1)),datenum(dT(end)),2e3);
    gdd1 = interp1(datenum(pheno.date),pheno.GDDc,dn1);        
    funcCCI = @(modi,UPVTCi) GF.CCI(modi.CC0,modi.CGC,modi.CCx,modi.CDC,modi.UPVT0,modi.UPVTSEN,UPVTCi);
        
    MC = [];
    for i =1:size(dFCOVER,1)
        
        dFCOVER99 = prctile(dFCOVER((idxK(i,:))',:),99,1);
        dFCOVER01 = prctile(dFCOVER((idxK(i,:))',:),1,1);

        dFCOVERstd = std(dFCOVER((idxK(i,:))',:),1,'omitnan');

        nDays = size(dFCOVER01,2);
        dFCOVERMC = nan(nRep,nDays);
        for j = 1:nDays
            dFCOVERMC(:,j) = unifrnd(dFCOVER01(j),dFCOVER99(j),[nRep,1]);
        end
    
    
        vn = {'CCx';'CGC';'CC0';'CDC';'UPVT0';'UPVTSEN';'TIC';'R2';'RMSE'};
        varFCOVER = table('Size',[nRep,numel(vn)],...
            'VariableTypes',repmat("double",[numel(vn),1]),...
            'VariableNames',vn);
        fn = varFCOVER.Properties.VariableNames;
        Yield = nan(nRep,1);
        for k = 1:nRep
            disp([num2str(i),'/',num2str(size(dFCOVER,1)),'-', num2str(k),'/',num2str(nRep)])
            FCOVERk = dFCOVERMC(k,:);
            [DTk,dFCOVERk,isNFilled] = padZeros(dT,FCOVERk,pheno,70);
            gddk = interp1(datenum(pheno.date),pheno.GDDc,datenum(DTk));   
            gddk(isnan(gddk)) = linspace(3000,50000,nnz(isnan(gddk)));
            w = ones(size(gddk));
            
            w(isNFilled)  = exp(-dFCOVERstd/0.25);
  
            
            modIk = estimateInitialOpt(gddk,dFCOVERk,pheno,GF);
            modOk = optimizeFCOVER(gddk,dFCOVERk,modIk,GF,w);
            modOk.UPVT50 = modIk.UPVT50; 
            modOk = optimizeFCOVER(gddk,dFCOVERk,modOk,GF,w); 
            for field = fn
                varFCOVER.(field{1})(k) = modOk.(field{1});            
            end
            AquaCropOS_RUN
            (modOk)
            [fink,~] = getAquacropOutputs(pathOuputAC);
            Yield(k) = fink.Yield;
        end    
        varFCOVER.Yield = Yield;

        dFCOVERMC1 = nan(nRep,numel(gdd1));
        for k = 1:nRep
            dFCOVERMC1(k,:) = funcCCI(varFCOVER(k,:),gdd1);
        end


   
    
    
        


        paramPERCi = paramPERC;    
        fnMCi = paramPERCi.Properties.VariableNames;
        for l=1:numel(fnMCi)
            paramPERCi{:,fnMCi{l}} = prctile([varFCOVER{:,fnMCi{l}}],perc,1);                    
        end
        paramPERCi.perc = perc';
        FCOVERui = prctile(dFCOVERMC1,perc,1);

        curvePERCi.perc= perc;
        curvePERCi.unc = FCOVERui;
        curvePERCi.gddc = gdd1;    
        curvePERCi.DNUM = dn1;    
        curvePERCi.DT = datetime(dn1,'ConvertFrom','datenum'); 
        
        MCi = [];
        MCi.curvePERC = curvePERCi;
        MCi.paramPERC =  paramPERCi;
        MCi.dFCOVERMC = dFCOVERMC1;
        MCi.dataFCOVERMC = dFCOVERMC;
        MCi.dataDT = dFCOVERMC;
        MC =  [MC;MCi];
    end
end
    


function mod = estimateInitialOpt(UPVTS2,FCOVERd,phenoWheat,GF)
    UPVTi = linspace(0,max(UPVTS2),1e4);
    %FCOVERd = FCOVER(150,:);

    FCOVERi = interp1(UPVTS2,FCOVERd,UPVTi);
    mod.UPVT0 = phenoWheat.GDDc(find(cellstr(phenoWheat.BBCHS)=="BBCH 10",1,'first'));
    mod.CCx =  prctile(FCOVERi,99);
    i50 = find(FCOVERi>mod.CCx/2,1,'first');
    mod.UPVT50 = UPVTi(i50);
    
    if nnz((abs(UPVTS2-mod.UPVT50)<=150))>4
        UPVTLIN = UPVTS2(abs(UPVTS2-mod.UPVT50)<=150);
        FCOVERLIN = FCOVERd(abs(UPVTS2-mod.UPVT50)<=150);
    else
        UPVTLIN = UPVTi(abs(UPVTi-mod.UPVT50)<=150);
        FCOVERLIN = FCOVERi(abs(UPVTi-mod.UPVT50)<=150);        
    end            
        
        
    lin = polyfit(UPVTLIN,FCOVERLIN,1);
    mod.CGC = lin(1);
    mod.CC0 = 0.5*mod.CCx*exp(-(mod.UPVT50-mod.UPVT0)*mod.CGC);
    if mod.CC0>0.3
        mod.CC0 =min(FCOVERd(FCOVERd>0)) ;
    end
    iSend = find(FCOVERd>0.95*mod.CCx,1,'last');
    
    mod.UPVTSEN = UPVTS2(iSend);
    FCOVERdSen = FCOVERd;
    
    FCOVERdSen(UPVTS2<mod.UPVTSEN) = [];
    CCsenInd = @(CDC) GF.CCsen(mod.CCx,CDC,mod.CCx,mod.UPVTSEN,UPVTS2(UPVTS2>=mod.UPVTSEN));
    TICsen = @(CDC) getTIC(CCsenInd(CDC),FCOVERdSen);
    TICsen(0.1)
    mod.CDC = fminbnd(TICsen,1e-4,1e-2);
    
    
    
    %mod.CC0 = 0.0675;
    
    mod.UPVT0 = phenoWheat.UPVTc(find(cellstr(phenoWheat.BBCHS)=="BBCH 10",1,'first'));
    mod.UPVTMAX = nan;
    

    
    
    
 end


 function mod1= optimizeFCOVER(UPVTS2,FCOVERd,mod0,GF,w)
%     UPVTGF =  UPVTS2(UPVTS2<UPVTsen);
%     FCOVERdF = FCOVERd(UPVTS2<UPVTsen);
    % CCGF= @(UPVT0i,CGCi,CC0i,CCxi) CCEG(CC0i,CGCi,CCxi,UPVT0i,UPVTS2) + ...
    %         CCG(CC0i,CGCi,CCxi,UPVT0i,UPVTsen,UPVTGF);

    CCGF = @(UPVT0i,CGCi,CC0i,CCxi,CDCi,UPVTseni) GF.CCI(CC0i,CGCi,CCxi,CDCi,UPVT0i,UPVTseni,UPVTS2);
    
    %mdl = fitnlm(tbl,modelfun,beta0)
    if (~isempty(w) || nargin == 4)
        TICGF = @(UPVT0i,CGCi,CC0i,CCxi,CDCi,UPVTseni) getTIC(CCGF(UPVT0i,CGCi,CC0i,CCxi,CDCi,UPVTseni),FCOVERd);
    else
        TICGF = @(UPVT0i,CGCi,CC0i,CCxi,CDCi,UPVTseni) getTIC(CCGF(UPVT0i,CGCi,CC0i,CCxi,CDCi,UPVTseni),FCOVERd,w);
    end
    
    TICGFx = @(x) TICGF(x(1),x(2),x(3),x(4),x(5),x(6));
    
    x0 = [mod0.UPVT0;mod0.CGC;mod0.CC0;mod0.CCx;mod0.CDC;mod0.UPVTSEN];
    lowerBound = [0;0.25*mod0.CGC;3e-2;prctile(FCOVERd,95);0.005*mod0.CDC;mod0.UPVTSEN];
    upperBound=[0.5*mod0.UPVT50;5*mod0.CGC;0.25*mod0.CCx;1;5e-3;1.5*mod0.UPVTSEN];
    
    opts = optimoptions('fmincon','Display','off');    
    x = fmincon(TICGFx,x0,[],[],[],[],lowerBound,upperBound,[],opts);
    


    
    
    
    %CCGF1 = @(x1,UPVT1) GF.CCI(x1(1),x1(2),x1(3),x1(4),x1(5),x1(6),UPVT1);

    %mdl = fitnlm(UPVTS2,FCOVERd,CCGF1,x0)
    
    mod1.UPVT0 = x(1);
    mod1.CGC = x(2);
    mod1.CC0 = x(3);
    mod1.CCx = x(4);
    mod1.CDC = x(5);
    mod1.UPVTSEN = x(6);
    mod1.UPVT50 = nan;
    mod1.UPVTMAX = nan;
    mod1.TIC = TICGFx(x);

    R2GF = @(UPVT0i,CGCi,CC0i,CCxi,CDCi,UPVTseni) getR2(CCGF(UPVT0i,CGCi,CC0i,CCxi,CDCi,UPVTseni),FCOVERd);
    R2GFx = @(x) R2GF(x(1),x(2),x(3),x(4),x(5),x(6));    
    RMSEGF = @(UPVT0i,CGCi,CC0i,CCxi,CDCi,UPVTseni) getRMSE(CCGF(UPVT0i,CGCi,CC0i,CCxi,CDCi,UPVTseni),FCOVERd);
    RMSEGFx = @(x) RMSEGF(x(1),x(2),x(3),x(4),x(5),x(6));       
    mod1.R2 = R2GFx(x);
    mod1.RMSE = RMSEGFx(x);

end





function TIC = getTIC(OMod,OMeas,w)
    if nargin ==3
        iocN = calcIoC(OMod,OMeas,w);
    else
        iocN = calcIoC(OMod,OMeas,ones(size(OMeas)));

    end
    TIC = iocN.TIC;
end


function RMSE = getRMSE(OMod,OMeas,w)
    if nargin ==3
        iocN = calcIoC(OMod,OMeas,w);
    else
        iocN = calcIoC(OMod,OMeas,ones(size(OMeas)));

    end
    RMSE = iocN.RMSE;

end

function R2 = getR2(OMod,OMeas,w)
    if nargin ==3
        iocN = calcIoC(OMod,OMeas,w);
    else
        iocN = calcIoC(OMod,OMeas,ones(size(OMeas)));

    end
    R2 = iocN.R2;

end


function [fin,dyn] = getAquacropOutputs(pathOuput)
   CG = readtable([pathOuput,'Sample_CropGrowth.txt']);
   WG = readtable([pathOuput,'Sample_WaterFluxes.txt']);      
   WC = readtable([pathOuput,'Sample_WaterContents.txt']);
   WC.Properties.VariableNames(6:end) = {'WC_Ap','WC_E','WC_Bt'};
   dT = datetime(CG.Year,CG.Month,CG.Day);
   dT = table(dT,'VariableNames',{'Date'});
   dyn = [dT,CG(:,6:end), WG(:,6:end),WC(:,6:end)];
   dyn = standardizeMissing(dyn,-999);
   fin = importFinalOutput([pathOuput,'Sample_FinalOutput.txt']);
   fin = fin(:,5:end);
  

end

function  SampleFinalOutput = importFinalOutput(pathFinalOutput)
    opts = fixedWidthImportOptions("NumVariables", 8);

    % Specify range and variable widths
    opts.DataLines = [2, Inf];
    opts.VariableWidths = [6, 12, 14, 17, 19, 18, 11, 19];

    % Specify column names and types
    opts.VariableNames = ["Season", "CropType", "PlantDate", "PlantSimDate", "HarvestDate", "HarvestSimDate", "Yield", "TotIrr"];
    opts.VariableTypes = ["double", "string", "datetime", "double", "datetime", "double", "double", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, "CropType", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "CropType", "EmptyFieldRule", "auto");
    opts = setvaropts(opts, "PlantDate", "InputFormat", "dd/MM/yyyy");
    opts = setvaropts(opts, "HarvestDate", "InputFormat", "dd/MM/yyyy");

    % Import the data
    SampleFinalOutput = readtable(pathFinalOutput, opts);

end

