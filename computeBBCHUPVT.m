function phenoWheat = computeBBCHUPVT(BBCHTAB,clim,dateS,dateE)
  DNUMS = datenum(dateS,'yyyy-mm-dd');
  DNUME = datenum(dateE,'yyyy-mm-dd');

  TCULT = 0.5*(clim.TMIN + clim.TMAX);
  %TCULT = clim.TMEAN;
  TOD = clim.TOD;
  iStart = find(clim.DNUM<DNUMS,1,'last');
  iEnd = find(clim.DNUM>DNUME,1,'first');
  %clim.date(iStart)
  UPVT = calcUPVT(TOD,TCULT,iStart,iEnd);
  
  UPVT(1:iStart) = 0;
  UPVT(iEnd:end) = 0;  
  
  UPVTc = cumsum(UPVT);
  UPVTc = UPVTc - max(UPVT(clim.DNUM<DNUMS));
  
  UPVTc(1:iStart) = nan;
  UPVTc(iEnd:end) = nan;  
  
  BBCHS = calcBBCH(BBCHTAB,UPVTc);
  iFloraison = find(BBCHS=="BBCH 65",1,'first');
  UPVT = calcUPVT(TOD,TCULT,iStart,iFloraison);
  UPVT(1:iStart) = 0;
  UPVT(iEnd:end) = 0;    
  UPVTc = cumsum(UPVT);
  UPVTc = UPVTc - max(UPVT(clim.DNUM<DNUMS));
  
  UPVTc(1:iStart) = nan;
  UPVTc(iEnd:end) = nan;   
  
  GDD = TCULT;
  GDD(1:iStart) = 0;
  GDD(iEnd:end) = 0;     
  GDDc = cumsum(GDD);
  GDDc = GDDc - max(GDD(clim.DNUM<DNUMS));
  GDDc(1:iStart) = nan;
  GDDc(iEnd:end) = nan;  
  BBCHS = calcBBCH(BBCHTAB,UPVTc);
  
  
  results0.DNUM = clim.DNUM(~isnan(UPVTc));
  results0.date = clim.date(~isnan(UPVTc));
  results0.UPVT = UPVT(~isnan(UPVTc));
  results0.UPVTc = UPVTc(~isnan(UPVTc));
  results0.GDD = GDD(~isnan(UPVTc));
  results0.GDDc = GDDc(~isnan(UPVTc));

  matBBCH = cell2mat(BBCHS);
  matBBCH(results0.DNUM<DNUMS,:) = repmat("NOT DEF",[nnz(results0.DNUM<DNUMS),1]); 
  matBBCH(results0.DNUM>DNUME,:) = repmat("NOT DEF",[nnz(results0.DNUM>DNUME),1]); 
  
  results0.BBCHS = matBBCH(~isnan(UPVTc),:);  
  
  results.date = (results0.date(1):results0.date(end))';
  results.DNUM = datenum(results.date);
  results.UPVT = interp1(results0.DNUM,results0.UPVT,results.DNUM,'spline');
  results.UPVTc = interp1(results0.DNUM,results0.UPVTc,results.DNUM,'spline');
  results.GDD = interp1(results0.DNUM,results0.GDD,results.DNUM,'spline');
  results.GDDc = interp1(results0.DNUM,results0.GDDc,results.DNUM,'spline');
  
  [BBCHpW,~,ipW] = unique(string([results0.BBCHS]));

  results.BBCHS = BBCHpW(interp1(results0.DNUM,ipW,results.DNUM,'nearest'));

  results

  phenoWheat = struct2table(results);

  
  
end




function UPVT = calcUPVT(TOD,TCULT,iStart,iFloraison)
  % Photoperiod effect
  
  SENSIPHOT = 0;
  PHOBASE = 6;
  PHOSAT = 20;  
  RFPI = 1 - (1-SENSIPHOT)*(PHOSAT-TOD)/(PHOSAT-PHOBASE) ;
  RFPI =  max(min(1,RFPI),0);
  % Vernalisation
  TFROID = 6.5;
  AMPFROID = 10;
  JVC_MINI = 5;
  JVI = min(max((1 - ((TFROID-TCULT)/AMPFROID).^2),0),1);
  JVI(1:iStart) = 0;
  JVC_V = 50;
  RFVI = min(max((cumsum(JVI)-JVC_MINI)./((JVC_V-JVC_MINI)),0),1);
  % Limitation in temperature
  TDMIN = 0;
  TDMAX = 28;
  TCXSTOP = 45;
  UDEV = zeros(size(TCULT));

  UDEV((TCULT >= TDMIN) & (TCULT <= TDMAX)) = TCULT((TCULT >= TDMIN) & (TCULT <= TDMAX)) -TDMIN;
  UDEV((TCULT >= TDMAX) & (TCULT <= TCXSTOP)) =  (TCXSTOP-TCULT((TCULT >= TDMAX) & (TCULT <= TCXSTOP)))*(TDMAX-TDMIN)./(TCXSTOP-TDMAX);
  UDEV(1:iStart) = 0;
  RFPI(1:iStart) = 0;
  RFPI(iFloraison:end) = 1;
  % Computation of UPVT
  UPVT = UDEV.*RFPI.*RFVI;
%   plot(UPVT,'b')
%   hold on
%   plot(UDEV,'r')


end



function BBCHS = calcBBCH(BBCH,DJc)
    BBCHS = repmat("BBCH XX",[length(DJc),1]);
    nHS = length(BBCH.BBCH_STAGE);
    for i = 1:nHS
        if(i < nHS)
            Ii = (DJc - BBCH.UPVTC(i))/(BBCH.UPVTC(i+1) - BBCH.UPVTC(i));
            BBCHS(Ii <= 1.0 & Ii >= 0) = BBCH.BBCH_STAGE{i};
        else
            Ii = DJc/BBCH.UPVTC(i);
            BBCHS(Ii >= 1.0) = BBCH.BBCH_STAGE{i};            
        end
    end
    
end