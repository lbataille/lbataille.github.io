function ioc = calcIoC(OMod,OMeas,w)
    
    if nargin == 3
        ioc.MSE = mean(w(:).*(OMod(:)-OMeas(:)).^2,'Omitnan');
    else
        ioc = calcIoC(OMod,OMeas,ones(size(OMeas)));
    end
% size(OMod)
% size(OMeas)
    
    ioc.RMSE = sqrt(ioc.MSE);
    ioc.rRMSE = ioc.RMSE/sqrt(mean(OMeas(:).^2,'Omitnan'));
    ioc.TIC = ioc.RMSE/(sqrt(mean(OMeas(:).^2,'Omitnan'))+sqrt(mean(OMod(:).^2,'Omitnan')));
    ioc.R2 = (corr(OMeas(:),OMod(:), 'rows','complete')).^2;
%     N = numel(isfinite(OMod-OMeas));
%     p = ?
%     ioc.AIC = N*sum(OMod-OMeas,'Omitnan') + 2*p;
end