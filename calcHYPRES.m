function mod = calcHYPRES(C,S,OM,D,topsoil)
    % Input
    % C :  pourcentage d'argile
    % S : pourcentage de limon
    % OM : pourcentage de matiere organique
    % D : densite apparente
    % topsoil : 1 si horizon de surface, 0 autrement...
    
    % Output : 
    % thetaS : teneur a saturation
    % alphas : 
    % ns
    % ls
    % Kss
    
    S(S<=0) = nan;
    C(C<=0) = nan;
    OM(OM<=0) = nan;    
    
    thetaS = 0.7919+0.001691.*C-0.29619.*D-0.000001491.*S.^2+0.0000821.*OM.^(2)+0.02427.*C.^(-1)+0.01113.*S.^(-1) ...
         +0.01472.*reallog(S)-0.0000733.*OM.*-0.000619.*D.*C-0.001183.*D.*OM-0.0001664.*topsoil.*S;

     alphas = -14.96+0.03135.*C+0.0351.*S+0.646.*OM+15.29.*D-0.192.*topsoil...
        -4.671.*D.^(2)-0.000781.*C.^(2)-0.00687.*OM.^2+0.0449.*OM.^(-1)+0.0663.*reallog(S)...
        +0.1482.*reallog(OM)-0.04546.*D.*S-0.4852.*D.*OM+0.00673.*topsoil.*C;

    ns = -25.23-0.02195.*C+0.0074.*S-0.1940.*OM+45.5.*D-7.24.*D.^2+0.0003658.*C.^2 ... 
        +0.002885.*OM.^2-12.81.*D.^(-1)-0.1524.*S.^(-1)-0.01958.*OM.^(-1)-0.2876.*reallog(S) ...
        -0.0709.*reallog(OM)-44.6.*reallog(D)-0.02264.*D.*C+0.0896.*D.*OM+0.00718.*topsoil.*C;

    ls = 0.0202+0.0006193.*C.^2-0.001136.*OM.^2-0.2316.*reallog(OM)-0.03544.*D.*C ...
        +0.00283.*D.*S+0.0488.*D.*OM;


    Kss = 7.755+0.0352.*S+0.93.*topsoil-0.967.*D.^(2)-0.000484.*C.^2-0.000322.*S.^2 ...
        +0.001.*S.^(-1)-0.0748.*OM.^(-1)-0.643.*reallog(S)-0.01398.*D.*C-0.1673.*D.*OM ...
        +0.02986.*topsoil.*C-0.03305.*topsoil.*S;
    
    mod.Ks = (exp(Kss));
    mod.alpha = (exp(alphas));
    mod.n = 1 + (exp(ns));
    mod.l = (10.*tanh(0.5.*ls)); 
    mod.thetaR = zeros(size(mod.Ks));
    mod.thetaS = (thetaS);
    
end