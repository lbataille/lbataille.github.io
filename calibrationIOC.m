function IOC = calibrationIOC(iocName,tMeas,OMeas,funcODE,varargin)
   arginODE = varargin{1};
   [OMod0,tMod0] = funcODE(arginODE{:}); 
   OMod = interp1(tMod0,OMod0,tMeas);

   iocNum = calcIoC(OMod,OMeas);
%    switch iocName
% 		case 'AIK'
% 			namesP = getArgumentsFunc(funcODE)
% 			for i = 1:numel(namesP)
% 				IOC.AIC.(namesP{i}) = numel(tMes)*(Omod - Omes) + 2*arginODE{i};
% 			end
% 		otherwise
			IOC = iocNum.(iocName)
   %end
end