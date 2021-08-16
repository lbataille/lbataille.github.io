function [RGBu8,maskAlpha] = collectRGBSentinel(S2DAT,DTi)

    dat2path = @(s2di) strcat(s2di.folder,"/",s2di.name);
    B02 = GRIDobj(dat2path(S2DAT([S2DAT.DTJ]==DTi & ([S2DAT.BAND]=="B02"))));
    B03 = GRIDobj(dat2path(S2DAT([S2DAT.DTJ]==DTi & ([S2DAT.BAND]=="B03"))));
    B04 = GRIDobj(dat2path(S2DAT([S2DAT.DTJ]==DTi & ([S2DAT.BAND]=="B04"))));
    RGB = zeros([size(B02.Z),3]);
    RGB(:,:,1) = B02.Z;
    RGB(:,:,2) = B03.Z;
    RGB(:,:,3) = B04.Z;
    prx = prctile(RGB,[5,95],[1,2]);
    repmat(prx,[size(RGB,1:2),1])
    RGBp5 = repmat(prx(1),[size(RGB,1:2),1]);
    RGBp95 = repmat(prx(2),[size(RGB,1:2),1]);
    RGBu8 = uint8(255*(RGB-RGBp5)./(RGBp95-RGBp5));
    maskAlpha = B02.Z;
end