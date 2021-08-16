function addNorth(go,Location,Height,ax)
    if nargin==3
        ax = gca;
    end
    xyh1= [0,1.2;sqrt(0.25),0;0,1-sqrt(0.5);0,1.2];
    xyh1= Height*xyh1/max(xyh1(:,2));

    xyh2= xyh1;
    xyh2(:,1) = -xyh2(:,1);
    W=2*max(xyh1(:,1));

    k=0.05;
    ext = go.getextent;
    if(contains(Location,'South'))
        y0 = ext(3) + k*(ext(4)-ext(3));

    elseif(contains(Location,'North'))

        y0 = ext(4) - k*(ext(4)-ext(3)) - 1.25*Height;
    else
        y0 = 0.5*(ext(3)+ext(4));
    end

    if(contains(Location,'East'))
        x0 = ext(2) - k*(ext(2)-ext(1))- 1.1*W;
    elseif(contains(Location,'West'))
        x0 = ext(1) + k*(ext(2)-ext(1))+ 1.1*W;
    else
        x0 = 0.5*(ext(2)+ext(1));

    end



    hold(ax,'on')

    fill(ax,x0+xyh1(:,1),y0+xyh1(:,2),'k','EdgeColor','k','LineWidth',1.5)
    fill(ax,x0+xyh2(:,1),y0+xyh2(:,2),[1,1,1],'EdgeColor','k','LineWidth',1.5)
    hold(ax,'off')  
    
end