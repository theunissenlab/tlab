function testbuttondown(src, evnt)

global cellselect
global selectcolors

pt=get(src,'CurrentPoint');

dat=get(src,'UserData');

x=pt(1,1);
y=pt(1,2);


thisxy=[NaN NaN NaN NaN NaN];

if ~isempty(dat.xy)
    thisxy=dat.xy(end,:);
    if ~isnan(thisxy(3))
        thisxy=[NaN NaN NaN NaN NaN];
    end
end

if isnan(thisxy(1))
    thisxy = [x y NaN NaN cellselect];
    dat.xy(end+1,:)=thisxy;
    set(src,'UserData',dat);
else
    thiscolor = selectcolors(cellselect+1,:);
    thisxy(3:4) = [x y];
    dat.xy(end,:)=thisxy;
    f=zeros(1,size(dat.snips,2));
    s=size(dat.snips,1);

    for ii=1:size(dat.snips,2)
        [xi yi]=polyxpoly(1:s,dat.snips(:,ii),thisxy([1 3]),thisxy([2 4]));
        if ~isempty(xi)
            f(ii)=1;
        end
    end

    
    dat.sortcode(logical(f))=cellselect;
    plot([thisxy(1) x], [thisxy(2) y],'color',thiscolor,'linewidth',2);
    set(dat.h(find(f)),'color',thiscolor);
    
    set(src,'UserData',dat);
    
end


function xout = ppwrap(s,t,xy1,xy2)

    [xout, ~] = polyxpoly(s,t,xy1,xy2);
    if isempty(xout)
        xout=0;
    end

    