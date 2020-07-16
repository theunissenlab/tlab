res=[];
i=0;
for alpha0=-1.0:0.005:.0;
    for beta0=-1.0:0.005:1.0;
        ttime=linspace(0,0.15,10);
        alpha=alpha0+0*ttime;
        beta=beta0+0*ttime;        
        f=findFun(smBGAs(0.15,24000,alpha',ttime',beta',ttime'));
        
        res=[res; alpha0 beta0 f];
    %    plot3(res(:,1),res(:,2),res(:,3),'+');
     %   drawnow
     i=i+1;
     if i>1000;
         i=0;
         [i alpha0 beta0 f]
     end;
    end   
end
