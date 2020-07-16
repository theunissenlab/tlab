function mysuperplot(Score)
n=16;
hold on 
col=[[1.,0,0];[0 1 0];[0 0 1];[0 0 0];[1 1 0];
    [1 0 1];[0 1 1];[0.5 0.5 0.5];[0.8 0.2 0.5];
    [0.2 0.8 0.5];[0.2 0.5 0.8];[0.8 0.8 0.1];[0.1 0.8 0.8];[0.8 0.1 0.8];[0.2 0.2 0.5];[1 0.5 0.5]];

for i=1:n;
    s=(i-1)*16+1;
    e=s+15;
    [256*i/16 0 256*i/16]
    z=colormap;
    z(10)
    myS=[mean(Score(s:e,1)),mean(Score(s:e,2)),mean(Score(s:e,3))];
    plot3(Score(s:e,1),Score(s:e,2),Score(s:e,3),'.',myS(1),myS(2),myS(3),'+','Color',col(i,:),'MarkerSize',20);
    str=sprintf('%d',i);
    text(myS(1),myS(2),myS(3),str,'FontSize',20)
end;