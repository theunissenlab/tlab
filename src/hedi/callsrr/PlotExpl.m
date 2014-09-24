function [dat1,dat2]=PlotExpl(expl)
ns=length(expl);
colo=[ 1 0 0; 0 1 0 ; 0 0 1 ;0 0 0;1 0 1; 0 1 1];
close all
titles={'2m: Spec2m/Moto2m','16m: Spec16m/Moto2m', '16m: Spec2m/Moto2m','16m: Spec16m/Spec2m'};
for i=1:4;
 subplot(2,2,i)
 t = 0 : .01 : 2 * pi;
 P = polar(t, 2 * ones(size(t)));
 set(P, 'Visible', 'off')
 title(titles(i));
end;
for site=1:ns;
    ne=length(expl{site});
    offset=(site-1)*16;
    dat1=[];
    dat2=[];
    dat3=[];
    dat4=[];
    dat5=[];
   
    for el=1:ne;
        d1=expl{site}{el}(1);
        d2=expl{site}{el}(2);
        d3=expl{site}{el}(3);
        d4=expl{site}{el}(4);
        d5=expl{site}{el}(5);
        dat1=[dat1;d1];
        dat2=[dat2;d2];
        dat3=[dat3;d3];
        dat4=[dat4;d4];
        dat5=[dat5;d5];
    end;
    
    if length(dat5)>0;
    subplot(2,2,1)    
    hold on
    h=polar([2*pi*((1:ne))/ne,2*pi/ne],([dat1./dat2; dat1(1)./dat2(1)])');
    set(h,'LineWidth',2.0,'Color',colo(site,:));
    %axis([0 2 0 360]);
    
    subplot(2,2,2)
    hold on
    h=polar([2*pi*((1:ne))/ne,2*pi/ne],([dat3./dat4; dat3(1)./dat4(1)])');%,'LineWidth',2.0,'Color',colo(site,:));
    set(h,'LineWidth',2.0,'Color',colo(site,:));
    
    subplot(2,2,3)
    hold on
    h=polar([2*pi*((1:ne))/ne,2*pi/ne],([dat5./dat4; dat5(1)./dat4(1)])');%,'LineWidth',2.0,'Color',colo(site,:));
     set(h,'LineWidth',2.0,'Color',colo(site,:));
    
    subplot(2,2,4)    
    hold on
    h=polar([2*pi*((1:ne))/ne,2*pi/ne],([dat3./dat5; dat3(1)./dat5(1)])');%,'LineWidth',2.0,'Color',colo(site,:));
     set(h,'LineWidth',2.0,'Color',colo(site,:));
    
    end
end;
% hold off
% close all
% bar(5*dat1(:,1),dat1(:,2),'r')
% hold on
% bar(5*dat2(:,1)+2.5,dat2(:,2),'k');
        
% nsp=length(find(dat1(:,2)>dat2(:,2)))
% ne=length(find(dat1(:,2)<dat2(:,2)))
% ne+nsp
% 
% nsp16=length(find(dat5(:,2)>dat3(:,2)))
% ne16=length(find(dat5(:,2)<dat3(:,2)))
% ne16+nsp16
% 
% nsp16v2=length(find(dat3(:,2)>dat4(:,2)))
% nm16v2=length(find(dat3(:,2)<dat4(:,2)))
% nsp16v2+nm16v2
% 
% 
% 
% subplot(2,2,1)
% bar(dat1(:,1),dat1(:,2)-dat2(:,2))
% title('spec 2m vs moto 2m at 2m');
% subplot(2,2,2)
% bar(dat3(:,1),dat3(:,2)-dat4(:,2))
% title('spec 16m vs moto 2m at 16m');
% subplot(2,2,3)
% bar(dat1(:,1),dat3(:,2)-dat5(:,2))
% title('spec 2m vs moto at 16m')
% subplot(2,2,4);
% plot(dat3(:,2)-dat5(:,2),dat3(:,2)-dat4(:,2),'+',-0.5:0.1:0.5,0*(-0.5:0.1:0.5),'k',0*(-0.5:0.1:0.5),-0.5:0.1:0.5,'k')