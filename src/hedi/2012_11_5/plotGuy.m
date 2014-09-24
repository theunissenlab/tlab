function plotGuy(birds,sex,id,call);


os=birds(sex).birds(id).calls(call).o_spec;
fs=birds(sex).birds(id).calls(call).f_spec;
syr=birds(sex).birds(id).calls(call).syr;
fil=birds(sex).birds(id).calls(call).fil;
msex=['m','f'];
ttname=sprintf('%cbird%dcall%d',msex(sex),id,call);

fname=sprintf('../output/image/im%cbird%dcall%d.pdf',msex(sex),id,call);


subplot(2,2,1)
imagesc(os.t-os.t(1),os.f,os.spec);
axis xy
title('Original Spectrum')
subplot(2,2,2)
imagesc(fs.t,fs.f,fs.spec);
axis xy
title('Fitted Spectrum');
subplot(2,2,3)
plot(1:length(syr),syr);
title('Syrinx Contraction');
subplot(2,2,4)
plot(1:length(fil),fil);
title('Mean Freq');

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

text(0.5, 1,ttname,'HorizontalAlignment','center','VerticalAlignment', 'top')

%print('-dpdf',fname);