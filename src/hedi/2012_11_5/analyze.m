% get all vocal stuff 
close all;

col=[[1.,0,0];[0 1 0];[0 0 1];[0 0 0];[1 1 0];
    [1 0 1];[0 1 1];[0.5 0.5 0.5];[0.8 0.2 0.5];
    [0.2 0.8 0.5];[0.2 0.5 0.8];[0.8 0.8 0.1];[0.1 0.8 0.8];[0.8 0.1 0.8];[0.2 0.2 0.5];[1 0.5 0.5]];
data_dir='../output';
for bird=[1:11,13:16];
    ldir=sprintf('%s/res.*v4*fbird%dcall*.mat',data_dir,bird);
    d=dir(ldir);
    n=length(d);
%    c=sprintf('%c',col(bird));
%    cc=sprintf('%c+',col(bird));
    max_len=-1;
    for i=1:n;
        fname=sprintf('%s/%s',data_dir,d(i).name);
        load(fname)
        if max_len<length(res.data(:,1));
            max_len=length(res.data(:,1));
            avg_res=res.data(:,1);
        end;
    end;
    avg_data=zeros(max_len,1);
    avg_filter=zeros(max_len,1);
    for i=1:n;
        fname=sprintf('%s/%s',data_dir,d(i).name);
        load(fname)       
        hold on;
     %   subplot(2,1,1)                 
          y=lowP(res.data(:,2),res.data(:,1),5);
          if length(y)<max_len;
              y=[y;,y(end)*ones(max_len-length(y),1)];
          end;
          avg_data = avg_data + y;
        fi=getFilter(res.data(:,1),res.data(:,3:end));        
        yf=lowP(fi(:,1),res.data(:,1),20);
        if length(yf)<max_len;
              yf=[yf;,yf(end)*ones(max_len-length(yf),1)];
          end;
          avg_filter = avg_filter + yf;
  %      plot3(res.data(:,1),y,yf,c);
  
%         
%         subplot(2,1,2)
%          hold on;
%         
%         fi=getFilter(res.data(:,1),res.data(:,3:end));
%         %errorbar(fi(:,1),fi(:,2),fi(:,2)-fi(:,1),c);
%         yf=lowP(fi(:,2),fi(:,1),5);
%         sdf=lowP(fi(:,3),fi(:,1),5);
%        
%         plot(fi(:,1),yf,c,fi(:,1),sdf,cc);
%         drawnow
    end;
    bird
    col(bird,:)
    plot3(avg_data/n,avg_filter/n,avg_res,'LineWidth',2.0,'Color',col(bird,:));
      %  axis([0 300 -0.50 -0.35])               
        drawnow
end;
    