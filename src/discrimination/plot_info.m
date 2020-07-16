% Plot info data
clear all;
load fieldL/zfnormal_fieldL

% Fix some of the data format
number_cells = length(all_outputs);

for is=1:nstim
    fprintf(1,'Stimulus %s\n', stim_type{is});

    for nc=1:number_cells
        output_db(nc) = all_outputs{nc, is+1};           % Remove the +1... 
    end

    % Measures of mutual information
    gamma_mi = vertcat(output_db(:).gamma_mutual_info);
    rate_mi_biased = vertcat(output_db(:).rate_info_biased);
    rate_mi = vertcat(output_db(:).rate_info_bcorr);
    confusion_mi = vertcat(output_db(:).mi_confusionB);            % The B version is with all other trials as templates
    zdist_mi = vertcat(output_db(:).mizdistSB);
    zdist_ncmi = vertcat(output_db(:).mizdistncSB);

    % Measures of distances
    dprime = vertcat(output_db(:).avgdprime);
    zdist = vertcat(output_db(:).zdistSB);

    % Measures of percent corred
    percorrect = vertcat(output_db(:).percorrectB);
    pzdist = vertcat(output_db(:).pzdistSB);

    winSize = output_db(1).winsize;

    %Make scatter plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1);
    % Rate Information vs. Gamma Information at 4 different time window
    
    plotwinid = [3 5 6 7];
    for ip=1:4
        subplot(2,2,ip);
        
        plot(gamma_mi, rate_mi(:,plotwinid(ip)), '+');
        axis([0 20 0 30]);
        axis_auto = axis();
        hold on;
        plot([0 axis_auto(2)], [0 axis_auto(2)], 'k');
        [b,bint,r,rint,stats] = regress(rate_mi(:,plotwinid(ip)),[gamma_mi ones(length(gamma_mi),1)]);
        plot(gamma_mi, gamma_mi*b(1) + b(2),'k:');
        hold off;
        title(sprintf('Slope %.2f R2 %.2f P %.4f',b(1), stats(1), stats(3)));
        xlabel('Gamma Inf (bits/s)');
        ylabel(sprintf('Rate Inf %d ms (bits/s)', winSize(plotwinid(ip))));
    end

    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    figure(2);
    % Confusion Information vs. Gamma Information
    subplot(1,1,1);
    bestwinip = 3;
    plot(gamma_mi, confusion_mi(:,bestwinip), '+');    
    axis([0 20 0 20]);
    axis_auto = axis();
    hold on
    plot([0 axis_auto(2)], [0 axis_auto(2)], 'k');
    plot([0 axis_auto(2)], [log2(20) log2(20)], 'k--');
    [b,bint,r,rint,stats] = regress(confusion_mi(:,bestwinip),[gamma_mi ones(length(gamma_mi),1)]);
    plot(gamma_mi, gamma_mi*b(1) + b(2),'k:');
    hold off;
    title(sprintf('Slope %.2f R2 %.2f P %.4f',b(1), stats(1), stats(3)));
    xlabel('Gamma MI (bits/s)');
    ylabel(sprintf('Confusion MI Tau = %d ms (bits/s)', winSize(bestwinip)));


     figure(3);
    % Z-dist Information vs. Gamma Information 
    subplot(1,1,1);
    bestwinip = 3;
    plot(gamma_mi, zdist_mi(:,bestwinip), '+');
    axis([0 20 0 20]);
    axis_auto = axis();
    hold on;
    plot([0 axis_auto(2)], [0 axis_auto(2)], 'k');
    [b,bint,r,rint,stats] = regress(zdist_mi(:,bestwinip),[gamma_mi ones(length(gamma_mi),1)]);
    plot(gamma_mi, gamma_mi*b(1) + b(2),'k:');
    hold off;
    title(sprintf('Slope %.2f R2 %.2f P %.4f',b(1), stats(1), stats(3)));
    xlabel('Gamma MI (bits/s)');
    ylabel(sprintf('Z MI Anthropic  Tau = %d ms (bits/s)',winSize(bestwinip)));

  
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(4);
    % Rate Information vs. Rate Information (without Bias Correction)
    plot(rate_mi(:,2), rate_mi_biased(:,2), '+');
    axis_auto = axis();
    hold on;
    plot([0 axis_auto(2)], [0 axis_auto(2)], 'k');
    [b,bint,r,rint,stats] = regress(rate_mi_biased(:,2),[rate_mi(:,2) ones(length(gamma_mi),1)]);
    plot(rate_mi(:,2), rate_mi(:,2).*b(1) + b(2),'k:');
    hold off;
    title(sprintf('Bias Corrections: Slope %.2f R2 %.2f P %.4f', b(1), stats(1), stats(3)));
    xlabel('Rate Information Corrected (bits/s)');
    ylabel('Rate Information Biased (bits/s)');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    figure(5);
    % Rate Information vs. Window Size
    errorbar(winSize, mean(rate_mi,1), 2*std(rate_mi)./sqrt(number_cells));   
    hold on;
    xlabel('Window (ms)');
    ylabel('Average Rate Information (bit/s)');
    title('Effect of window on Rate Information');
    axis_auto = axis();
    plot([0 axis_auto(2)],[mean(gamma_mi) mean(gamma_mi)], 'k');
    hold off;
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(6);
    % Confusion Information vs. Window Size
    errorbar(winSize, mean(confusion_mi,1), 2*std(confusion_mi)./sqrt(number_cells));   
    xlabel('Window (ms)');
    ylabel('Average Confusion Information (bit/s)');
    title('Effect of window on Confusion Information');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(7);
    % Z-dist Information vs. Window Size
    errorbar(winSize, mean(zdist_mi,1), 2*std(zdist_mi)./sqrt(number_cells));   
    xlabel('Window (ms)');
    ylabel('Average Z-dist Information (bit/s)');
    title('Effect of window on Z-dist Information');
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   figure(8);
    % Percent Correct vs. Window Size
    errorbar(winSize, 100*mean(percorrect,1), 200*std(percorrect)./sqrt(number_cells));
    xlabel('Tau Exponential Decay (ms)');
    ylabel('Average Percent Correct (%)');
    title('Effect of Tau on Percent Correct');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   figure(9);
    % Percent Correct Z-dist vs. Window Size
    errorbar(winSize, 100*mean(pzdist,1), 2*std(pzdist)./sqrt(number_cells));
    xlabel('Tau Exponential Decay (ms)');
    ylabel('Average Percent Correct (%)');
    title('Effect of Tau on P Correct Z-dist');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   figure(10);
    % Percent Correct Z-dist vs. Percent Correct
    bestwinip = 3;
    plot(percorrect(:,bestwinip), pzdist(:,bestwinip), '+');
    xlabel(sprintf('Percent Correct (%d ms)', winSize(bestwinip)));
    ylabel(sprintf('Z-dist Correct (%d ms)', winSize(bestwinip)));
    axis_auto = axis();
    hold on;
    plot([0 axis_auto(2)],[0 axis_auto(2)], 'k');
    [b,bint,r,rint,stats] = regress(pzdist(:,bestwinip),[percorrect(:,bestwinip) ones(length(gamma_mi),1)]);
    plot(percorrect(:,1), percorrect(:,1)*b(1) + b(2),'k:');
    hold off;
    title('Z-dist Correct vs Percent Correct');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    figure(11);
    % Percent Correct vs. Gamma Information
     bestwinip = 3;
    plot(gamma_mi, percorrect(:,bestwinip), '+');
    axis([0 20 0 1]);
    ylabel('P Correct');
    xlabel('Gamma MI (bits/s)');
    hold on;
    [b,bint,r,rint,stats] = regress(percorrect(:,bestwinip), [gamma_mi ones(length(gamma_mi),1)]);
    plot(gamma_mi, gamma_mi*b(1) + b(2),'k:');  
    plot([0 20],[1/20 1/20], 'k--');   % Chance level
    title(sprintf('Percent Correct %d ms vs Gamma MI R2=%.2f', winSize(bestwinip), stats(1)));
    hold off;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     figure(12);
    % Percent Z-dist Correct vs. Gamma Information
   bestwinip = 3;
    plot(gamma_mi, pzdist(:,bestwinip), '+');
    title(sprintf('Percent Correct Z-Dist(Tau= %d ms) vs Gamma Info', winSize(bestwinip)));
    xlabel('Gamma Information');
    ylabel('PC Z-dist');
%     hold on
%     [B,BINT,R,RINT,STATS] = regress(info,[percorrect(:,1) ones(length(percorrect),1)]);
%     plot(percorrect(:,1), percorrect(:,1)*B(1) + B(2),'k');
%     hold off;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     figure(13);
    % D-prime vs. Gamma Information

    plot(gamma_mi, dprime, '+');
    xlabel('Gamma Information');
    ylabel('D Prime');
    hold on;
    [b,bint,r,rint,stats] = regress(dprime,[gamma_mi ones(length(gamma_mi),1)]);
    plot(gamma_mi, gamma_mi*b(1) + b(2),'k:');
    axis([0 20 0 2]);
    hold off;
    title(sprintf('Slope %.2f R2 %.2f P %.4f',b(1), stats(1), stats(3)));
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     figure(14);
    % Z-dist vs. Gamma Information
   bestwinip = 3;
    plot(gamma_mi, zdist(:,bestwinip), '+');
    title('Z dist vs Gamma Info');
    xlabel('Gamma Information');
    ylabel(sprintf('Z Dist (%d ms)', winSize(bestwinip)));
   
%     hold on
%     [B,BINT,R,RINT,STATS] = regress(info,[percorrect(:,1) ones(length(percorrect),1)]);
%     plot(percorrect(:,1), percorrect(:,1)*B(1) + B(2),'k');
%     hold off; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     figure(15);
    % Z-dist vs. Dprime 
   bestwinip = 3;
    plot(dprime, zdist(:,bestwinip), '+');
    title('D prime vs zdist');
    xlabel('D Prime');
    ylabel(sprintf('Z Dist (%d ms)', winSize(bestwinip)));
   
%     hold on
%     [B,BINT,R,RINT,STATS] = regress(info,[percorrect(:,bestwinip) ones(length(percorrect),1)]);
%     plot(percorrect(:,1), percorrect(:,1)*B(1) + B(2),'k');
%     hold off; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     figure(16);
    % Confusion Information vs. Percent Correct 
   bestwinip = 3;
    plot(percorrect(:,bestwinip), confusion_mi(:,bestwinip), '+');
    axis([0 1 0 5]);
    xlabel('P Correct');
    ylabel('Confusion MI (bits/s)');
    hold on;
    [b,bint,r,rint,stats] = regress(confusion_mi(:,bestwinip), [percorrect(:,bestwinip) ones(length(confusion_mi(:,bestwinip)),1)]);
    plot(percorrect(:,1), percorrect(:,1)*b(1) + b(2),'k:');  
    plot([0 1],[log2(20) log2(20)], 'k--');   % Saturation level
    title(sprintf('Confusion MI vs Percent Correct Tau = %d ms R2=%.2f', winSize(bestwinip), stats(1)));
    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     figure(17);
    % Z-dist not corrected Information vs Confusion Information 
   bestwinip = 3;
    plot(confusion_mi(:,bestwinip), zdist_ncmi(:,bestwinip), '+');
    axis([0 5 0 5]);
    ylabel('Z MI (bits/s)');
    xlabel('Confusion MI (bits/s)');
    hold on;
    [b,bint,r,rint,stats] = regress(zdist_ncmi(:,bestwinip), [confusion_mi(:,bestwinip) ones(length(confusion_mi(:,bestwinip)),1)]);
    plot(confusion_mi(:,1), confusion_mi(:,1)*b(1) + b(2),'k:'); 
    plot([0 5], [0 5], 'k');

    title(sprintf('Z MI vs Confusion MI Tau = %d R2=%.2f', winSize(bestwinip), stats(1)));
    hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     figure(18);
    % Z-dist corrected vs Z-dis not corrected 
   bestwinip = 3;
    plot(zdist_ncmi(:,bestwinip), zdist_mi(:,bestwinip), '+');
    axis([0 10 0 10]);
    ylabel('Z Anthro MI (bits/s)');
    xlabel('Z MI (bits/s)');
    hold on;
    [b,bint,r,rint,stats] = regress(zdist_mi(:,bestwinip), [zdist_ncmi(:,bestwinip) ones(length(zdist_ncmi(:,bestwinip)),1)]);
    % plot(zdist_ncmi(:,1), zdist_ncmi(:,1)*b(1) + b(2),'k:');  
    plot([log2(20) log2(20)], [0 10], 'k--');
    plot([0 10], [log2(20) log2(20)], 'k--');
    plot([0 10], [0 10], 'k');
    title(sprintf('Z Antrho MI vs Z MI Tau = %d ms R2=%.2f', winSize(bestwinip), stats(1)));
    hold off;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     figure(19);
    % Z-dist corrected vs Z-dis not corrected 
   bestwinip = 3;
    plot(confusion_mi(:,bestwinip), zdist_mi(:,bestwinip), '+');
    axis([0 10 0 10]);
    ylabel('Z Anthro MI (bits/s)');
    xlabel('Confusion MI (bits/s)');
    hold on;
    [b,bint,r,rint,stats] = regress(zdist_mi(:,bestwinip), [confusion_mi(:,bestwinip) ones(length(zdist_ncmi(:,bestwinip)),1)]);
    % plot(zdist_ncmi(:,1), zdist_ncmi(:,1)*b(1) + b(2),'k:');  
    plot([log2(20) log2(20)], [0 10], 'k--');
    plot([0 10], [log2(20) log2(20)], 'k--');
    plot([0 10], [0 10], 'k');
    title(sprintf('Z Antrho MI vs Z MI Tau = %d ms R2=%.2f', winSize(bestwinip), stats(1)));
    hold off;
    
pause();


end