addpath /auto/fhome/anne/sarah_thane_info/scripts
areadir={'fieldL', 'MLD', 'cHV'};
stimdir={'con', 'songrip', 'flatrip'};
%stimdir={'conwn','flatrip100', 'conflatrip100', 'noise'};
paramfilenames={'infosnr.dat', 'infoolw4.dat', 'infoolw12.dat', 'infogamma.dat', 'infoparams.dat', 'z_scores.dat', 'cc_twohalves.dat'};

     
        filename=['database' areadir{ar} stimdir{st} '.txt']; 
        outputdir=['/auto/fhome/anne/sarah_thane_info/' areadir{ar} '/' stimdir{st}] ;
           eval(['cd ' outputdir]);
            for pf=1:length(paramfilenames)
               eval(['fid' num2str(pf) '=fopen(paramfilenames{' num2str(pf) '}, ''w'');']);
            end
            fid = fopen(filename, 'r');     
            databaseline=fgetl(fid);
            goflag=1;
            if or(isempty(databaseline), databaseline==-1)
                goflag=0;
            end
            onset=200;
            
            while goflag
                
                birdname=getnextword(databaseline, 1);
                neuronnumber=getnextword(databaseline,2);
                stimtype=getnextword(databaseline,3);
                nfiles=getnextword(databaseline,4);
                strfpath=getnextword(databaseline, 5);
                spiketrain=[];
                
                eval(['cd /auto/fdata/anne/sarah_thane_info/' birdname]);
                for nf=1:str2num(nfiles)
                    eval(['spike=load(''spike' neuronnumber stimtype num2str(nf) '.dat'');']);
                    if ~isempty(spiketrain)
                        if size(spike,1) ~=size(spiketrain,1);
                            display(['spike' neuronnumber stimtype num2str(nf) ' skipped']);
                        else
                            spiketrain=[spiketrain spike(:,onset+1:end)];
                        end
                    else
                         spiketrain=[spiketrain spike(:,onset+1:end)];
                    end
                   
                end
                %if stim is conwn or conflatrip100, the strfdir is just
                %con.
                if or(strcmp(stimtype, 'conwn')==1, strcmp(stimtype, 'conflatrip100')==1)
                    stimtype='con';
                end
                %for thane's STRF directory, noise=whitenoise
                if and(strmatch('/auto/fdata/thane/STRF', strfpath), strcmp(stimtype, 'noise')==1)
                    stimtype='whitenoise';
                end
                
                eval(['cd ' strfpath neuronnumber '/' stimtype ]);
                fres=read_newresults2('newresults2');
                zscorenoonset=fres.zscore_noonset(1);
                fprintf(fid6, '%f\n', zscorenoonset);
                    
                [cc_two_halves_corrval, cc_two_halves_best, cc_two_halves_constant]=cc_twohalves(spiketrain, 21);
                fprintf(fid7, '%f %f\n', cc_two_halves_best, cc_two_halves_constant);

                
                eval(['cd ' outputdir]);
                
                
                [info, infoup, infodown]=make_snrinfo(spiketrain); 
                fprintf(fid1, '%f %f %f\n', info, infoup, infodown);
                
                [info, noiseentropy, totalentropy]=make_psth_onelw4noonset(spiketrain);
                fprintf(fid2, '%f %f %f\n', info, noiseentropy, totalentropy);
                
                [info, noiseentropy, totalentropy]=make_psth_onelw12noonset(spiketrain);
                fprintf(fid3, '%f %f %f\n', info, noiseentropy, totalentropy);
                
                [info, noiseentropy, totalentropy]=make_psth_gammanoonset(spiketrain); 
                fprintf(fid4, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f  \n', info, noiseentropy, totalentropy);
                
                
                [gammaconst, gammalb, gammaub, gammaerror, rategamma, rategammalb, rategammaub, rategammaerror, bandwidth, ratemean, avgff]=make_gamma_noonset(spiketrain);
                fprintf(fid5, '%f %f %f %f %f %f %f %f %f %f %f\n',gammaconst, gammalb, gammaub, gammaerror, rategamma, rategammalb, rategammaub, rategammaerror, bandwidth, ratemean, avgff);
                
                databaseline=fgetl(fid);
                %skip spaces between lines.
                while isempty(databaseline)
                    databaseline=fgetl(fid);
                end
                if databaseline==-1
                    goflag=0;
                end
            end
            
            for pf=1:length(paramfilenames)
                eval(['fclose(fid' num2str(pf) ');']);
            end
       
