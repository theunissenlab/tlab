function [Data]=loadfromTdrive_savio(From_file,Storage_path, KeepF)
%% This function transfer a file from tdrive to savio, load it as Data and erase the file
% Note that for this function to work properly you have to configure
% correctly your connexion between savio and tdrive networks by:
% 1.creating an RSA set of keys on savio in ~/.ssh (e.g. id_rsasavio and
% id_rsasavio.pub) and placing the public key (id_rsasavio.pub) in your 
% home directory on Tdrive cluster. 
% 2.configure the connexion by adding to the config file (~/.ssh/config) on
% Savio:
% Host TheunissenLab
%	User username
%	Port 22
%	IdentityFile ~/.ssh/id_rsasavio
%	HostName jlg.berkeley.edu

% Input variables:   From_file: path and name of the file to be loaded
%                    Storage_path: folder where the file should temporary
%                    be stored. Default is home directory (~)
%                    KeepF: set to 1 to keep the file localy, set to 0 to
%                    remove it
    
if nargin<3
    KeepF = 0;
end

if nargin<2
    InitialPath = pwd;
    cd ~
    Storage_path = pwd; %Use your home directory as a default temporary storage place
    cd(InitialPath);
end

%% Transfer file locally
% Here we're using a data transfer node to transfer the data in the home
% directory
[~,name,ext] = fileparts(From_file);
TempFile = fullfile(Storage_path, [name ext]);
Filecmd0 = ['ssh dtn.brc scp TheunissenLab:' From_file ' ' TempFile];
[status,cmdout]=system(Filecmd0);


%% Load file
% Because synching of the home directory is not optimized in savio, it can
% takes time for a file to be seen by another node (like computing one)
if ~status % the ssh and scp commands were successful
    Wait4F=1;
    tic
    while Wait4F==1
        try
            if strcmp(ext, '.mat')
                Data=load(TempFile);
                Wait4F=0;
            elseif strcmp(ext, '.h5')
                Data = read_unit_h5file(TempFile, 'r');
            end
        catch ME
        end
        ElapsedTime = toc;
        if ElapsedTime>10*60
            error('Times out!! It takes more than 10 min to see the file\nsomething must be wrong with data transfer\nThe error message from load cmd is:\n%s', ME);
        end
    end
else
    error('The file %s could not successfully be transfered to savio\nThe error is: %s\n', From_file, cmdout);
end

%% remove file 
if ~KeepF
    Filecmd1 = ['rm ' TempFile];
    system(Filecmd1);
end

end
