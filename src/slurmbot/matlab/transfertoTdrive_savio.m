function [Status]=transfertoTdrive_savio(From_file, To_file, KeepF)
%% This function transfer a file from savio to tdrive and erase the file under savio
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

% Input variables:   From_file: path and name of the file to be transfered
%                    To_file: path and name of the file location on Tdrive
%                    KeepF: set to 1 to keep the file localy, set to 0 to
%                    remove it
    
if nargin<3
    KeepF = 0;
end

Filecopied=0;
while Filecopied==0
    CommandSystemtransfer1=['ssh dtn.brc scp ' From_file ' TheunissenLab:' To_file];
    fprintf(1,'Attempting transfer 1 with command: %s\n',CommandSystemtransfer1);
    [Status]=system(CommandSystemtransfer1);
    if (~Status) || (~KeepF) 
        fprintf(1,'file correctly transfered');
        fprintf(1,'Now remove file from Savio');
        system(['rm ' From_file])
        Filecopied=1;
    end
end

end
