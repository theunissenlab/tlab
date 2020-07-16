function recover_best_strf(calc_dir,best_tol_val,best_std_val)

%% Output file
outfile = fullfile(calc_dir,'best_strf.mat');
if exist(outfile)
    error('File %s already exists.',outfile)
end

%% Load calculation data
g = load(fullfile(calc_dir,'Global_Variables.mat'));
i = load(fullfile(calc_dir,'info_r_result.mat'));

%% Load best STRF
best_tol = find(g.Tol_val==best_tol_val);
best_std = find(g.Std_val==best_std_val);
s = load(fullfile(calc_dir,sprintf('strfResult_Tol%i.mat',best_tol)));

%% Assign values to save
cc_ratio_const = i.cc_ratio_const{best_tol}{best_std};
cc_ratio_max = i.cc_ratio_max{best_tol}{best_std};
cc_spike_pre_best = i.cc_spike_pre_best{best_tol}{best_std};
cc_spike_pre_constant = i.cc_spike_pre_constant{best_tol}{best_std};
cc_two_halves_constant = i.cc_two_halves_constant{best_tol}{best_std};
cc_two_halves_tmax = i.cc_two_halves_tmax{best_tol}{best_std};
info = i.info{best_tol}{best_std};
infodown = i.infodown{best_tol}{best_std};
infodownpre = i.infodownpre{best_tol}{best_std};
infopre = i.infopre{best_tol}{best_std};
infoup = i.infoup{best_tol}{best_std};
infouppre = i.infoup{best_tol}{best_std};
psth_smoothconst = 0;
smoothVect = i.smoothVect;
strf = s.STRF_Cell;
tmax_pre = i.tmax_pre{best_tol}{best_std};
tmax_ratio = i.tmax_ratio{best_tol}{best_std};

%% Save output

save(outfile,'cc_ratio_const','cc_ratio_max','cc_spike_pre_best',...
     'cc_spike_pre_constant','cc_two_halves_tmax','info','infodown',...
     'infodownpre','infopre','infoup','infouppre','psth_smoothconst',...
     'smoothVect','strf','tmax_pre','tmax_ratio');