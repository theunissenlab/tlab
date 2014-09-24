%%% Testing the Read all spikes  


h5Path = '/auto/fdata/fet/julie/h5/BlaBro09xxF/Site1_L1500R1500_18.h5';
snipPath = '/auto/fdata/julie/Exported_data_tank/BlaBro09xxF/Site1_depth1500_Call1 waves.f32';


h5_plot_all_spikes(h5Path, 'Call1', 0.05, 1, snipPath);

