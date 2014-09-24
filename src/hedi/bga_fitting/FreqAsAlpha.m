function [freq_a,slope,intercept]=FreqAsAlpha(alpha0);
addpath ../bga_model;
load ../bga_model/data/stf.fun.dat;
[slope,intercept,freq_a]=FindIso(stf_fun,alpha0);

