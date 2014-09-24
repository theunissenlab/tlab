function launch_bird_series(n);
addpath('../sound');
data_dir='../StressCalls/';
start;
zf=length(Sbirds);
for i=1:(zf/4);
    offset=(n-1)*(zf/4)+i
    sbird=Sbirds(offset)
    launch_stress_birds(sbird,data_dir);
end;