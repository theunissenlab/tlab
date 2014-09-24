read_anova_file = function(filePath)
{
  ds = read.table(filePath, sep=',', header=TRUE)
  ds$thresh = factor(ds$thresh)
  
  return(ds)
}

plot_marginals = function(ds)
{
  plot(ds$thresh, ds$score)
  plot(ds$model, ds$score) 
  plot(ds$region, ds$score)
  plot(ds$preproc, ds$score)
}

runtime_by_threshold = function(ds)
{
  par(mfrow=c(5, 1))
  rt_all = ds$runtime / 60.0
  rtmax = max(rt_all)
  rtmin = min(rt_all)
  for (l in levels(ds$thresh)) {
    indx = ds$thresh == l
    rt = rt_all[indx]
    rtmean = mean(rt)
    rtstd = sd(rt)
    hist(rt, xlim=c(rtmin, rtmax),
         xlab='Runtime (min)', ylab='Count',
         main=sprintf('Thresh=%s, mean=%0.0f, std=%0.0f', l, rtmean, rtstd))    
  }  
}

iters_by_threshold = function(ds, eindx=TRUE)
{
  par(mfrow=c(5, 1))
  iter_all = ds$num_iters
  itermax = max(iter_all)
  itermin = min(iter_all)
  for (l in levels(ds$thresh)) {
    indx = ds$thresh == l & eindx
    iter = iter_all[indx]
    itermean = mean(iter)
    iterstd = sd(iter)
    hist(iter, xlim=c(itermin, itermax),
         xlab='# of Iterations', ylab='Count',
         main=sprintf('Thresh=%s, mean=%0.0f, std=%0.0f', l, itermean, iterstd))    
  }  
}



perf_by_threshold = function(ds, eindx=TRUE)
{  
  tlvls = levels(ds$thresh)
  par(mfrow=c(length(tlvls), 1))
  perf_all = ds$score
  perfmax = max(perf_all)
  perfmin = min(perf_all)
  for (l in tlvls) {
    indx = ds$thresh == l & eindx
    perf = perf_all[indx]
    perfmean = mean(perf)
    perfstd = sd(perf)
    hist(perf, xlim=c(perfmin, perfmax),
         xlab='Performance Ratio', ylab='Count',
         main=sprintf('Thresh=%s, mean=%0.2f, std=%0.2f', l, perfmean, perfstd))    
  }
}

perf_by_region = function(ds, eindx=TRUE)
{
  rlvls = c('MLd', 'OV', 'L', 'CM')
  par(mfrow=c(length(rlvls), 1))
  perf_all = ds$score
  perfmax = max(perf_all)
  perfmin = min(perf_all)
  for (l in rlvls) {
    indx = ds$region == l & eindx
    perf = perf_all[indx]
    perfmean = mean(perf)
    perfstd = sd(perf)
    hist(perf, xlim=c(0.0, 1.0),
         xlab='Performance Ratio', ylab='Count',
         main=sprintf('Region=%s, mean=%0.2f, std=%0.2f', l, perfmean, perfstd))    
  }  
}

perf_by_model = function(ds, eindx=TRUE)
{
  rlvls = c('linear', 'sepnl_spline', 'binomial', 'poisson', 'leglm')
  par(mfrow=c(length(rlvls), 1))
  perf_all = ds$score
  perfmax = max(perf_all)
  perfmin = min(perf_all)
  for (l in rlvls) {
    indx = ds$model == l & eindx
    perf = perf_all[indx]
    perfmean = mean(perf)
    perfstd = sd(perf)
    hist(perf, xlim=c(0.0, 1.0),
         xlab='Performance Ratio', ylab='Count',
         main=sprintf('Model=%s, mean=%0.2f, std=%0.2f', l, perfmean, perfstd))    
  }  
}

perf_by_preproc = function(ds, eindx=TRUE)
{
  rlvls = c('rawspectrogram', 'lyons', 'surprise')
  par(mfrow=c(length(rlvls), 1))
  perf_all = ds$score
  perfmax = max(perf_all)
  perfmin = min(perf_all)
  for (l in rlvls) {
    indx = ds$preproc == l & eindx
    perf = perf_all[indx]
    perfmean = mean(perf)
    perfstd = sd(perf)
    hist(perf, xlim=c(0.0, 1.0),
         xlab='Performance Ratio', ylab='Count',
         main=sprintf('Preproc=%s, mean=%0.2f, std=%0.2f', l, perfmean, perfstd))    
  }  
}


do_anova = function(ds)
{  
  model = lm(score ~ model + preproc + thresh + model*thresh + model*preproc + preproc*thresh + 0,
             data=ds)
  model_anova = anova(model)
  return(model_anova)
}