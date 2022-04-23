corr.vs.diff_coarsen <- function(coarsen.paths = NULL, inversion.out.ext = 'inversion/out',
                                 obs.error = 0, output.path = 'Out/Coarse_Results') {
  
  library(ggallin); library(Metrics)
  
  #First, let's take a look at the XCO2 error
  for(i in 1:length(coarsen.paths)) {
    
    #obtain the aggregation factor
    agg.factor <- as.numeric(gsub('F', '', basename(coarsen.paths[i])))
    
    #list the available SAMs
    SAM.list <- list.files(coarsen.paths[i], full.names = TRUE)
    xco2.path <- list.files(file.path(SAM.list, inversion.out.ext,
                                      'Analysis'), full.names = TRUE)
    
    #list the xco2.files
    all.xco2 <- do.call(rbind, lapply(xco2.path, read.csv))
    all.xco2$agg.factor <- paste0(agg.factor, 'km x ', agg.factor, 'km')
    
    if(i == 1)
      all.factors_xco2 <- all.xco2
    if(i > 1)
      all.factors_xco2 <- rbind(all.factors_xco2, all.xco2)
    
  }
  
  disp.values <- subset(all.factors_xco2, obs > obs.error & obs <= 5)
  disp.values <- subset(all.factors_xco2, obs <= 5)
  omitted.data <-
    100*(1 - nrow(subset(all.factors_xco2, obs <= 5))/nrow(all.factors_xco2))
  
  lm.data <- data.frame(matrix(NA, nrow = 0, ncol = 6))
  names(lm.data) <- c('agg.factor', 'R2', 'slope',
                      'y.int', 'y.sign', 'RMSE')
  for(i in 1:length(unique(all.factors_xco2$agg.factor))) {
    
    tmp <- subset(all.factors_xco2,
                  agg.factor == unique(all.factors_xco2$agg.factor)[i])
    
    lm.sum <- summary(lm(data = tmp,
                         formula = abs(posterior - prior) ~
                           abs(obs - prior)))
    RMSE <- rmse(abs(tmp$obs - tmp$prior),
                 abs(tmp$posterior - tmp$prior))
    
    add.line <- data.frame(agg.factor =
                             unique(all.factors_xco2$agg.factor)[i],
                           R2 = lm.sum$r.squared,
                           slope = lm.sum$coefficients[2,1],
                           y.int = lm.sum$coefficients[1,1],
                           y.sign = NA,
                           RMSE = RMSE)
    lm.data <- rbind(lm.data, add.line)
  }
  
  lm.data$y.sign[lm.data$y.int < 0] <- '-'
  lm.data$y.sign[lm.data$y.int >= 0] <- '+'
  
  
  Corr.Vs.Diff.plot <- ggplot() +
    ggtitle(expression(paste('XCO'[2],
                             ' Correction vs. Initial Difference'))) +
    labs(subtitle = 'Grid Resolution',
         caption = paste0('z>5ppm removed for plotting.', '\n',
                          round(omitted.data, 2), '% of data removed.')) +
    xlab(expression('|'*italic(z)-bold(H)*italic(s[p])*'| [ppm]')) +
    ylab(expression('|'*bold(H)*(hat(italic(s)) - italic(s[p]))*'| [ppm]')) +
    geom_point(data = disp.values, shape = 1, alpha = 0.5,
               color = 'black',
               aes(x = (abs(obs - prior)),
                   y = (abs(posterior - prior)))) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    geom_smooth(data = disp.values, se = FALSE,
                color = 'gold', linetype = 'dashed',
                aes(x = (abs(obs - prior)),
                    y = (abs(posterior - prior)))) +
    geom_abline(data = lm.data, linetype = 'solid',
                color = 'blue',
                aes(slope = slope, intercept = y.int)) +
    geom_text(data = lm.data, hjust = 0, vjust = 1,
              size = 3,
              aes(x = -Inf, y = Inf,
                  label = paste0('\n',
                                 ' y=', round(slope, 2), 'x',
                                 y.sign, round(abs(y.int), 2), '\n',
                                 ' R^2=', round(R2, 2)))) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(. ~ agg.factor, ncol = 2)
  
  ggsave(Corr.Vs.Diff.plot, device = 'jpg',
         height = 4, width = 6, units = 'in',
         filename = file.path(output.path, 'Corr.vs.Diff.jpg'))
  
}