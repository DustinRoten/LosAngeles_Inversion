collect.bias.correction <- function(bias.paths = NULL, gg.map = NULL, plot.output = NULL) {
  
  library(stringr)
  
  #get the category titles; remove any special '_' categories
  cat.titles <- read.csv('ext/defined_vulcan_sectors.csv')[,1]
  cat.titles <- grep(x = cat.titles, pattern = '_',
                     value = TRUE, invert = TRUE)
  cat.titles <- gsub(' ', '.', cat.titles)
  
  #create the dataframe to store category-based information in
  lambda_df <- data.frame(matrix(NA, nrow = 0, ncol = 5))
  names(lambda_df) <- c('LPS', 'bias', 'SAMs', 'category', 'lambda')
  
  #create the dataframe to store total-based information in 
  lambda_df.tot <- data.frame(matrix(NA, nrow = 0, ncol = 4))
  names(lambda_df.tot) <- c('LPS', 'bias', 'SAMs', 'lambda')
  for(i in 1:length(bias.paths)) {
    
    path <- bias.paths[i]
    
    SAM.list <- list.files(path, full.names = TRUE)
    
    params <- str_split_fixed(basename(path),
                              pattern = '_',
                              n = 2)
    
    bias <- as.numeric(params[1]); LPS <- params[2]
    
    for(j in 1:length(SAM.list)) {
      
      #read in the values of K (Hs_prior)
      xco2.path <- list.files(file.path(SAM.list[1:j],
                                        'inversion/out/Analysis'),
                              pattern = 'all_xco2.csv',
                              full.names = TRUE)
      all.xco2 <- do.call(rbind, lapply(xco2.path, read.csv))
      
      #get the appropriate columns from the K dataframe
      #K contains the category-based prior XCO2 data
      #K.tot is the total prior
      col.idx <- which(names(all.xco2) %in% cat.titles == TRUE)
      K <- as.matrix(all.xco2[,col.idx])
      K.tot <- as.matrix(all.xco2$prior)
      
      #Obtain the observation values
      z <- all.xco2$obs
      
      #prepare new aggregated R matrix
      agg.R <- matrix(0, nrow = dim(K)[1], ncol = dim(K)[1])
      
      #read in the R matrix
      R.path <- list.files(file.path(SAM.list[1:j], 'inversion/out'),
                           pattern = 'R.rds',
                           full.names = TRUE)
      
      row.end <- col.end <- 0
      for(k in 1:length(R.path)) {
        R <- readRDS(R.path[k])
        rows <- nrow(R); cols <- ncol(R)
        
        agg.R[(row.end+1):(row.end+rows),
              (col.end+1):(col.end+cols)] <- R
        
        row.end <- (row.end+rows)
        col.end <- (col.end+cols)
      } #closes the aggregation of R matrices
      
      #agg matrix becomes R
      R <- agg.R
      
      #build the s_p matrix (with categories)
      lambda <- matrix(bias, nrow = 4, ncol = 1)
      zero.matrix <- matrix(0, nrow = 4, ncol = 4)
      diag(zero.matrix) <- 1-bias
      s_p <- zero.matrix
      
      #build the s_p matrix (with total)
      lambda.tot <- bias
      s_p.tot <- 1-bias
      
      #fancy math (with categories)
      lambda_hat <-
        lambda + s_p %*% t(K) %*% solve(K %*% s_p %*% t(K) + R) %*% (z - K %*% lambda)
      
      #fancy math (with total)
      lambda_hat.tot <-
        lambda.tot + s_p.tot %*% t(K.tot) %*%
        solve(K.tot %*% s_p.tot %*% t(K.tot) + R) %*% (z - K.tot %*% lambda.tot)
      
      #add to the category-based dataframe
      add.line <-
        data.frame(LPS = LPS, bias = bias, SAMs = j,
                   lambda = cbind(cat.titles,
                                  as.numeric(lambda_hat)))
      names(add.line) <- names(lambda_df)
      lambda_df <- rbind(lambda_df, add.line)
      
      #add to the total-based dataframe
      add.line <-
        data.frame(LPS = LPS, bias = bias, SAMs = j,
                   lambda = lambda_hat.tot)
      names(add.line) <- names(lambda_df.tot)
      lambda_df.tot <- rbind(lambda_df.tot, add.line)
      
    } #closes the number of SAMs to use
  } #closes the bias loop
  
  #fix a stupid formatting error
  #this can be done more elegantly later
  lambda_df$lambda <- as.numeric(lambda_df$lambda)
  lambda_df$category <- gsub('\\.', ' ', lambda_df$category)
  lambda_df$LPS[lambda_df$LPS == 'LPS'] <- 'Large Point Source Error'
  lambda_df$LPS[lambda_df$LPS == 'noLPS'] <- 'No Large Point Source Error'
  lambda_df.tot$LPS[lambda_df.tot$LPS == 'LPS'] <- 'Large Point Source Error'
  lambda_df.tot$LPS[lambda_df.tot$LPS == 'noLPS'] <- 'No Large Point Source Error'
  
  ### Plot the percentage of each category ###
  lps.idx <- grep(pattern = '_LPS', bias.paths)
  
  SAM.list <- list.files(bias.paths[lps.idx], full.names = TRUE)
  xco2.path <- list.files(file.path(SAM.list, 'inversion/out/Analysis'),
                          full.names = TRUE,
                          pattern = 'all_xco2.csv')
  all.xco2 <- do.call(rbind, lapply(xco2.path, read.csv))
  
  #generate xco2 percentage plot
  col.idx <- which(names(all.xco2) %in% cat.titles == TRUE)
  category.percentages <- as.data.frame(matrix(NA, nrow = 0, ncol = 4))
  names(category.percentages) <- c('lon', 'lat', 'Percent', 'Category')
  for(i in 1:length(col.idx)) {
    add.lines <- cbind(all.xco2$lon, all.xco2$lat,
                       100*all.xco2[,col.idx[i]]/all.xco2$prior,
                       cat.titles[i])
    names(add.lines) <- names(category.percentages)
    category.percentages <- rbind(category.percentages, add.lines)
  }
  
  #clean up the dataframe
  category.percentages[,1] <- as.numeric(category.percentages[,1])
  category.percentages[,2] <- as.numeric(category.percentages[,2])
  category.percentages[,3] <- as.numeric(category.percentages[,3])
  category.percentages[,4] <- gsub('\\.', ' ', category.percentages[,4])
  
  #remove the NaN crap
  category.percentages <- subset(category.percentages,
                                 !is.nan(category.percentages[,3]) &
                                   category.percentages[,3] > 0)
  
  #renaming just in case
  names(category.percentages) <- c('lon', 'lat', 'Percent', 'Category')
  
  #grid the uneven data with a custom function.
  custom_category.percentages <-
    custom.grid(lons = category.percentages$lon,
                lats = category.percentages$lat,
                values = category.percentages$Percent,
                layers = category.percentages$Category,
                resolution = 0.05)
  names(custom_category.percentages) <- names(category.percentages)
  
  ### Plots ###
  #generate the category box plot
  xco2.breakdown <- ggplot() +
    ggtitle(expression(paste('Composition of XCO'[2], ' Soundings'))) +
    labs(subtitle = '(OSSE Pseudo-Observations)') +
    xlab('Percent [%]') +
    ylab('Emission Category') +
    geom_boxplot(data = category.percentages,
                 aes(x = Percent, y = Category,
                     fill = Category)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'none')
  
  ggsave(xco2.breakdown,
         device = 'jpg', height = 4, width = 6, units = 'in',
         filename = file.path(plot.output, 'XCO2_Distribution.jpg'))
  
  #generate lambda plot  
  lambda.plot <- ggplot() +
    ggtitle('Bias Correction') +
    xlab('Number of SAMs') +
    ylab(expression(hat(lambda))) +
    geom_hline(aes(yintercept = 4), color = 'gray') + #hard-coded value of 4
    geom_line(data = lambda_df.tot,
              linetype = 'dashed',
              color = 'black',
              aes(x = SAMs, y = lambda,
                  group = 1)) +
    geom_line(data = lambda_df,
              aes(x = SAMs, y = lambda,
                  color = category)) +
    geom_point(data = lambda_df,
               aes(x = SAMs, y = lambda,
                   color = category)) +
    scale_color_discrete(name = 'Category:') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom') +
    facet_grid(. ~ LPS)
  
  ggsave(lambda.plot, device = 'jpg',
         filename = file.path(plot.output, 'Lambda_Plot.jpg'),
         height = 3, width = 6, units = 'in')
  
  #generate the spatially aggregated data
  custom_category.percentages <-
    subset(custom_category.percentages, Percent > 0)
  
  XCO2.spatial <- ggmap(gg.map) +
    ggtitle(expression(paste('Averaged Contributions to XCO'[2], ' Soundings'))) +
    xlab('Longitude') +
    ylab('Latitude') +
    geom_tile(data = custom_category.percentages, alpha = 0.7,
               aes(x = lon, y = lat, fill = Percent)) +
    scale_fill_gradient2(low = 'yellow', mid = 'orange', high = 'red',
                        midpoint = 50, name = 'Percent [%]') +
    guides(fill = guide_colourbar(ticks.colour = "black",
                                  ticks.linewidth = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 1),
           name = 'Percent [%]') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(0.5, 'in'),
          axis.text.x = element_text(hjust = 1, angle = 45)) +
    facet_wrap(. ~ Category, ncol = 2)
  
  ggsave(XCO2.spatial, device = 'jpg',
         filename = file.path(plot.output, 'XCO2_Spatial.jpg'),
         height = 6, width = 6, units = 'in')
  
} #closes the function
