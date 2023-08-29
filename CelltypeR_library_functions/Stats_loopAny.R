run_stats2 <- function(input_df, group_cols = c("Sample", "CellType", "Marker"),
                       value_col = "value",
                       stat_type = "ANOVA",
                       id1 = "Genotype",
                       id2 = NULL,
                       use_means = TRUE,
                       loop_by = "none") {
  aov.l <- list() # for ANOVA output
  tuk.l <- list() # for TUKEY output
  if (use_means) {
    df_means <- get_means(input_df, group_cols, value_col)
    
  } else {
    df_means <- input_df
    names(df_means)[names(df_means) == 'value'] <- 'expression'
  }
  
  if (stat_type == "ANOVA") {
    if(loop_by== "none"){
    one_way <- aov(expression ~ df_means[[id1]], data = df_means)
    aov.l <- summary(one_way)
    # now the posthoc test
    tukey <- TukeyHSD(one_way)
    tuk.l <- tukey
    aov_df <- data.frame(
      IndependentVariable = id1,
      Fvalue = aov.l[[1]][["F value"]][1],
      pValue = aov.l[[1]][["Pr(>F)"]][1],
      Df = aov.l[[1]][["Df"]][1],
      SumSqEffect = aov.l[[1]][["Sum Sq"]][1],
      MeanSqEffect = aov.l[[1]][["Mean Sq"]][1],
      SumSqError = aov.l[[1]][["Sum Sq"]][2],
      MeanSqError = aov.l[[1]][["Mean Sq"]][2]
    )
    print("ANOVA results for all celltypes and markers combined")
    rownames <- row.names(tuk.l[[1]])
    tuk_df  <- data.frame(
      IndependentVariable = id1,
      Contrast = rownames,
      tuk.l[[1]], row.names = NULL)
    print("TukeyHSD results for all celltypes and markers combined")
  }
    else{
      var_list <- unique(df_means[[loop_by]])
      # to store outputs and format in a readable way
      for (i in var_list) {
        df <- df_means %>% filter(loop_by == i)
        one_way <- aov(expression ~ df[[id1]], data = df)
        output <- summary(one_way)
        aov.l[[as.character(i)]] <- output # Append output to list
        # now the posthoc test
        tukey <- TukeyHSD(one_way)
        tuk.l[[as.character(i)]] <- tukey
      }
      aov_df <- do.call(rbind, lapply(seq_along(aov.l), function(i) {
        data.frame(Celltype = names(aov.l)[i],
                   IndependentVariable = id1,
                   Fvalue = aov.l[[i]][[1]][["F value"]][1],
                   pValue = aov.l[[i]][[1]][["Pr(>F)"]][1],
                   Df = aov.l[[i]][[1]][["Df"]][1],
                   SumSqEffect = aov.l[[i]][[1]][["Sum Sq"]][1],
                   MeanSqEffect = aov.l[[i]][[1]][["Mean Sq"]][1],
                   SumSqError = aov.l[[i]][[1]][["Sum Sq"]][2],
                   MeanSqError = aov.l[[i]][[1]][["Mean Sq"]][2]
        )
      }
      ))
      print(paste("ANOVA results for each cell type comparing ", id1))
      tuk_df <- do.call(rbind, lapply(seq_along(tuk.l), function(i){
        rownames <- row.names(tuk.l[[i]][[1]])
        data.frame(Celltype = names(tuk.l)[i],
                   IndependentVariable = id1,
                   Contrast = rownames,
                   tuk.l[[i]][[1]], row.names = NULL)
      }))
      print(paste("TukeyHSD results for each ", loop_by," comparing ", id1))
    } 
      aov_df <- do.call(rbind, lapply(seq_along(aov.l), function(i) {
        data.frame(Marker = names(aov.l)[i],
                   IndependentVariable = id1,
                   Fvalue = aov.l[[i]][[1]][["F value"]][1],
                   pValue = aov.l[[i]][[1]][["Pr(>F)"]][1],
                   Df = aov.l[[i]][[1]][["Df"]][1],
                   SumSqEffect = aov.l[[i]][[1]][["Sum Sq"]][1],
                   MeanSqEffect = aov.l[[i]][[1]][["Mean Sq"]][1],
                   SumSqError = aov.l[[i]][[1]][["Sum Sq"]][2],
                   MeanSqError = aov.l[[i]][[1]][["Mean Sq"]][2]
        )
      }
      ))
      print(paste("ANOVA results for each marker comparing ", id1))
      tuk_df <- do.call(rbind, lapply(seq_along(tuk.l), function(i){
        rownames <- row.names(tuk.l[[i]][[1]])
        data.frame(Marker = names(tuk.l)[i],
                   IndependentVariable = id1,
                   Contrast = rownames,
                   tuk.l[[i]][[1]], row.names = NULL)}))
      print(paste("Tukey results for each marker comparing ", id1))
  } 
  
  
  }else if (stat_type == "ANOVA2") {
    if (!is.null(loop_by)) {
      var_list <- unique(df_means$[[loop_by]])
      # to store outputs and format in a readable way
      for (i in var_list) {
        df <- df_means %>% filter(.data[[loop_by]] == i)
        formula <- as.formula(paste0("expression ~ ", id1, "*", id2))
        two_way <- aov(formula, data = df)
        output <- summary(two_way)
        aov.l[[as.character(i)]] <- output # Append output to list
        # now the posthoc test
        tukey <- TukeyHSD(two_way)
        tuk.l[[as.character(i)]] <- tukey
      }
      # make a dataframe of the 2 way anova results
      aov_df <- do.call(rbind, lapply(seq_along(aov.l), function(i) {
        data.frame(
          rownames <- row.names(aov.l[[i]][[1]]),
          data.frame(Celltype = names(aov.l)[i],
                     Contrast = rownames,
                     aov.l[[i]][[1]], row.names = NULL
                     
          )
        )
      }))
      # get ride of the residual which will confuse things and the extra row added
      aov_df <- aov_df[2:8]
      # fix the column names
      colnames(aov_df) <- c("Celltype","Contrast","Df","SumSq","MeanSq",
                            "Fvalue","Pvalue")
      # remove extra spaces
      aov_df <- aov_df %>% mutate_all(trimws)
      # remove residuals
      aov_df <- aov_df[aov_df$Contrast != "Residuals", ]
      
      print(paste("ANOVA 2way results for each cell type comparing ", id1,
                  "and ", id2))
      # make dataframes from each of the tukey results and put them into a list
      tuk_summary.l <- list()
      
      for(J in 1:length(tuk.l[[1]])){
        tuk_df <- do.call(rbind, lapply(seq_along(tuk.l), function(i){
          rownames <- row.names(tuk.l[[i]][[J]])
          data.frame(Celltype = names(tuk.l)[i],
                     Contrast = names(tuk.l[[i]])[J],
                     Subgroups = rownames,
                     tuk.l[[i]][[J]], row.names = NULL)
        }))
        # put the dataframe in a list
        tuk_summary.l[[as.character(names(tuk.l[[i]])[J])]] <- tuk_df
      }
      # filter the interaction dataframe to have only df1 contrasts
      dft <- tuk_summary.l[[3]]
      filtered_df <- dft %>%
        filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(2, 4)) %>%
                 apply(2, function(x) all(x == x[1])))
      # add the filtered dataframe for only id2 pairs
      tuk_summary.l[[paste("Interactions_",id1)]] <- filtered_df
      # filter the interaction dataframe
      filtered_df2 <- dft %>%
        filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(1, 3)) %>%
                 apply(2, function(x) all(x == x[1])))
      # add the filtered dataframe to have matching id1 or id2 contrasts
      tuk_summary.l[[paste("Interactions_",id2)]] <- filtered_df2
      
      print(paste("TukeyHSD results for each cell type comparing ", id1,
                  "and ", id2))
      output_list <- list(ANOVA = aov_df,TukeyHSD = tuk_summary.l)
      
    } else if (loop_by == "Marker") {
      var_list <- unique(df_means$Marker)
      # to store outputs and format in a readable way
      for (i in var_list) {
        df <- df_means %>% filter(Marker == i)
        formula <- as.formula(paste0("expression ~ ", id1, "*", id2))
        two_way <- aov(formula, data = df)
        output <- summary(two_way)
        aov.l[[as.character(i)]] <- output # Append output to list
        # now the posthoc test
        tukey <- TukeyHSD(two_way)
        tuk.l[[as.character(i)]] <- tukey
      }
      # make a dataframe of the 2 way anova results
      aov_df <- do.call(rbind, lapply(seq_along(aov.l), function(i) {
        data.frame(
          rownames <- row.names(aov.l[[i]][[1]]),
          data.frame(Marker = names(aov.l)[i],
                     Contrast = rownames,
                     aov.l[[i]][[1]], row.names = NULL
                     
          )
        )
      }))
      # get ride of the residual which will confuse things and the extra row added
      aov_df <- aov_df[2:8]
      # fix the column names
      colnames(aov_df) <- c("Marker","Contrast","Df","SumSq","MeanSq",
                            "Fvalue","Pvalue")
      # remove extra spaces
      aov_df <- aov_df %>% mutate_all(trimws)
      # remove residuals
      aov_df <- aov_df[aov_df$Contrast != "Residuals", ]
      
      print(paste("ANOVA 2way results for each Marker comparing ", id1,
                  "and ", id2))
      # make dataframes from each of the tukey results and put them into a list
      tuk_summary.l <- list()
      
      for(J in 1:length(tuk.l[[1]])){
        tuk_df <- do.call(rbind, lapply(seq_along(tuk.l), function(i){
          rownames <- row.names(tuk.l[[i]][[J]])
          data.frame(Marker = names(tuk.l)[i],
                     Contrast = names(tuk.l[[i]])[J],
                     Subgroups = rownames,
                     tuk.l[[i]][[J]], row.names = NULL)
        }))
        # put the dataframe in a list
        tuk_summary.l[[as.character(names(tuk.l[[i]])[J])]] <- tuk_df
      }
      # filter the interaction dataframe to have only df1 contrasts
      dft <- tuk_summary.l[[3]]
      filtered_df <- dft %>%
        filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(2, 4)) %>%
                 apply(2, function(x) all(x == x[1])))
      # add the filtered dataframe for only id2 pairs
      tuk_summary.l[[paste("Interactions_",id1)]] <- filtered_df
      # filter the interaction dataframe
      filtered_df2 <- dft %>%
        filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(1, 3)) %>%
                 apply(2, function(x) all(x == x[1])))
      # add the filtered dataframe to have matching id1 or id2 contrasts
      tuk_summary.l[[paste("Interactions_",id2)]] <- filtered_df2
      
      print(paste("TukeyHSD results for each Marker comparing ", id1,
                  "and ", id2))
      output_list <- list(ANOVA = aov_df,TukeyHSD = tuk_summary.l)
    } else {
      df <- df_means
      formula <- as.formula(paste0("expression ~ ", id1, "*", id2))
      two_way <- aov(formula, data = df)
      aov.l <- summary(two_way)
      # make the df of the ANOVA two way results
      rownames <- row.names(aov.l[[1]])
      aov_df <- data.frame(
        Contrast = rownames,
        aov.l[[1]], row.names = NULL
        
      )
      # get ride of the residual which will confuse things and the extra row added
      # fix the column names
      colnames(aov_df) <- c("Contrast","Df","SumSq","MeanSq",
                            "Fvalue","Pvalue")
      # remove extra spaces
      aov_df <- aov_df %>% mutate_all(trimws)
      # remove residuals
      aov_df <- aov_df[aov_df$Contrast != "Residuals", ]
      
      # now the posthoc test
      tukey <- TukeyHSD(two_way)
      tuk.l <- tukey
      print("ANOVA results for all celltypes and markers combined")
      # get the Tukey outputs
      for(J in 1:length(tuk.l)){
        rownames <- row.names(tuk.l[[J]])
        tuk_df <- data.frame(
          Contrast = names(tuk.l[J]),
          Subgroups = rownames,
          tuk.l[[J]], row.names = NULL)
        # put the dataframe in a list
        tuk_summary.l[[as.character(names(tuk.l)[J])]] <- tuk_df
        dft <- tuk_summary.l[[3]]
        filtered_df <- dft %>%
          filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(2, 4)) %>%
                   apply(2, function(x) all(x == x[1])))
        # add the filtered dataframe for only id2 pairs
        tuk_summary.l[[paste("Interactions_",id1)]] <- filtered_df
        # filter the interaction dataframe
        filtered_df2 <- dft %>%
          filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(1, 3)) %>%
                   apply(2, function(x) all(x == x[1])))
        # add the filtered dataframe to have matching id1 or id2 contrasts
        tuk_summary.l[[paste("Interactions_",id2)]] <- filtered_df2
      }
      output_list <- list(ANOVA = aov_df,TukeyHSD = tuk_summary.l)
    }
    
    return(output_list)
  }
}

#

# original get means function

get_means <- function(df, group_cols, value_col) {
  df_means <- df %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(expression = mean(value)) %>%
    distinct(across(all_of(group_cols)), expression, .keep_all = TRUE) %>%
    select(-value)
  return(df_means)
}

