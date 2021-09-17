# 3_1_SOLAR Circos globe plot

# Read in SOLAR and CHS mz/cpd key and mwas results
chs_single_matches_w_mwas <- readRDS(fs::path(dir_temp, 
                                                "chs_mz_cpd_pathway_key_w_mwas.rds"))

solar_single_matches_w_mwas <- readRDS(fs::path(dir_temp, 
                                                  "solar_mz_cpd_pathway_key_w_mwas.rds"))

single_matches_w_mwas_final <- inner_join(solar_single_matches_w_mwas, 
                                          chs_single_matches_w_mwas, 
                                          by = c("cpd", "exposure", "exposure2"), 
                                          suffix = c("_solar", "_chs")) %>% 
  mutate(sig = if_else((p_value_solar < 0.05) & 
                         (p_value_chs < 0.05), 
                       "Sig", 
                       "Not-sig"))




############################################
## EXAMPLE FROM DOM
############################################

p_circos_omics <- function(pcor_list = pcor_tar_feat, omics_layer){
  
  # if there's no target features, exit the function early with message
  if(nrow(pcor_list[[paste0("pcor_tar_", omics_layer)]]) == 0){
    
    message("There are no target features for this combination of liver injury outcome and multi-omics layer")
    
  }else{
    
    # omics layers in order
    omics_layer_prefix <- c("pro", "metser", "meturi")
    omics_layer_order <- c("Proteins", "Serum Metabolites", "Urine Metabolites")
    metser_class_order <- sort(unique(lst_metabolome_serum$names_metabolome_serum$mtb_metabolites_serum$class))
    group_order <- c(omics_layer_order[1], metser_class_order, omics_layer_order[3])
    
    # serum metabolite classes for each feature
    metser_class <- lst_metabolome_serum$names_metabolome_serum$mtb_metabolites_serum$class
    names(metser_class) <- clean_punct(clean_feat_aa_ae(clean_feat_labs(lst_metabolome_serum$names_metabolome_serum$mtb_metabolites_serum$feat_label)))
    
    # select parameters for single omics layer used in function
    i_params <- tibble(prefix = omics_layer_prefix, string = omics_layer_order) %>% dplyr::filter(prefix == omics_layer)
    
    # names of target features
    names_tar_feat <- clean_punct(clean_feat_aa_ae(clean_feat_labs(unique(pcor_list[[paste0("pcor_tar_", i_params$prefix)]]$tar_feat))))
    
    # single dataframe for circlize plot
    df_messy <- pcor_list[[paste0("pcor_tar_", i_params$prefix)]] %>% 
      # add metabolite classes to target dataframe
      left_join(by = "all_feat", y = dplyr::rename(select(lst_metabolome_serum$names_metabolome_serum$mtb_metabolites_serum, feat_label, class, common_name), all_feat = feat_label)) %>% 
      # add grouping variable identifying omics layer
      group_feature_prefix(all_feat) %>% 
      # remove correlations within the target omics layer (too distracting)
      dplyr::filter(omics_group != i_params$string) %>%
      # add identity rows for target features
      bind_rows(data.frame(
        tar_feat = names_tar_feat, 
        all_feat = names_tar_feat, 
        estimate = rep(0, length(names_tar_feat)), 
        omics_group = rep(i_params$string, length(names_tar_feat)),
        # assign class if serum metabolite analysis, otherwise leave blank
        class = ifelse(omics_layer == "metser", metser_class[names_tar_feat], "")
      ))
    
    # format data for circlize
    df_details <- df_messy %>% 
      dplyr::mutate(
        # format omics grouping in consistent order, drop unused levels
        omics_group_class = ifelse(omics_group == "Serum Metabolites", class, omics_group) %>% 
          factor(levels = group_order, ordered = TRUE) %>% 
          droplevels(),
        # remove long prefixes from names
        tar_feat = clean_punct(clean_feat_aa_ae(clean_feat_labs(tar_feat))),
        all_feat = clean_punct(clean_feat_aa_ae(clean_feat_labs(all_feat))),
        # indicator for target features, helps with ordering
        target_feat = ifelse(all_feat %in% names_tar_feat, 1, 2) %>% factor(levels = 1:2, labels = c("Target", "Correlation"))
      ) %>% 
      # arrange target features first, then magnitudes of association
      dplyr::arrange(omics_group_class, target_feat, estimate)
    
    # select variables for chord diagram
    df <- df_details %>% select(tar_feat, all_feat, estimate, omics_group, omics_group_class)
    
    # vector of features by group, formatted for circlize
    vec_group <- df$omics_group_class
    names(vec_group) <- df$all_feat
    
    # dataframe of feature colors
    df_grid_col <- df %>% select(all_feat, omics_group, omics_group_class) %>% 
      mutate(grid.col = case_when(
        # replace target features with dark colors
        all_feat %in% clean_feat_aa_ae(clean_feat_labs(names_tar_feat)) ~ colors_fct_omics_dark[i_params$string],
        # add serum metabolites light green color if serum metabolite and class does not match broad omics groups
        omics_group == "Serum Metabolites" & !(omics_group_class %in% unique(.data$omics_group)) ~ colors_fct_omics_light[["Serum Metabolites"]],
        # add light colors for all omics features
        TRUE ~ colors_fct_omics_light[as.character(omics_group_class)]
      ))
    # vector of colors for grid blocks, formatted for circlize
    grid.col <- df_grid_col$grid.col
    names(grid.col) <- df_grid_col$all_feat
    
    # color pal indicating strength of correlation
    #   always set extreme values at ±1 for blue-white-red colors
    col_pal_estimate <- colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red"))
    
    # vector of link colors, formatted for circlize
    link.col <- col_pal_estimate(df$estimate)
    names(link.col) <- df$all_feat
    
    # data frame of group colors
    df_group <- unique(df %>% select(omics_group, omics_group_class)) %>% 
      dplyr::mutate(color = colors_fct_omics_light[.data$omics_group])
    
    # initialize chord diagram
    par(xpd = NA)
    circos.par(start.degree = 90)
    chordDiagram(df, annotationTrack = c("grid"), group = vec_group, grid.col = grid.col, col = link.col, preAllocateTracks = 1, )
    
    # text for serum metabolomics subgroups
    group_names <- function(n){
      # text for serum groups
      group_vec <- vec_group[vec_group == df_group$omics_group_class[n]]
      highlight.sector(names(group_vec), track.index = 1, text = str_to_title(as.character(group_vec[n])), col = df_group$color[n],
                       text.col = "black", facing = "bending", niceFacing = TRUE, text.vjust = "0mm", cex = 0.5, 
                       padding = c(0.075, 0, -0.9, 0))
    }
    lapply(1:length(unique(vec_group)), group_names)
    
    # https://stackoverflow.com/questions/31943102/rotate-labels-in-a-chorddiagram-r-circlize
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      
      circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.15, CELL_META$sector.index, cex = 0.5,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    }, bg.border = NA)
    
    # legend omics layers
    legend(title = "Omics Layer", bty = "n", x = "topright", cex = 0.8, pch = 15, pt.cex = 1,
           col = c(df_group$color), 
           legend = c(str_to_title(as.character(df_group$omics_group_class)))
    )
    
    # legend links exposures
    # legend is responsive, so shows colors for maximum values
    est_max <- max(df$estimate)
    est_med_u <- median(dplyr::filter(df, estimate >= 0)$estimate)
    est_med_l <- median(dplyr::filter(df, estimate < 0)$estimate)
    est_min <- min(df$estimate)
    range_scale <- c(1, est_max, est_med_u, 0, est_med_l, est_min, -1)
    legend(title = "Correlation", bty = "n", x = "bottomright", cex = 0.8, pch = 15, pt.cex = 1, xjust = 1,
           col = col_pal_estimate(range_scale),
           legend = paste0(c(rep(" ", 4), rep("", 3)),
                           paste(sep = "   ",
                                 formatC(range_scale, digits = 2, format = "f"), 
                                 c("", "Maximum", "Median Positive", "None", "Median Negative", "Minimum", ""))
           )
    )
    
    # info details
    # circos.info()
    
    # close commands of the circos plot
    circos.clear()
    
  }
  
}
