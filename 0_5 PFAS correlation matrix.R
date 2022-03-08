# Correlation plots
library(correlation)
library(GGally)
source(here::here("!functions.R"))
sol_g <- grid.grabExpr(print(sol_corplot))
chs_g <- grid.grabExpr(print(chs_corplot))

# Function:
cor_func <- function(data, mapping, method, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor(x, y, method=method, use='complete.obs')
  
  
  ggally_text(
    label = formatC(corr, format = "f", 2), 
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black',
    ...
  )
}

# SOLAR  ------------------------------------------------------------------
## Create SOLAR exposure data frame --------------
exposures_continuous <-  c("pfos", "pfhxs", "pfhps", "pfoa", "pfna", "pfda")
# sol_chs
# Change 0/1s to factor and remove PFAS with 100% Below LOD (in solar, removed PFHxA)
solar_exposure <- exposure_outcome$solar %>% 
  select(wave, all_of(exposures_continuous))%>% 
  rename_all(toupper) %>% 
  rename(PFHxS = PFHXS, 
         PFHpS = PFHPS)


chs_exposure <- exposure_outcome$chs %>% 
  select(all_of(exposures_continuous))%>% 
  rename_all(toupper) %>% 
  rename(PFHxS = PFHXS, 
         PFHpS = PFHPS)

# Make Final Plots -----------------------------------------------------
# Solar 
sol_corplot <- ggpairs(solar_exposure, columns = 2:7,
                       upper = list(continuous = wrap(cor_func,
                                                      method = "spearman", 
                                                      size = 5)), 
                       lower = list(continuous = wrap("points",alpha = .5, size=.75))) + 
  xlab("Plasma Concentration (µg/L)")  +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12), 
        axis.text.y = element_text(size=10))

# CHS
chs_corplot <- ggpairs(chs_exposure, 
                       upper = list(continuous = wrap(cor_func,
                                                      method = "spearman", 
                                                      size = 5)), 
                       lower = list(continuous = wrap("points",alpha = .5, size=.75))) + 
  xlab("Plasma Concentration (µg/L)") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        axis.text.y = element_text(size=10))



# Combine SOLAR and CHS figs ------------------
sol_g <- grid.grabExpr(print(sol_corplot))
chs_g <- grid.grabExpr(print(chs_corplot))


(fig1 <- plot_grid(NULL, NULL, 
                   sol_g, 
                   chs_g,
                   labels = c("A. SOLAR Cohort", 
                              "B. CHS Cohort", 
                              "", ""),
                   rel_widths = c(1,1),
                   rel_heights = c(.02,1),
                   hjust = 0, vjust = 1,
                   label_size = 14,
                   # align = "vh",
                   axis = "lr",
                   nrow = 2))


## Save Correlation Plot All PFAS -------------------------
ggsave(plot = fig1,
       filename = path(dir_reports, 
                       "SOLAR CHS Correlation Plot All PFAS.jpeg"), 
       width = 16, 
       height = 7)
  
