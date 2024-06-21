#' """ FL Bay SAV state and stability using centroid distance
#'     @author: W. Ryan James, Rolando O. Santos, Benjamin Jones, Jennifer Rehage, Marianna Coppola, Gina Badlowski, Jonathan Rodemann
#'     date: 6/21/24"""

# load libraries
library(tidyverse)
library(hypervolume)

# load sav monitoring data 
# BASIN = Basin sampled
# YEAR = year of monitoring
# STATION = monitoring station 
# TT = Thalassia testudium percent cover
# HW = Halodule wrightii percent cover
# SF = Syringodium filiforme percent cover
# TMA = total macroalgae percent cover
# TDR = total drift algae percent cover
# sg_rich = seagrass species richness

df = read_csv('data/FLbay_SAV.csv') 
head(df)

# z-score and nest data to make hypervolume
set.seed(14)
df = df |> 
  # z score data across all sites and years
  mutate(across(c(TT:sg_rich), scale), 
  # add tiny amount so when all values the same can make hv       
         across(c(TT:sg_rich), 
                ~map_dbl(., ~. + rnorm(1, mean = 0, sd = 0.0001)))) |> 
  # remove station from dataset
  select(-STATION) |> 
  # nest data by basin and year
  group_by(BASIN, YEAR) |> 
  nest() 
head(df)

# generate hypervolumes
df = df |> 
  mutate(hv = map(data, \(data) hypervolume_gaussian(data, name = paste(BASIN,YEAR,sep = '_'),
                                                     samples.per.point = 1000,
                                                     kde.bandwidth = estimate_bandwidth(data), 
                                                     sd.count = 3, 
                                                     quantile.requested = 0.95, 
                                                     quantile.requested.type = "probability", 
                                                     chunk.size = 1000, 
                                                     verbose = F)),
         centroid = map(hv, \(hv) get_centroid(hv)))
# do not try to open dataframe with hv column it will hang because it is too big
head(df) 

# save output as .rds 
#saveRDS(df, 'data/SAV_hvs.rds') 
# df = readRDS('data/hvAll.rds')

# plot hypervolumes
hvj = hypervolume_join(df$hv[[1]], df$hv[[2]])

plot(hvj, pairplot = T, colors=c('goldenrod','blue'),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE)

# comparison of across each year
df_y= tibble(y1 = unique(df$YEAR),
             y2 = unique(df$YEAR)) |> 
  expand(y1,y2)

# make all unique year comparisons 
df_y = df_y[!duplicated(t(apply(df_y,1,sort))),] %>% 
  filter(!(y1 == y2))

# make two df to join all unique comparisons  
df1 = df |> 
  select(BASIN, y1 = YEAR, hv1 = hv, cent1 = centroid)

df2 = df |> 
  select(BASIN, y2 = YEAR, hv2 = hv, cent2 = centroid)


# create data frame of all data and make yearly comparisons
df_cd = tibble(BASIN = rep(unique(df$BASIN),
                           each = nrow(df_y)),
               y1 = rep(df_y$y1, times = length(unique(df$BASIN))),
               y2 = rep(df_y$y2, times = length(unique(df$BASIN)))) |> 
  inner_join(df1, by = c('BASIN', 'y1')) |> 
  inner_join(df2, by = c('BASIN', 'y2')) |> 
  mutate(ychange = y2-y1,
  # join hypervolumees in a set for centroid distance
         set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
  # calculate centroid distance 
         dist_cent = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F)),
  # calculate the difference of centroid of each axis
         dif = map2(cent1, cent2, \(cent1,cent2) cent2 - cent1)) |> 
  #unnest centroid differences
  unnest_wider(dif) |> 
  # select only metrics of interest
  select(BASIN, y1, y2, ychange,  
         dist_cent, TT, HW, SF, sg_rich, TMA, TDR)
df_cd

# save output
write_csv(df_cd, 'data/SAV_centDist.csv')

# plot centroid distance
df_cd = read_csv('data/SAV_centDist.csv') |> 
  mutate(BASIN = factor(BASIN, levels = 
                          c('JON', 'RKB', 'RAN', 'EAG')))

ggplot(df_cd, aes(BASIN, dist_cent, fill = BASIN))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed', linewidth = 1)+
  geom_point(aes(color = BASIN), size = 1, 
             position=position_jitterdodge(dodge.width = 1, jitter.width = 1))+
  # geom_errorbar(aes(ymin = lc, ymax = uc), linewidth = 2, width = 0)+
  geom_boxplot(alpha = 0.6, outliers = F)+
  labs(x = 'Basin', y = 'Centroid distance')+
  scale_fill_viridis_d(option = 'turbo')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# ggsave('boxCentDist.png',
#        units="in", width=8, height=5, dpi=600)

# trends in stability 
library(MuMIn)
df_cd = df_cd |> 
  group_by(BASIN) |>
  nest() |> 
  # fit intercept, linear, and quadratic model
  mutate(m_int = map(data, \(df)lm(dist_cent~1, data = df)),
         m_lin = map(data, \(df)lm(dist_cent~ychange, data = df)),
         m_quad = map(data, \(df)lm(dist_cent~ychange + I(ychange^2), data = df)),
         AICc_int = map_dbl(m_int, \(x) AICc(x)),
         AICc_lin = map_dbl(m_lin, \(x) AICc(x)),
         AICc_quad = map_dbl(m_quad, \(x) AICc(x)),
         model = case_when(
           AICc_int - min(c(AICc_int,AICc_lin,AICc_quad)) <= 4 ~ 'Intercept',
           AICc_lin < AICc_quad ~ 'Linear',
           AICc_quad < AICc_lin ~ 'Quadratic'))

# unnest data 
d = df_cd |> 
  select(BASIN, data, model) |> 
  unnest(cols = c(data)) |>  
  mutate(BASIN = factor(BASIN, levels = 
                          c('JON', 'RKB', 'RAN', 'EAG')))

ggplot(d, aes(ychange, dist_cent, color = BASIN))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(data = d |> filter(model == 'Intercept'),
              method = 'lm', formula = y~1, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Linear'),
              method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Quadratic'),
              method = 'lm', formula = y~x+I(x^2), 
              linewidth = 1, color = 'black')+
  facet_wrap(~BASIN,  nrow = 1)+
  labs(x = 'Years between comparison', y = 'Centroid distance')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# ggsave('stabCentDist.png',
#        units="in", width=8, height=5, dpi=600)

# variable importance 
# pivot data longer 
df_c = read_csv('data/SAV_centDist.csv') |> 
  mutate(across(TT:TDR, \(x) x^2)) |> 
  pivot_longer(TT:TDR, names_to = 'axis', values_to = 'dist')

# create vector of unique axes
ax = unique(df_c$axis)

# for loop to calculate variable importance of each axis
for(i in 1:length(ax)){
  d = df_c |> 
    # remove axis 
    filter(axis != ax[i]) |> 
    group_by(BASIN,y1,y2,ychange,dist_cent) |> 
    #calculate euclidean distance without axis 
    summarise(cd = sqrt(sum(dist))) |> 
    # calculate importance of axis
    mutate(imp = (dist_cent/cd) - 1,
           axis = ax[i])
  
  # bind data into new data frame to store 
  if(i == 1){
    df_imp = d
  }else{
    df_imp = bind_rows(df_imp, d)
  }
}

# calculate relative importance across all years for each basin
df_cdi = df_imp |> 
  #filter(ychange == 1) |>
  group_by(BASIN,axis) |>
  summarize(imp = mean(imp)) |>
  group_by(BASIN) |> 
  mutate(s_imp = imp/max(imp))|> 
  mutate(BASIN = factor(BASIN, levels = 
                          c('JON', 'RKB', 'RAN', 'EAG')),
         axis = factor(axis, levels = c('TDR', 'TMA', 'sg_rich',
                                        'SF', 'HW', 'TT')))
# labeler function for plotting
y_label_formatter = function(x) {
  ifelse(x %% 1 == 0, formatC(x, format = "f", digits = 0), formatC(x, format = "f", digits = 2))
}

ggplot(df_cdi, aes(axis, s_imp, fill = BASIN))+
  geom_col()+
  labs(x = 'Variable', y = 'Centroid distance variable importance')+
  coord_flip()+
  theme_bw()+
  facet_wrap(~BASIN,  nrow = 2)+
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter) +
  scale_x_discrete(labels = c('Total drift algae', 'Total Macroalgae', 'Seagrass richness',
                              expression(italic('Syringodium filiforme')), 
                              expression(italic('Halodule wrightii')),
                              expression(italic('Thalassia testudium'))))+
  scale_fill_viridis_d(option = 'turbo')+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('s_impCentDist.png', 
       units="in", width=10, height=6, dpi=600)
