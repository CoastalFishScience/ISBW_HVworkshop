#' """ Puerto Rico coral benthic community trait diversity 
#'     @author: W. Ryan James
#'     date: 6/21/24"""

# load libraries
library(tidyverse)
library(hypervolume)
library(truncnorm)

# load trait data
df_tr = read_csv('data/coral_traits.csv') 
df_tr

# load benthic community percent cover data across regions
# pc = percent cover 0-100
# pc_sd = standard deviation in percent cover
df_ben = read_csv('data/coral_PC.csv')
df_ben

# randomly generate percent cover data for each species and region
# based on mean and sd of percent cover for the designated number of reps
reps = 100

set.seed(14)
df = df_ben |> 
  # join trait data to benthic data
  left_join(df_tr, by = 'species') |> 
  # duplicate the data to create the number of reps needed
  slice(rep(1:n(), each=reps))|> 
  # randomly generate percent cover
  mutate(i = rep(1:reps, times=nrow(df_ben)),
         percentcover = truncnorm::rtruncnorm(1, a = 0.001, b = 100,
                                              mean = pc, sd = pc_sd)) |> 
  select(region, `corallite diameter`:`percentcover`)

df 

# z-score across all regions for each rep and generate hypervolumes
df = df |> 
  # z-score across regions for each rep
  group_by(i) |> 
  mutate(across(`corallite diameter`:`colony maximum diameter`, scale)) |> 
  # nest data for each region and rep to make hypervolume
  group_by(region, i) |> 
  # create a column for the percent cover to weight hypervolume as well as input data
  nest(weight = percentcover, data = `corallite diameter`:`colony maximum diameter`) |> 
  # create community weighted hypervolumes 
  mutate(hv = map2(data,weight, \(data,weight) hypervolume_gaussian(data, 
                                                                    name = paste(region,i,sep = '_'),
                                                                    weight = weight$percentcover,
                                                                    samples.per.point = 1000,
                                                                    kde.bandwidth = estimate_bandwidth(data), 
                                                                    sd.count = 3, 
                                                                    quantile.requested = 0.95, 
                                                                    quantile.requested.type = "probability", 
                                                                    chunk.size = 1000, 
                                                                    verbose = F)),
  # extrace size for each hypervolume 
          hv_size = map_dbl(hv, \(hv) get_volume(hv)))

# do not try to open df it will freeze your r since it is too big
head(df)
# save output
saveRDS(df, 'data/coral_region_hvs.rds')

# plot hypervolume size
d = df |> 
  group_by(region) |> 
  summarize(mean = mean(hv_size),
            upper = quantile(hv_size, 0.975),
            lower = quantile(hv_size, 0.025)) |> 
  mutate(region = factor(region, 
                         levels = c('North/Northeast', 'Vieques/Culebra',
                                    'Southeast', 'South', 'Southwest',
                                     'West', 'Mona/Desecheo'),
                         labels = c('North/\nNortheast', 'Vieques/\nCulebra',
                                    'Southeast', 'South', 'Southwest',
                                    'West', 'Mona/\nDesecheo')))


ggplot(d, aes(region, mean, color = region))+
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 1.5, linewidth = 2)+
  labs(x = 'Region', y = 'Trait diversity', color = 'Region')+
  scale_y_log10()+
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

ggsave('coral_size.png', 
       units="in", width=10, height=6, dpi=600)
