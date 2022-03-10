#################################################
# Global temperature maps R script from Valerio #
#################################################

# set working directory
rm(list = ls())
getwd()
setwd("C:/Users/wcepv/surfdrive/papers in progress/pcrit fish/global maps")

library(foreach); library(dplyr); library(raster)

# get average water temperature for land from fishsuit

climate_models <- c('gfdl','ipsl','hadgem','miroc','noresm')

ls_files <- paste0(climate_models,'/Tma','_hist.tif')

r <- foreach(i = ls_files) %do% raster(i) %>% brick()

rav<-r
rav <- calc(r,mean,na.rm = T)

#convert from kelvin to celcius

rav <- rav-273.15

# there are < 0.02% of values with temperatures of 0 degrees or lower --> set them to NA

sum(rav[]<=0,na.rm=T)/sum(!is.na(rav[]))*100

rav[][rav[] <= 0] <- NA

# check

plot(rav)



# get marine water temperature data

mar <- raster('mar/longterm_avg_C.grd')

#check

mar[][mar[] <= 0] <- 10

mar2<-mar
mar2[][is.na(mar[])] <- 0
plot(mar2)


# read the model
m <- readRDS("C:/Users/wcepv/surfdrive/papers in progress/pcrit fish/best_model.rds")
data <- read.csv("C:/Users/wcepv/surfdrive/papers in progress/pcrit fish/Pcrit data for 195 fish species ( 20210823 ).csv")

summary(m)$estimates
fixedEstimates                       <-  brms::fixef(m, estimate = 'mean')
intercept                            <-  fixedEstimates['Intercept', 'Estimate']
c_value_log_Slope                    <-  fixedEstimates['c_value_log', 'Estimate']
temp_test_Slope                      <-  fixedEstimates['temp_test', 'Estimate']
mass_log_Slope                       <-  fixedEstimates['mass_log', 'Estimate']
sal_test_Slope                       <-  fixedEstimates['sal_test', 'Estimate']
res_metab_rate_Slope                 <-  fixedEstimates['res_metab_rate', 'Estimate']
c_value_log_and_temp_test_Slope      <-  fixedEstimates['c_value_log:temp_test', 'Estimate']
mass_log_and_temp_test_Slope         <-  fixedEstimates['temp_test:mass_log', 'Estimate']

summary(10^data$logMass)
summary(log10(data$Estimated_C.Value))
1 pg, 10 grams
4 pg, 1 kg
# use temperature data to calculate Pcrit for a 1 kg fish with 4pg genome size
mar_large<-(intercept+c_value_log_Slope*log10(4)+mass_log_Slope*log10(1000)+sal_test_Slope*35+
  temp_test_Slope*mar+
  c_value_log_and_temp_test_Slope*log10(4)*mar+
  mass_log_and_temp_test_Slope*log10(1000)*mar)
fw_large<-(intercept+c_value_log_Slope*log10(4)+mass_log_Slope*log10(1000)+
              temp_test_Slope*rav+
              c_value_log_and_temp_test_Slope*log10(4)*rav+
              mass_log_and_temp_test_Slope*log10(1000)*rav)

# use temperature data to calculate Pcrit for a 1 kg fish with 4pg genome size
mar_small<-(intercept+c_value_log_Slope*log10(0.5)+mass_log_Slope*log10(5)+sal_test_Slope*35+
              temp_test_Slope*mar+
              c_value_log_and_temp_test_Slope*log10(0.5)*mar+
              mass_log_and_temp_test_Slope*log10(5)*mar)
fw_small<-(intercept+c_value_log_Slope*log10(0.5)+mass_log_Slope*log10(5)+
              temp_test_Slope*rav+
              c_value_log_and_temp_test_Slope*log10(0.5)*rav+
              mass_log_and_temp_test_Slope*log10(5)*rav)


# merging the layers

# first make sure extent, projection and resolution are exactly the same
fw_small <- extend(fw_small,extent(mar_small))
res(fw_small) <- res(mar_small)
crs(fw_small) <- crs(mar_small)

fw_large <- extend(fw_large,extent(mar_large))
res(fw_large) <- res(mar_large)
crs(fw_large) <- crs(mar_large)

# now substitute cells in mar with rav
small <- mar_small
small[][!is.na(fw_small[])] <- fw_small[][!is.na(fw_small[])]

large <- mar_large
large[][!is.na(fw_large[])] <- fw_large[][!is.na(fw_large[])]

plot(small)
plot(large)
# there are values with Pcrit below 2 --> set them to 2
large[][large[] <= 2] <- 2
small[][small[] <= 2] <- 2

par(mfrow=c(2,1))

plot(21/small)
title("small fish with a small genome",cex=0.3)
plot(21/large)
title("large fish with a large genome",cex=0.3)

combo_lg<-21/large
combo_sm<-21/small


# plot

library(sf); library(ggplot2)

# base layers

# crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

crs_custom <- 4326

world <- rnaturalearth::ne_countries(returnclass = "sf")[,1] %>%
  
  st_transform(crs_custom)

bb <- rnaturalearth::ne_download(type = "wgs84_bounding_box", category = "physical",returnclass = "sf") %>%
  
  st_transform(crs_custom)

graticules <- rnaturalearth::ne_download(type = "graticules_30", category = "physical",returnclass = "sf") %>%
  
  st_transform(crs_custom)



# convert raster to data frame for ggplot

df_lg <- as(combo_lg %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
  
  as.data.frame(.) %>%
  
  dplyr::select(x,y,value = 'layer')

df_sm <- as(combo_sm %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
  
  as.data.frame(.) %>%
  
  dplyr::select(x,y,value = 'layer')


# and draw

p <- ggplot() +
  
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  
  # geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  
  geom_sf(data = world, fill = "grey90", lwd = NA) +
  
  geom_raster(data = df_sm, aes(x=x, y=y, fill=value)) +
  
  scale_fill_gradientn(
    
    colors = inlmisc::GetColors(length(seq(0,20,5)),scheme='inferno'),
    
    breaks = seq(0,20,2),
    
    labels = paste0(seq(0,20,2)),
    
    na.value = 'transparent') +
  
  theme_minimal() +
  
  theme(text = element_text(size = 12),
        
        panel.grid.major = element_line(color=NA),
        
        axis.text = element_blank(),
        
        axis.title = element_blank(),
        
        legend.position = 'bottom',
        
        legend.key.width = unit(3,'line')
        
      # legend.title = element_text('dd')
      #  legend.title = element_blank()
  )

sm<-p
sm

p <- ggplot() +
  
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  
  # geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  
  geom_sf(data = world, fill = "grey90", lwd = NA) +
  
  geom_raster(data = df_lg, aes(x=x, y=y, fill=value)) +
  
  scale_fill_gradientn(
    
    colors = inlmisc::GetColors(length(seq(0,20,5)),scheme='inferno'),
    
    breaks = seq(0,20,2),
    
    labels = paste0(seq(0,20,2)),
    
    na.value = 'transparent') +
  
  theme_minimal() +
  
  theme(text = element_text(size = 9),
        
        panel.grid.major = element_line(color=NA),
        
        axis.text = element_blank(),
        
        axis.title = element_blank(),
        
        legend.position = 'bottom',
        
        legend.key.width = unit(2,'line')
        
        # legend.title = element_text('dd')
        #  legend.title = element_blank()
  )

lg<-p
lg

sm<- sm+ labs(title = "Small fish with a small genome")+ theme(legend.position = "none")+
  theme(plot.title = element_text(color = "black", size = 12,hjust = 0.5))

lg<-lg + labs(title = "Large fish with a large genome", 
         fill = "Factorial Aerobic Scope \n")+
theme(plot.title = element_text(color = "black", size = 12,hjust = 0.5),
      legend.title = element_text(color = "black", size = 10,vjust = 0.5)
      )

sav<-ggarrange(sm, lg, heights = c(1,1))

library(grid)
ggsave(paste0('FAS_map.pdf'),sav,
       width = 130,height = 210,dpi = 600,units = 'mm')
ggsave(paste0('FAS_map.jpg'),sav,
       width = 130,height = 210,dpi = 600,units = 'mm')
