################################################ 
#### 1. Occurrence records
################################################ 

# Read ccurrence records
records <- read.csv(
  "gbif.epiita.csv",
  header = TRUE,
  sep = ";",
  dec = ".")

# Plot occurrence records
library(maptools)
data(wrld_simpl)
projection(wrld_simpl)<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
plot(wrld_simpl, axes=T)
points(records[,"lon"], records[,"lat"], pch= 20, col="blue", cex=1)
summary(records)


################################################ 
#### 2. Environmental variables
################################################ 

environmental.variables <- raster::stack("biooracle_current_variables.tif")
plot(environmental.variables)


studyarea200_raster <- raster::stack("studyarea200.tif")
e1<-extent(c(-60,-25, -40, 5)) #Atlantico Sul
studyarea200_raster_atlantic <- crop(studyarea200_raster,e1)
plot(studyarea200_raster_atlantic)

plot(studyarea200_raster_atlantic$sst_mean)

# Plot continent
e2<-extent(c(-90,-25, -40, 5)) #Atlantico Sul
atlantic <- crop(wrld_simpl,e2)
plot(atlantic, axes=T, add=T)

# Download states shapefile
bra_states <- geodata::gadm(
  country = "BRA",
  level = 1,
  path = tempdir())

# Add Brazilian states
plot(
  bra_states,
  add = TRUE,
  border = "black",
  lwd = 0.7)

################################################ 
#### 3. Ecological niche modelling
################################################ 


# Import model
model_ensemble_epiita <- raster::stack("model_epiita_tss.tif")


colfunc <- colorRampPalette(
  c("red4", "red", "yellow", "forestgreen"))

# Plot model
plot(
  model_ensemble_epiita,
  axes = TRUE,
  col = colfunc(20),
  zlim = c(0, 1),
  main = "Environmental suitability")

# Plot continent
plot(atlantic, axes=T, add=T)

# Add Brazilian states
plot(
  bra_states,
  add = TRUE,
  border = "black",
  lwd = 0.7)

# Add occurrence records
points(
  records$lon,
  records$lat,
  pch = 20,
  col = "blue",
  cex = 1.5)



################################################ 
#### 4.3. Variable importance
################################################ 


vi_model_var <- c(
  sst_mean = 0.3198467,
  sst_range = 0.2249733,
  ph_mean = 0.1754267,
  dissolved_oxygen = 0.1913767
)

# Ordenar do maior para o menor
vi_model_var <- sort(vi_model_var, decreasing = TRUE)
vi_model_var

library(ggplot2)

vi_df <- data.frame(
  variable = names(vi_model_var),
  importance = as.numeric(round(vi_model_var,2)))

# Barplot
ggplot(vi_df, aes(x = reorder(variable, importance),y = importance)) +
  geom_col(  fill = "#5B8DB8",
             width = 0.7  ) +
  coord_flip() +
  geom_text(aes(label = round(importance, 3)),
            hjust = -0.15,
            size = 4  ) +  labs(x = NULL, y = "Variable importance"  ) +
  ylim(0, max(vi_df$importance) * 1.15) +
  theme_minimal(base_size = 13) +
  theme(  panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank()  )



################################################ 
#### 4.4. Model evaluation
################################################ 

evaluation.model <- getEvaluation(
  model_epiita,
  wtest = "test.dep",
  stat = c("AUC", "TSS", "threshold", "deviance"),
  opt = 2)

# Data frame 
perf_df <- data.frame( metric = c( rep("AUC", length(evaluation.model$AUC)), 
                                   rep("TSS", length(evaluation.model$TSS)) ), 
                       value = c( evaluation.model$AUC, evaluation.model$TSS ) )
# Plot 
ggplot( perf_df, aes( x = metric, y = value ) ) +
  geom_boxplot( fill = "#5B8DB8", width = 0.5, outlier.alpha = 0.6 ) + 
  labs( x = NULL, y = "Model performance" ) + ylim(0, 1) + theme_minimal(base_size = 13) +
  theme( panel.grid.minor = element_blank() )


################################################ 
#### 4.6. Create binary map
################################################ 

# Arbitrary threshold
threshold_mean <- 0.25

# Binary model
binary.model <- model_ensemble_epiita >= threshold_mean

plot(
  binary.model,
  col = c("white", "red"),
  legend = FALSE,
  main = "Binary suitability map")

# Plot continent
e2<-extent(c(-90,-25, -40, 5)) #Atlantico Sul
atlantic <- crop(wrld_simpl,e2)
plot(atlantic, axes=T, add=T)

# Add Brazilian states
plot(
  bra_states,
  add = TRUE,
  border = "black",
  lwd = 0.7)

# Add occurrence records
points(
  data.epiita$lon,
  data.epiita$lat,
  pch = 20,
  col = "blue",
  cex = 1.5)


