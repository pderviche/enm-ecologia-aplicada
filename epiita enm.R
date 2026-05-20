################################################ 
#### ECOLOGIA APLICADA
#### Modelos de Nicho Ecológico e Distribuição de Espécies
#### Exemplo de caso com o mero
#### Professor Patrick Derviche
################################################ 

# Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/LEC/Aulas/ENMs/Exemplo mero")

# Packages
# install.packages("pacman")
if(!require(pacman)) install.packages("pacman")

pacman::p_load(
  terra, raster, corrplot, sp, dismo, kernlab, rJava,
  rgbif, sdm, biooracler, usdm)


################################################ 
#### 1. Occurrence records
################################################ 

# 1.1. Import data from GBIF

# Species Epinephelus itajara
gbif.epiita <- gbif("Epinephelus", species = "itajara")

# Select records by year
gbif.epiita <- gbif.epiita[gbif.epiita$year > 1950, ]

# Select columns 
gbif.epiita <- gbif.epiita[, c("lon", "lat", "species", "datasetKey")]

# Eliminate NAs
gbif.epiita <- na.omit(gbif.epiita)

# Check for duplicate coordinates
dups2 <- duplicated(gbif.epiita[, c("lon", "lat")])
sum(dups2)

gbif.epiita <- gbif.epiita[!dups2, ]

summary(gbif.epiita)

# Add reference column
reference <- rep("GBIF", nrow(gbif.epiita))
gbif.epiita <- cbind(gbif.epiita, reference)

summary(gbif.epiita)

# Save occurrence records
write.table(
  gbif.epiita,
  "gbif.epiita.csv",
  sep = ";",
  dec = ".",
  row.names = FALSE)

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

# 2.1. Download Bio-ORACLE v3 variables

constraints <- list(
  latitude  = c(-35, 45),
  longitude = c(-125, 145))

# SST mean and SST range
sst <- download_layers(
  dataset_id = "thetao_baseline_2000_2019_depthsurf",
  variables = c("thetao_mean", "thetao_range"),
  constraints = constraints)

# pH
ph <- download_layers(
  dataset_id = "ph_baseline_2000_2018_depthsurf",
  variables = "ph_mean",
  constraints = constraints)

# Dissolved oxygen
dissox <- download_layers(
  dataset_id = "o2_baseline_2000_2018_depthsurf",
  variables = "o2_mean",
  constraints = constraints)

# Bathymetry
bathy <- download_layers(
  dataset_id = "terrain_characteristics",
  variables = "bathymetry_mean",
  constraints = constraints)


################################################ 
#### 2.2. Calculate current mean
################################################ 

# Check layer names
names(sst)
names(ph)
names(dissox)
names(bathy)

# SST mean: average temporal layers
sst_mean_current <- mean(
  sst[[grep("thetao_mean", names(sst))]],
  na.rm = TRUE)

# SST range: average temporal layers
sst_range_current <- mean(
  sst[[grep("thetao_range", names(sst))]],
  na.rm = TRUE)

# pH: average temporal layers
ph_current <- mean(ph, na.rm = TRUE)

# Dissolved oxygen: average temporal layers
dissox_current <- mean(dissox, na.rm = TRUE)

# Bathymetry: static layer
bathy_current <- bathy


################################################ 
#### 2.3. Combine environmental variables
################################################ 

environmental.variables <- c(
  sst_mean_current,
  sst_range_current,
  ph_current,
  dissox_current,
  bathy_current)

names(environmental.variables) <- c(
  "sst_mean",
  "sst_range",
  "ph_mean",
  "dissolved_oxygen",
  "bathymetry_mean")

environmental.variables <- raster::stack("biooracle_current_variables.tif")
plot(environmental.variables)


################################################ 
#### 2.4. Save final current variables
################################################ 

writeRaster(
  environmental.variables,
  "biooracle_current_variables.tif",
  overwrite = TRUE)


################################################ 
#### 2.5. Pearson correlation
################################################ 

env_values <- values(environmental.variables, na.rm = TRUE)

cor_matrix <- cor(
  env_values,
  method = "pearson",
  use = "complete.obs")

cor_matrix

# Plot Pearson correlation using absolute values for colors
cor_abs <- abs(cor_matrix)

corrplot(
  cor_abs,
  method = "color",
  type = "upper",
  col = colorRampPalette(c(
    "#1A9850",
    "#FFFFBF",
    "#B2182B"
  ))(200),
  is.corr = FALSE,
  col.lim = c(0, 1),
  addCoef.col = "black",
  tl.cex = 0.9,
  number.cex = 0.8,
  tl.col = "black",
  number.digits = 2)


################################################ 
#### 2.6. Variance Inflation Factor (VIF)
################################################ 

# Convert SpatRaster to RasterStack for usdm

env <- environmental.variables
env_df <- terra::as.data.frame(env, xy = FALSE, na.rm = TRUE)
env_df <- env_df %>%
  dplyr::select(where(is.numeric))
env_df <- env_df[complete.cases(env_df), ]
env_df <- env_df[apply(env_df, 1, function(x) all(is.finite(x))), ]
set.seed(123)
if (nrow(env_df) > 10000) {
  env_df <- env_df[sample(seq_len(nrow(env_df)), 10000), ]}

# Calculate VIF
vif_result <- usdm::vif(env_df)
vif_result

# Save VIF results
write.table(
  vif_result,
  "vif.csv",
  sep = ";",
  dec = ".",
  row.names = FALSE)


################################################ 
#### 2.7. Extract values from each pixel
################################################ 

# Randomly sample pixels from raster stack
env.var.table <- terra::spatSample(
  environmental.variables,
  size = 1000000,      # number of sampled pixels
  method = "random",
  na.rm = TRUE,
  xy = TRUE,
  as.df = TRUE)

# Rename coordinate columns
colnames(env.var.table)[1:2] <- c("lon", "lat")

# Check number of valid sampled pixels
nrow(env.var.table)

# Save as CSV
write.table(
  env.var.table,
  "env.var_sample.csv",
  sep = ";",
  dec = ".",
  row.names = FALSE)


################################################ 
#### 2.8. Extract environmental values at occurrence points
################################################ 

# Occurrence records
records <- read.csv(
  "gbif.epiita.csv",
  header = TRUE,
  sep = ";",
  dec = ".")

# Check columns
head(records)

# Convert records to spatial vector
records_vect <- vect(
  records,
  geom = c("lon", "lat"),
  crs = crs(environmental.variables))

# Extract environmental values at occurrence points
values.records <- extract(
  environmental.variables,
  records_vect,
  cells = TRUE)

# Combine records and environmental values
data.epiita <- cbind(records, values.records)

# Remove duplicated cells
data.epiita <- data.epiita[!duplicated(data.epiita$cell), ]

# Remove NAs
data.epiita <- na.omit(data.epiita)

# Add presence column
data.epiita$sp <- 1

# Check final occurrence data
dim(data.epiita)
head(data.epiita)

# Save final occurrence data
write.table(
  data.epiita,
  "data.epiita.csv",
  sep = ";",
  dec = ".",
  row.names = FALSE)



################################################ 
#### 2.9. Plot environmental variables and occurrence records
################################################ 

plot(environmental.variables[["bathymetry_mean"]])

points(
  data.epiita$lon,
  data.epiita$lat,
  pch = 20,
  col = "blue",
  cex = 1.5)


################################################ 
#### 2.10. Clip bathymetry shallower than 200 m
################################################ 

# Isolate bathymetry layer
depth <- environmental.variables[["bathymetry_mean"]]

# Keep only marine areas shallower than 200 m
depth_200 <- depth

depth_200[depth_200 < -200] <- NA
depth_200[depth_200 > 0] <- NA

plot(depth_200)

points(
  data.epiita$lon,
  data.epiita$lat,
  pch = 20,
  col = "blue",
  cex = 2)

# Mask all environmental variables by bathymetry
studyarea200 <- mask(environmental.variables, depth_200)

plot(studyarea200)

# Save study area
writeRaster(
  studyarea200,
  "studyarea200.tif",
  overwrite = TRUE)


################################################ 
#### 3. Ecological niche modelling
################################################ 

# My memory use
gc()


################################################ 
#### 3.1. Import files
################################################ 

# Import raster as RasterStack for sdm package
studyarea200_raster <- raster::stack("studyarea200.tif")


e1<-extent(c(-60,-25, -40, 5)) #Atlantico Sul
studyarea200_raster_atlantic <- crop(studyarea200_raster,e1)

# Import occurrence data
data.epiita <- read.table(
  "data.epiita.csv",
  header = TRUE,
  sep = ";",
  dec = ".")

# Criar objeto espacial com coordenadas
coordinates(data.epiita) <- ~ lon + lat

proj4string(data.epiita) <- CRS(
  "+proj=longlat +datum=WGS84 +no_defs"
)

# Clear previous plots
dev.off()

# Plot occurrence records
plot(studyarea200_raster_atlantic[["bathymetry_mean"]])

e2<-extent(c(-90,-25, -40, 5)) #Atlantico Sul
atlantic <- crop(wrld_simpl,e2)
plot(atlantic, axes=T, add=T)

points(
  data.epiita$lon,
  data.epiita$lat,
  pch = 20,
  col = "blue",
  cex = 2)

summary(studyarea200)


################################################
#### 3.2. Preparing data and generating background
################################################

pre_model <- sdmData(
  formula = sp ~ sst_mean + sst_range + ph_mean + dissolved_oxygen,
  train = data.epiita,
  predictors = studyarea200_raster_atlantic,
  bg = list(
    n = 1000,
    remove = TRUE))


################################################
#### 3.3. Adjusting the models
################################################

model_epiita <- sdm(
  formula = sp ~ sst_mean + sst_range + ph_mean + dissolved_oxygen,
  data = pre_model,
  methods = c("maxent", "rf", "gam"),
  replication = "boot",
  test.percent = 30,
  n = 20)


################################################
#### 3.4. Ensemble
################################################

model_ensemble_epiita <- ensemble(
  model_epiita,
  newdata = studyarea200_raster_atlantic,
  filename = "model_epiita_tss.tif",
  setting = list(
    method = "weighted",
    stat = "TSS",
    opt = 2),
  overwrite = TRUE)


################################################
#### 3.5. Plot ensemble model
################################################

library(geodata)

# Import model
model_ensemble_epiita <- raster::stack("model_epiita_tss.tif")

# Download states shapefile
bra_states <- geodata::gadm(
  country = "BRA",
  level = 1,
  path = tempdir())

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
  data.epiita$lon,
  data.epiita$lat,
  pch = 20,
  col = "blue",
  cex = 1.5)

################################################ 
#### 4. Assessing models
################################################ 

################################################ 
#### 4.1. Model information
################################################ 

# Model information
model.info <- getModelInfo(model_epiita)

# Check successful model runs
model.info

# Save model information
write.table(
  model.info,
  "model_info_epiita.csv",
  sep = ";",
  dec = ".",
  row.names = FALSE
)


################################################ 
#### 4.2. Plot ensemble model
################################################ 

colfunc <- colorRampPalette(
  c("red4", "red", "yellow", "forestgreen"))

plot(
  model_ensemble_epiita,
  axes = TRUE,
  col = colfunc(20),
  zlim = c(0, 1),
  main = "Environmental suitability")

points(
  data.epiita$lon,
  data.epiita$lat,
  pch = 20,
  col = "blue",
  cex = 0.8)


################################################ 
#### 4.3. Variable importance
################################################ 

fun_model <- function(x) {
  vi_model <- getVarImp(
    model_epiita,
    id = x,
    wtest = "test.dep")
    return(as.data.frame(vi_model@varImportance[[3]]))}

# Model IDs
V_model <- model.info$modelID

# Variable importance for each model
M_model <- lapply(V_model, fun_model)

# Combine results
vi_model <- do.call(cbind, M_model)

# Transpose
vi_model <- t(vi_model)

# Rename columns
colnames(vi_model) <- c(
  "sst_mean",
  "sst_range",
  "ph_mean",
  "dissolved_oxygen"
)

# Mean importance
vi_model_var <- colMeans(
  vi_model,
  na.rm = TRUE
)

vi_model
vi_model_var

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

evaluation.model

summary(evaluation.model)

# Save evaluation
write.table(
  evaluation.model,
  "evaluation_model_epiita.csv",
  sep = ";",
  dec = ".",
  row.names = FALSE)

# Plot model performance
boxplot(
  evaluation.model$AUC,
  evaluation.model$TSS,
  names = c("AUC", "TSS"),
  ylab = "Model performance",
  ylim = c(0, 1))

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
#### 4.5. Predicted probability of occurrence
################################################ 

# Occurrence records
records <- read.csv(
  "gbif.epiita.csv",
  header = TRUE,
  sep = ";",
  dec = ".")

# Keep important columns
records <- records[, c("lon", "lat", "species", "reference")]
# Remove records without coordinates
records <- records[complete.cases(records[, c("lon", "lat")]), ]

# Convert predictors to SpatRaster
studyarea200_terra <- terra::rast(studyarea200_raster_atlantic)

# Convert occurrence records to SpatVector
records_vect <- terra::vect(
  records,
  geom = c("lon", "lat"),
  crs = "EPSG:4326")

# Reproject occurrence records to raster CRS, if needed
records_vect <- terra::project(
  records_vect,
  terra::crs(studyarea200_terra))

# Extract environmental variables
values.records <- terra::extract(
  studyarea200_terra,
  records_vect,
  cells = TRUE)

# Remove ID column
values.records <- values.records[, -1]

# Combine occurrence records and environmental values
data.epiita <- cbind(records, values.records)

# Remove duplicated raster cells
data.epiita <- data.epiita[!duplicated(data.epiita$cell), ]

# Remove NAs
data.epiita <- na.omit(data.epiita)

# Presence column
data.epiita$sp <- 1

# Check names
names(data.epiita)


#### Extract predicted suitability at occurrence points

# Convert ensemble model to SpatRaster, if needed
model_ensemble_terra <- terra::rast(model_ensemble_epiita)

# Convert cleaned occurrence records to SpatVector
data.epiita_vect <- terra::vect(
  data.epiita,
  geom = c("lon", "lat"),
  crs = "EPSG:4326")

# Reproject occurrence records to ensemble raster CRS, if needed
data.epiita_vect <- terra::project(
  data.epiita_vect,
  terra::crs(model_ensemble_terra))

# Extract predicted suitability
pred_values <- terra::extract(
  model_ensemble_terra,
  data.epiita_vect,
  cells = TRUE)

# Remove ID column
pred_values <- pred_values[, -1]

# Combine coordinates and predicted suitability
predicted.probability.of.occurrence <- cbind(
  data.epiita[, c("lon", "lat", "species", "reference")],
  pred_values)

# Remove NAs
predicted.probability.of.occurrence <- na.omit(
  predicted.probability.of.occurrence)

# Check result
predicted.probability.of.occurrence

# Save predicted suitability at occurrence points
write.table(
  predicted.probability.of.occurrence,
  "predicted_probability_of_occurrence_epiita.csv",
  sep = ";",
  dec = ".",
  row.names = FALSE)

# Boxplot
boxplot(predicted.probability.of.occurrence$ensemble_weighted,
  xlab = "Occurrence records",
  ylab = "Environmental suitability",
  ylim = c(0, 1))

summary(predicted.probability.of.occurrence$ensemble_weighted)
sd(predicted.probability.of.occurrence$ensemble_weighted)

# Data frame 
suitability_df <- data.frame( suitability = predicted.probability.of.occurrence$ensemble_weighted ) 

# Plot 
ggplot( suitability_df, aes( x = "Occurrence records", y = suitability ) ) + 
  geom_boxplot( fill = "#5B8DB8", width = 0.4, outlier.alpha = 0.6 ) + 
  labs( x = NULL, y = "Environmental suitability" ) + ylim(0, 1) + theme_minimal(base_size = 13) +
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


########## END