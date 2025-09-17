# contour_maps_inflow.R ######
# Updated: 11/30/2023
# Catarina Pien and Lisa Elliott (USBR)
# cpien@usbr.gov; lelliott@usbr.gov

#  This code uses zone of influence modeling results (DSM2) to create contours showing
#  how zone of influence changes from operational facilities based on pumping
#  Contour lines of a specific level are then compared between different flow levels
#  to indicate how changing OMR will influence the zone of influence.

# general
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(knitr)
library(tidyr)
library(sf)
library(stars)
library(readr)
library(janitor)
library(here)

#IDW
library(sp)
library(gstat)
library(raster)
library(tmap)
#library(rgdal)

# Visualization
library(ggmap)
library(ggspatial)
library(deltamapr)
library(viridis)

# This file contains many of the functions used to create contour maps
source(here("functions_zoi.R"))

# Read/Join data ---------------------------------------------------------
delta <- st_read(here("shapefiles/Bay_Delta_Poly_New.shp"))
# ZOI files are created in a Jupyter Notebook script for proportional overlap:
# https://github.com/FlowWest/lto-proportion-overlap-calculations
zoi_file_NAA = list.files("data_raw/zoi/", pattern = "NAA_.*csv$", full.names = TRUE)
zoi_file_Act5 = list.files("data_raw/zoi/", pattern = "Act5.*csv$", full.names = TRUE)


zoi_data_NAA <- lapply(zoi_file_NAA, read_csv) %>%
  bind_rows(.id = "id") %>%
  mutate(id = paste0("-", substr(zoi_file_NAA[as.numeric(id)],26, 29))) %>%
  rename(OMR_Flow = id) %>%
  mutate(OMR_Flow = if_else(OMR_Flow == "-sTha", "<-5500", OMR_Flow))%>%
  mutate(Alt = "NAA")
zoi_data_Act5 <- lapply(zoi_file_Act5, read_csv) %>%
  bind_rows(.id = "id") %>%
  mutate(id = paste0("-", substr(zoi_file_Act5[as.numeric(id)],28, 31))) %>%
  rename(OMR_Flow = id) %>%
  mutate(OMR_Flow = if_else(OMR_Flow == "-1300", "<-5500", OMR_Flow)) %>%
  mutate(Alt = "Act5")

# combine each individual file
zoi_data <- rbind(zoi_data_NAA, zoi_data_Act5)

# read in nodes, channels data
nodes <- st_read("shapefiles/nodes.shp") %>%
  dplyr::select(node)
nodes_4326 <- st_transform(nodes, crs = 4326) %>%
  mutate(points = "DSM2 nodes")
channels0 <- read_csv("data_raw/DSM2_Version822_Grid_20231102.csv") %>% # associates the down/up nodes and channel length between
  janitor::clean_names()  %>%
  rename(channel_number = chan_no) %>%
  dplyr::select(-manning, -dispersion)

# Drop nodes that are causing issues
dropNodes <- c(146, 147, 148, 206, 242, 246, 432, 433, 434)
# Drop duplicate channels
channels1 <- channels0[!channels0$upnode %in% dropNodes, ]
channels <- channels1[!channels1$downnode %in% dropNodes, ]
total_channel_length <- sum(channels$length)

# Join channel lengths with zoi data
zoi_channel <- left_join(zoi_data, channels)
#zoi_channel <- merge(zoi_data, channels, by = "channel_number")

# Change projections to 4326 (WGS)
delta_4326 <- st_transform(delta, crs = 4326) %>%
  mutate(line = "analysis boundary")
nodes_4326 <- st_transform(nodes, crs = 4326) %>%
  mutate(points = "DSM2 nodes")
WW_Delta_4326 <- st_transform(WW_Delta, crs = st_crs(delta_4326))
WW_Delta_crop <- st_crop(WW_Delta_4326,xmin = -122.2, xmax = -121, ymin = 37.5, ymax = 38.8) %>%
  filter(HNAME!= "SAN FRANCISCO BAY")

# Convert data frame to long
zoi_channel_long <- zoi_channel %>%
  pivot_longer(cols = c(lolo:hihi), names_to = "group", values_to = "overlap")

# Write data for channel length script
write_csv(zoi_channel_long, "data_export/prop_overlap_data_long_act5.csv")

# Look at raw data -----------------------------------
summary_vals <- zoi_channel_long %>%
  filter(overlap>=0) %>%
  group_by(group, OMR_Flow, Alt) %>%
  summarize(min = min(overlap),
            max = max(overlap),
            mean = mean(overlap)) %>%
  ungroup()

## QC plots ----------------------
ggplot(summary_vals) +
  geom_col(aes(x = Alt, y = mean, fill = Alt)) + facet_grid(OMR_Flow~group)+
  scale_fill_viridis_d()
ggplot(zoi_channel_long %>% filter(overlap>=0)) +
  geom_jitter(aes(x = Alt, y = overlap, color = node)) + facet_wrap(~group)
ggplot(zoi_channel_long %>% filter(overlap>=0)) +
  geom_violin(aes(x = Alt, y = overlap, fill = Alt)) + facet_grid(OMR_Flow~group) +
  scale_fill_viridis_d()

# Create data frames for each inflow-OMR group --------------------
# Look at which combinations are missing for the next exercise
#png("figures/allalts_missingcombos.png", units = "in", width = 7, height = 7, res = 300)
zoi_channel_long %>% filter(overlap>=0) %>%
  group_by(Alt, group, OMR_Flow) %>%
  summarize(n = n()) %>%
  mutate(group = factor(group, levels = inflow_order)) %>%
  ggplot()  + geom_tile(aes(x = OMR_Flow, y = group, fill = n), color = "black") + facet_wrap(~Alt) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
#dev.off()

# filtered to medium hydrologic overlap sample sizes
#png("figures/allalts_samplesizes_medhydro.png", units = "in", width = 7, height = 7, res = 300)
zoi_channel_long %>% filter(overlap>=0 & overlap <=0.75) %>%
  group_by(Alt, group, OMR_Flow) %>%
  summarize(n = n()) %>%
  mutate(group = factor(group, levels = inflow_order)) %>%
  ggplot()  +
  geom_tile(aes(x = OMR_Flow, y = group, fill = n), color = "black") +
  geom_text(aes(x = OMR_Flow, y = group, label = n), color = "gray65", size = 2.7) +
  facet_wrap(~Alt) +
  viridis::scale_fill_viridis(option = "plasma") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
#dev.off()


# Run contour functions -------------------------------------
inflow_order = c("lolo", "lomed", "lohi", "medlo", "medmed", "medhi", "hilo", "himed", "hihi")
alt_order = c("NAA","Act5")
delta_sp <- as(delta_4326, "Spatial")

# NAA
lolo_contour_NAA <- f_data_interp_contour(gpname = "lolo", altname = "NAA")
lomed_contour_NAA <- f_data_interp_contour_no5500(gpname = "lomed", altname = "NAA")
lohi_contour_NAA <- f_data_interp_contour_no5500(gpname = "lohi", altname = "NAA")
medlo_contour_NAA <- f_data_interp_contour(gpname = "medlo", altname = "NAA")
medmed_contour_NAA <- f_data_interp_contour(gpname = "medmed", altname = "NAA")
medhi_contour_NAA <- f_data_interp_contour(gpname = "medhi", altname = "NAA")
hilo_contour_NAA <- f_data_interp_contour_no2000(gpname = "hilo", altname = "NAA")
himed_contour_NAA <- f_data_interp_contour(gpname = "himed", altname = "NAA")
hihi_contour_NAA <- f_data_interp_contour(gpname = "hihi", altname = "NAA")

# Act5

lolo_contour_Act5 <- f_data_interp_contour(gpname = "lolo", altname = "Act5")
lomed_contour_Act5 <- f_data_interp_contour_no5500(gpname = "lomed", altname = "Act5")
lohi_contour_Act5 <- f_data_interp_contour_no5500(gpname = "lohi", altname = "Act5")
medlo_contour_Act5 <- f_data_interp_contour(gpname = "medlo", altname = "Act5")
medmed_contour_Act5 <- f_data_interp_contour(gpname = "medmed", altname = "Act5")
medhi_contour_Act5 <- f_data_interp_contour(gpname = "medhi", altname = "Act5")
hilo_contour_Act5 <- f_data_interp_contour_no2000(gpname = "hilo", altname = "Act5")
himed_contour_Act5 <- f_data_interp_contour(gpname = "himed", altname = "Act5")
hihi_contour_Act5 <- f_data_interp_contour(gpname = "hihi", altname = "Act5")

## Combine contours --------------------------

# 0.75 represents contour at which 75% overlap exists
# combine first by alternative
contours_all_NAA <- rbind(lolo_contour_NAA, lomed_contour_NAA, lohi_contour_NAA,
                      medlo_contour_NAA, medmed_contour_NAA, medhi_contour_NAA,
                      hilo_contour_NAA, himed_contour_NAA, hihi_contour_NAA) %>%
  mutate(OMR_flow = factor(flow, levels = c("-2000", "-3500", "-5000", "<-5500")))

contours_all_Act5 <- rbind(lolo_contour_Act5, lomed_contour_Act5, lohi_contour_Act5,
                          medlo_contour_Act5,
                          medmed_contour_Act5,
                          medhi_contour_Act5,
                          hilo_contour_Act5, himed_contour_Act5,
                          hihi_contour_Act5
                          ) %>%
  mutate(OMR_flow = factor(flow, levels = c("-2000", "-3500", "-5000", "<-5500")))


save(contours_all_NAA, contours_all_Act5, file = "contours_act5.Rdata")

# Can start here if desired
# load("contours_allalts.Rdata")

# Make one contour file for all
alt_order2 = c("NAA","Act5")

contourGroup <- rbind(contours_all_NAA, contours_all_Act5)%>%
  mutate(grouper = paste0(group, "_", flow, group2),
         label = paste0(group2, "_", flow)) %>%
  rename(Inflow = group2) %>%
  mutate(Inflow = factor(Inflow, levels = inflow_order),
         Alt = factor(Alt, levels = alt_order2))

# Make contour plots -------------------------------------------
plot_contours(alt = "NAA", cont = 0.75)
plot_contours(alt = "Act5", cont = 0.75)


# Make output table of prop overlap by OMR for Turner Cut, SJR Jer --------

# Turner Cut, SJR at Jersey Point, Old R at Middle River
# Realized that these are by inflow group which is different than the table
nodesofinterest <- c(26, 469, 52)
nodenames <- read.csv("data_raw/inputFile_withnode.csv") %>%
  filter(node.no %in% nodesofinterest)

prop_overlap_nodesofinterest <- zoi_data |>
  filter(node %in% nodesofinterest)
