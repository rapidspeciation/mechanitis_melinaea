##### Distribution graphs #####

##### load libraries #####
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(sf)
library(tmap)
library(tmaptools)
library(mapview)

sf_use_s2(TRUE)

##### Import data #####
Mael <- get(load("~/Ithomiini_final.RData"))
OurOwn <- read_excel("~/coordinate list melanitis.xlsx", sheet = "Sheet1")
OurOwn <- na.omit(OurOwn)
dat_sf_mech <- st_as_sf(OurOwn, coords = c("Longitude", "Latitude"), crs = 4326)
dat_mael <- st_as_sf(Mael, coords = c("Longitude", "Latitude"), crs = 4326)

##Mael data of our subspecies:
Mechanitis <- Ithomiini_final[Ithomiini_final$Genus=="Mechanitis",]
Melinaea <- Ithomiini_final[Ithomiini_final$Genus=="Melinaea",]
unique(Mechanitis$Subsp_full)
unique(Melinaea$Subsp_full)
dat_mael_mech <- st_as_sf(Mechanitis, coords = c("Longitude", "Latitude"), crs = 4326)
dat_mael_mel <- st_as_sf(Melinaea, coords = c("Longitude", "Latitude"), crs = 4326)

##### Phylogeny paper #####
tmap_mode("view")

## east = red
## Mazaeus
conv_maz <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "mazaeus")))), dist = 15000), dTolerance = 5000)
conv_mazM <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.mazaeus.fallax")))), dist = 15000), dTolerance = 5000)

maz <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_mazM)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "mazaeus"))+
  tm_dots(size = 0.06, col = "red", shape = 21)
lf_maz <- tmap_leaflet(maz)
mapshot(lf_maz, file="~/phylodist-maz2.jpg")

## Messenoides
conv_mess <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "messenoides")))), dist = 15000), dTolerance = 5000)
conv_mesM <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.messenoides.messenoides" | 
                                                                                Subsp_full == "Mechanitis.messenoides.deceptus")))), dist = 15000), dTolerance = 5000)

mess <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_mesM)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "messenoides"))+
  tm_dots(size = 0.06, col = "red", shape = 21)
lf_mess <- tmap_leaflet(mess)
mapshot(lf_mess, file="~/phylodist-mess.jpg")

## Lysimnia Ecuador
conv_lysE <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "lysimniaEcu")))), dist = 15000), dTolerance = 5000)
conv_lysE.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.menapis.roqueensis")))), dist = 15000), dTolerance = 5000)

lysE <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_lysE.M)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "lysimniaEcu"))+
  tm_dots(size = 0.06, col = "red", shape = 21)
lf_lysE <- tmap_leaflet(lysE)
mapshot(lf_lysE, file="~/phylodist-lysE.jpg")

## Polymnia East
conv_polE <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "polymniaEast")))), dist = 15000), dTolerance = 5000)
conv_polE.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.mazaeus.dorissides")))), dist = 15000), dTolerance = 5000)

polE <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_polE.M)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "polymniaEast"))+
  tm_dots(size = 0.06, col = "red", shape = 21)

lf_polE <- tmap_leaflet(polE)
mapshot(lf_polE, file="~/phylodist-polE.jpg")

## Polymnia Peru
conv_polP <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "polymniaPeru")))), dist = 15000), dTolerance = 5000)
conv_polP.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.mazaeus.proceriformis" |
                                                                                  Subsp_full == "Mechanitis.mazaeus.proceriformes")))), dist = 15000), dTolerance = 5000)

polP <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_polP.M)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "polymniaPeru"))+
  tm_dots(size = 0.06, col = "red", shape = 21)

lf_polP <- tmap_leaflet(polP)
mapshot(lf_polP, file="~/phylodist-polP.jpg")

## west = yellow
## Macrinus
conv_mac <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "macrinus")))), dist = 15000), dTolerance = 5000)
conv_macM <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.macrinus." | 
                                                                                Subsp_full == "Mechanitis.lysimnia.utemaia")))), dist = 15000), dTolerance = 5000)

mac <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_macM)+
  tm_fill(col = "yellow",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "macrinus"))+
  tm_dots(size = 0.06, col = "darkorange", shape = 21)

lf_mac <- tmap_leaflet(mac)
mapshot(lf_mac, file="~/phylodist-mac.jpg")

## Menapis
conv_men <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "menapis")))), dist = 15000), dTolerance = 5000)
conv_menM <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.menapis.caribensis" | 
                                                                                Subsp_full == "Mechanitis.menapis.dariensis" |
                                                                                Subsp_full == "Mechanitis.menapis.mantineus" |
                                                                                Subsp_full == "Mechanitis.menapis.menapis" |
                                                                                Subsp_full == "Mechanitis.menapis.occasiva" |
                                                                                Subsp_full == "Mechanitis.menapis.saturata")))), dist = 15000), dTolerance = 5000)

men <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_menM)+
  tm_fill(col = "yellow",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "menapis"))+
  tm_dots(size = 0.06, col = "darkorange", shape = 21)
lf_men <- tmap_leaflet(men)
mapshot(lf_men, file="~/phylodist-men.jpg")

## Polymnia West
conv_polW <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "polymniaWest")))), dist = 15000), dTolerance = 5000)
conv_polW.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.polymnia.caucaensis" |
                                                                                  Subsp_full == "Mechanitis.polymnia.chimborazona" |
                                                                                  Subsp_full == "Mechanitis.polymnia.isthmia")))), dist = 15000), dTolerance = 5000)

polW <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_polW.M)+
  tm_fill(col = "yellow",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "polymniaWest"))+
  tm_dots(size = 0.06, col = "darkorange", shape = 21)

lf_polW <- tmap_leaflet(polW)
mapshot(lf_polW, file="~/phylodist-polW.jpg")

## atlantic forest = blue
## Lysimnia Brazil
conv_lysB <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "lysimniaBraz")))), dist = 15000), dTolerance = 5000)
conv_lysB.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.lysimnia.lysimnia")))), dist = 15000), dTolerance = 5000)

lysB <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_lysB.M)+
  tm_fill(col = "blue",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "lysimniaBraz"))+
  tm_dots(size = 0.06, col = "blue", shape = 21)

lf_lysB <- tmap_leaflet(lysB)
mapshot(lf_lysB, file="~/phylodist-lysB.jpg")

## Nesaea
conv_nes <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "nesaea")))), dist = 15000), dTolerance = 5000)
conv_nesM <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.lysimnia.nesaea")))), dist = 15000), dTolerance = 5000)

nes <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_nesM)+
  tm_fill(col = "blue",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "nesaea"))+
  tm_dots(size = 0.06, col = "blue", shape = 21)

lf_nes <- tmap_leaflet(nes)
mapshot(lf_nes, file="~/phylodist-nes.jpg")

## Polymnia Atlantic forest
conv_polB <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "polymniaBraz")))), dist = 15000), dTolerance = 5000)
conv_polB.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.polymnia.casabranca")))), dist = 15000), dTolerance = 5000)

polA <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_polB.M)+
  tm_fill(col = "blue",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "polymniaBraz"))+
  tm_dots(size = 0.06, col = "blue", shape = 21)

lf_polA <- tmap_leaflet(polA)
mapshot(lf_polA, file="~/phylodist-polA.jpg")

## guyanan shield = pink
## Polymnia Guaianan Shield
conv_polG <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "polymniaRest")))), dist = 15000), dTolerance = 5000)
conv_polG.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mech, Subsp_full == "Mechanitis.polymnia.polymnia")))), dist = 15000), dTolerance = 5000)

polG <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_polG.M)+
  tm_fill(col = "pink2",alpha=0.8)+
  tm_shape(filter(dat_sf_mech, Category == "polymniaRest"))+
  tm_dots(size = 0.06, col = "pink2", shape = 21)

lf_polG <- tmap_leaflet(polG)
mapshot(lf_polG, file="~/phylodist-polG.jpg")

### Melinaea ####
## east = red
## maeonis = east
conv_mae <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "maeonis")))), dist = 15000), dTolerance = 5000)
conv_mae.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Subsp_full == "Melinaea.satevis.maeonis")))), dist = 15000), dTolerance = 5000)

maeo <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_mae.M)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "maeonis"))+
  tm_dots(size = 0.06, col = "red", shape = 21)
lf_maeo <- tmap_leaflet(maeo)
mapshot(lf_maeo, file="~/phylodist-maeo.jpg")

## isocomma = east
conv_iso <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "isocomma")))), dist = 15000), dTolerance = 5000)
conv_iso.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Species == "isocomma")))), dist = 15000), dTolerance = 5000)

iso <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_iso.M)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "isocomma"))+
  tm_dots(size = 0.06, col = "red", shape = 21)

lf_iso <- tmap_leaflet(iso)
mapshot(lf_iso, file="~/phylodist-iso.jpg")

## marsaeus = east
conv_mars <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "marsaeus")))), dist = 15000), dTolerance = 5000)
conv_mars.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Subsp_full == "Melinaea.marsaeus.clara" |
                                                                                  Subsp_full == "Melinaea.marsaeus.macaria" |
                                                                                  Subsp_full == "Melinaea.marsaeus.rileyi" |
                                                                                  Subsp_full == "Melinaea.marsaeus.phasiana")))), dist = 15000), dTolerance = 5000)

mars <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_mars.M)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "marsaeus"))+
  tm_dots(size = 0.06, col = "red", shape = 21)

lf_mars <- tmap_leaflet(mars)
mapshot(lf_mars, file="~/phylodist-mars.jpg")

## Menophilus = east (Ecuador, Peru) (Menophilus, ssp1, hicetas)
conv_meno <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "menophilus")))), dist = 15000), dTolerance = 5000)
conv_meno.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Subsp_full == "Melinaea.menophilus.menophilus" |
                                                                                  Subsp_full == "Melinaea.menophilus.hicetas" |
                                                                                  Subsp_full == "Melinaea.menophilus.n. ssp. [1]")))), dist = 15000), dTolerance = 5000)

meno <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_meno.M)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "menophilus"))+
  tm_dots(size = 0.06, col = "red", shape = 21)

lf_meno <- tmap_leaflet(meno)
mapshot(lf_meno, file="~/phylodist-meno.jpg")

#  tm_shape(conv_meno)+
#  tm_borders(col = "red", lwd = 2)+
#  tm_fill(col = "red", alpha = 0.6)

## mothone = east
conv_mot <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "mothone")))), dist = 15000), dTolerance = 5000)
conv_mot.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Species == "mothone")))), dist = 15000), dTolerance = 5000)

mot <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_mot.M)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "mothone"))+
  tm_dots(size = 0.06, col = "red", shape = 21)

lf_mot <- tmap_leaflet(mot)
mapshot(lf_mot, file="~/phylodist-mot.jpg")

#  tm_shape(conv_mot)+
#  tm_borders(col = "red", lwd = 2)+
#  tm_fill(col = "red", alpha = 0.6)

## tarapotensis = east
conv_tar <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "tarapotensis")))), dist = 15000), dTolerance = 5000)
conv_tar.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Subsp_full == "Melinaea.satevis.tarapotensis")))), dist = 15000), dTolerance = 5000)

tar <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_tar.M)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "tarapotensis"))+
  tm_dots(size = 0.06, col = "red", shape = 21)

lf_tar <- tmap_leaflet(tar)
mapshot(lf_tar, file="~/phylodist-tar.jpg")

#  tm_shape(conv_tar)+
#  tm_borders(col = "red", lwd = 2)+
#  tm_fill(col = "red", alpha = 0.6)

## zaneka = east
conv_zan <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "zaneka")))), dist = 15000), dTolerance = 5000)
conv_zan.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Subsp_full == "Melinaea.menophilus.zaneka" |
                                                                                 Subsp_full == "Melinaea.menophilus.ernestoi")))), dist = 15000), dTolerance = 5000)

zan <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_zan.M)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "zaneka"))+
  tm_dots(size = 0.06, col = "red", shape = 21)
lf_zan <- tmap_leaflet(zan)
mapshot(lf_zan, file="~/phylodist-zan.jpg")

## west = orange
## lilis = west
conv_lil <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "lilis")))), dist = 15000), dTolerance = 5000)
conv_lil.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Subsp_full == "Melinaea.lilis.lilis" | 
                                                                                 Subsp_full == "Melinaea.lilis.parallelis")))), dist = 15000), dTolerance = 5000)

lil <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_lil.M)+
  tm_fill(col = "yellow",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "lilis"))+
  tm_dots(size = 0.06, col = "darkorange", shape = 21)

lf_lil <- tmap_leaflet(lil)
mapshot(lf_lil, file="~/phylodist-lil.jpg")

## idae = west
conv_idae <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "idae")))), dist = 15000), dTolerance = 5000)
conv_idae.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Subsp_full == "Melinaea.idae.idae")))), dist = 15000), dTolerance = 5000)

idae <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_idae.M)+
  tm_fill(col = "yellow",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "idae"))+
  tm_dots(size = 0.06, col = "darkorange", shape = 21)

lf_idae <- tmap_leaflet(idae)
mapshot(lf_idae, file="~/phylodist-idae.jpg")

#broad region 
## satevis = east to middle
conv_sat <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "satevis")))), dist = 15000), dTolerance = 5000)
conv_sat.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Subsp_full == "Melinaea.satevis.cydon" | 
                                                                                 Subsp_full == "Melinaea.satevis.maelus" |
                                                                                 Subsp_full == "Melinaea.satevis.lamasi")))), dist = 15000), dTolerance = 5000)

sat <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_sat.M)+
  tm_fill(col = "red",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "satevis"))+
  tm_dots(size = 0.06, col = "red", shape = 21)
lf_sat <- tmap_leaflet(sat)
mapshot(lf_sat, file="~/phylodist-sat.jpg")


## mneme = Guyanan shield + east
conv_mne <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "mneme", Region == "east")))), dist = 15000), dTolerance = 5000)
conv_mne2 <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "mneme", Region == "guiananshield")))), dist = 15000), dTolerance = 5000)
conv_mne.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Species == "mneme")))), dist = 15000), dTolerance = 5000)

mne <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_mne.M)+
  tm_fill(col = "black",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "mneme", Region == "east"))+
  tm_dots(size = 0.06, col = "red", shape = 21)+
  tm_shape(filter(dat_sf_mech, Category == "mneme", Region == "guiananshield"))+
  tm_dots(size = 0.06, col = "pink2", shape = 21)
lf_mne <- tmap_leaflet(mne)
mapshot(lf_mne, file="~/phylodist-mne.jpg")

## ludovica = east, guyana, Atlantic
conv_lud <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "ludovica", Region == "atlantic")))), dist = 15000), dTolerance = 5000)
conv_lud2 <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "ludovica", Region == "east")))), dist = 15000), dTolerance = 5000)
conv_lud3 <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Category == "ludovica", Region == "guiananshield")))), dist = 15000), dTolerance = 5000)
conv_lud.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael_mel, Subsp_full == "Melinaea.ludovica.ludovica" | 
                                                                                 Subsp_full == "Melinaea.ludovica.paraiya" |
                                                                                 Subsp_full == "Melinaea.ludovica.parayia")))), dist = 15000), dTolerance = 5000)

ludo <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_lud.M)+
  tm_fill(col = "black",alpha=0.4)+
  tm_shape(filter(dat_sf_mech, Category == "ludovica", Region == "atlantic"))+
  tm_dots(size = 0.06, col = "blue", shape = 21)+
  tm_shape(filter(dat_sf_mech, Category == "ludovica", Region == "east"))+
  tm_dots(size = 0.06, col = "red", shape = 21)+
  tm_shape(filter(dat_sf_mech, Category == "ludovica", Region == "guiananshield"))+
  tm_dots(size = 0.06, col = "pink2", shape = 21)

lf_ludo <- tmap_leaflet(ludo)
mapshot(lf_ludo, file="~/phylodist-ludo.jpg")


#### overview based on our own data ####

conv_east <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Region == "east")))), dist = 15000), dTolerance = 5000)
conv_west <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Region == "west")))), dist = 15000), dTolerance = 5000)
conv_atlantic <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Region == "atlantic")))), dist = 15000), dTolerance = 5000)
conv_guia <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_sf_mech, Region == "guiananshield")))), dist = 15000), dTolerance = 5000)

tmap_mode("view")
overview <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_east)+
  tm_borders(col = "red", lwd = 2)+
  tm_fill(col = "red",alpha=0.8)+
  tm_shape(conv_west)+
  tm_borders(col = "orange", lwd = 2)+
  tm_fill(col = "orange",alpha=0.8)+
  tm_shape(conv_guia)+
  tm_borders(col = "pink2", lwd = 2)+
  tm_fill(col = "pink2")+
  tm_shape(conv_atlantic)+
  tm_borders(col = "blue", lwd = 2)+
  tm_fill(col = "blue",alpha=0.8)

lf_overview <- tmap_leaflet(overview)
mapshot(lf_overview, file="~/phylodist-overview.jpg")

#east
conv_east.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael, Subsp_full == "Mechanitis.mazaeus.fallax" | Subsp_full == "Mechanitis.messenoides.messenoides" | 
                                                                                  Subsp_full == "Mechanitis.messenoides.deceptus" | Subsp_full == "Mechanitis.menapis.roqueensis" | Subsp_full == "Mechanitis.mazaeus.dorissides" | Subsp_full == "Mechanitis.mazaeus.proceriformis" | 
                                                                                  Subsp_full == "Mechanitis.mazaeus.proceriformes" | Species == "isocomma" | Subsp_full == "Melinaea.marsaeus.clara" | Subsp_full == "Melinaea.marsaeus.macaria" | Subsp_full == "Melinaea.marsaeus.rileyi" |
                                                                                  Subsp_full == "Melinaea.marsaeus.phasiana" | Subsp_full == "Melinaea.menophilus.menophilus" | Subsp_full == "Melinaea.menophilus.hicetas" | Subsp_full == "Melinaea.menophilus.n. ssp. [1]" | 
                                                                                  Subsp_full == "Melinaea.satevis.maeonis" | Species == "mothone" | Subsp_full == "Melinaea.satevis.tarapotensis" | Subsp_full == "Melinaea.menophilus.zaneka" | Subsp_full == "Melinaea.menophilus.ernestoi" |
                                                                                  Subsp_full == "Melinaea.satevis.cydon" | Subsp_full == "Melinaea.satevis.maelus" | Subsp_full == "Melinaea.satevis.lamasi")))), dist = 15000), dTolerance = 5000)

#west
conv_west.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael, Subsp_full == "Mechanitis.macrinus." | Subsp_full == "Mechanitis.lysimnia.utemaia" | Subsp_full == "Mechanitis.menapis.caribensis" | 
                                                                                  Subsp_full == "Mechanitis.menapis.dariensis" | Subsp_full == "Mechanitis.menapis.mantineus" | Subsp_full == "Mechanitis.menapis.menapis" | Subsp_full == "Mechanitis.menapis.occasiva" | Subsp_full == "Mechanitis.menapis.saturata" |
                                                                                  Subsp_full == "Mechanitis.polymnia.caucaensis" | Subsp_full == "Mechanitis.polymnia.chimborazona" | Subsp_full == "Mechanitis.polymnia.isthmia" | Subsp_full == "Melinaea.lilis.lilis" | Subsp_full == "Melinaea.lilis.parallelis" |
                                                                                  Subsp_full == "Melinaea.idae.idae")))), dist = 15000), dTolerance = 5000)

#atlantic forest
conv_braz.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael, Subsp_full == "Mechanitis.lysimnia.lysimnia" | Subsp_full == "Mechanitis.lysimnia.nesaea" | Subsp_full == "Mechanitis.polymnia.casabranca")))), dist = 15000), dTolerance = 5000)

#Guianan shield
conv_gui.M <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(filter(dat_mael, country == "Guyana" | country == "French Guiana" | country == "Suriname" | Subsp_full == "Mechanitis.polymnia.polymnia")))), dist = 15000), dTolerance = 5000)

tmap_mode("view")
overview <- tm_basemap("Esri.WorldTerrain")+
  tm_shape(filter(dat_sf_mech, Category == "corner"))+
  tm_dots(size = 0.02, col = "black", shape = 21)+
  tm_shape(conv_east.M)+
  tm_fill(col = "red",alpha=0.8)+
  tm_shape(conv_west.M)+
  tm_fill(col = "orange",alpha=0.8)+
  tm_shape(conv_gui.M)+
  tm_fill(col = "pink2", alpha = 0.8)+
  tm_shape(conv_braz.M)+
  tm_fill(col = "blue",alpha=0.8)+
  tm_shape(filter(dat_mael, Subsp_full == "Mechanitis.mazaeus.fallax" | Subsp_full == "Mechanitis.messenoides.messenoides" | 
                    Subsp_full == "Mechanitis.messenoides.deceptus" | Subsp_full == "Mechanitis.menapis.roqueensis" | Subsp_full == "Mechanitis.mazaeus.dorissides" | Subsp_full == "Mechanitis.mazaeus.proceriformis" | 
                    Subsp_full == "Mechanitis.mazaeus.proceriformes" | Species == "isocomma" | Subsp_full == "Melinaea.marsaeus.clara" | Subsp_full == "Melinaea.marsaeus.macaria" | Subsp_full == "Melinaea.marsaeus.rileyi" |
                    Subsp_full == "Melinaea.marsaeus.phasiana" | Subsp_full == "Melinaea.menophilus.menophilus" | Subsp_full == "Melinaea.menophilus.hicetas" | Subsp_full == "Melinaea.menophilus.n. ssp. [1]" | 
                    Subsp_full == "Melinaea.satevis.maeonis" | Species == "mothone" | Subsp_full == "Melinaea.satevis.tarapotensis" | Subsp_full == "Melinaea.menophilus.zaneka" | Subsp_full == "Melinaea.menophilus.ernestoi" |
                    Subsp_full == "Melinaea.satevis.cydon" | Subsp_full == "Melinaea.satevis.maelus" | Subsp_full == "Melinaea.satevis.lamasi"))+
  tm_dots(size = 0.06, col = "red", shape = 21)+
  tm_shape(filter(dat_mael, Subsp_full == "Mechanitis.lysimnia.lysimnia" | Subsp_full == "Mechanitis.lysimnia.nesaea" | Subsp_full == "Mechanitis.polymnia.casabranca"))+
  tm_dots(size = 0.06, col = "blue", shape = 21)+
  tm_shape(filter(dat_mael, Subsp_full == "Mechanitis.macrinus." | Subsp_full == "Mechanitis.lysimnia.utemaia" | Subsp_full == "Mechanitis.menapis.caribensis" | 
                    Subsp_full == "Mechanitis.menapis.dariensis" | Subsp_full == "Mechanitis.menapis.mantineus" | Subsp_full == "Mechanitis.menapis.menapis" | Subsp_full == "Mechanitis.menapis.occasiva" | Subsp_full == "Mechanitis.menapis.saturata" |
                    Subsp_full == "Mechanitis.polymnia.caucaensis" | Subsp_full == "Mechanitis.polymnia.chimborazona" | Subsp_full == "Mechanitis.polymnia.isthmia" | Subsp_full == "Melinaea.lilis.lilis" | Subsp_full == "Melinaea.lilis.parallelis" |
                    Subsp_full == "Melinaea.idae.idae"))+
  tm_dots(size = 0.06, col = "orange", shape = 21)+
  tm_shape(filter(dat_mael, country == "Guyana" | country == "French Guiana" | country == "Suriname" | Subsp_full == "Mechanitis.polymnia.polymnia"))+
  tm_dots(size = 0.06, col = "pink2", shape = 21)

lf_overview <- tmap_leaflet(overview)
mapshot(lf_overview, file="~/phylodist-overview.pdf")
