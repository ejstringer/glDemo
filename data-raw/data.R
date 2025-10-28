
library(dartRverse)

#load all the data
# https://r-pkgs.org/data.html#sec-data-data-raw

#created via:


# load --------------------------------------------------------------------


tympos1 <- gl.read.dart("./data-raw/Report_DTym25-12345_10_moreOrders_SNP_1.csv", 
                        ind.metafile = "./data-raw/GED2025_corr.csv")
saveRDS(tympos1, './data-raw/tympos1.rds')

tympos1 <- readRDS('./data-raw/tympos1.rds')

plot(tympos1@other$loc.metrics$CallRate)
hist(tympos1@other$loc.metrics$CallRate)

# downsample loci ---------------------------------------------------------

nLoc(tympos1)
systematic_sample <- seq(1, nLoc(tympos1), length.out = 3210)

systematic_sample <- 1:5432
tympos <- gl.keep.loc(tympos1, loc.list = tympos1@loc.names[systematic_sample])
tympos


# wild pops only ----------------------------------------------------------

index <- tympos@other$ind.metrics$Colony=="Wild" 
index <- ifelse(is.na(index), FALSE, index)
tw <- tympos[index,]
index2 <- tw@other$ind.metrics$group !="Monaro"
tw2 <- tw[index2,]
tw2_5 <- gl.drop.ind(tw2, ind.list =  "AA61626")
tw2_5 <- gl.drop.ind(tw2, ind.list =  "CK1 hatchling")
tw_final <- tw2_5
nInd(tw_final)

# generic pop names -------------------------------------------------------

pop(tw_final) <- factor(ifelse(tw_final@other$ind.metrics$pop == 'Jerrabomberra West',
                               'W_Can',as.character(tw_final@other$ind.metrics$group)))

table(tw_final@pop)
tw_final@other$ind.metrics[tw_final@pop == 'unknown',]


# random lat lon ----------------------------------------------------------

gl.map.interactive(possums.gl)
poplocation <- cbind(id = possums.gl@ind.names, 
      pop = possums.gl@pop, 
      possums.gl@other$latlon) %>% 
  filter(pop %in% LETTERS[c(1:3,5:6)]) %>% 
  mutate(pop = case_when(
    pop == 'F' ~ 'A',
    pop == 'C' ~ 'B',
    .default = pop
  )) %>% 
  mutate(pop = case_when(
    pop == 'A' ~ 'N_Can',
    pop == 'E' ~ 'S_Can',
    pop == 'B' ~ 'W_Can'
  )) %>% 
  group_by(pop) %>% 
  summarise(min_lat = min(lat),
            max_lat = max(lat),
            min_lon = min(lon),
            max_lon = max(lon))

idpop <- data.frame(id = tw_final@ind.names, pop = tw_final@pop)

individuals<-left_join(idpop, poplocation)
table(tw_final@ind.names == individuals$id)
individuals<-individuals %>% 
  #filter(complete.cases(min_lat)) %>% 
  rowwise() %>% 
  mutate(#lat =runif(1, min_lat, max_lat),
         #lon = runif(1, min_lon, max_lon),
         lat = rnorm(1, mean = mean(c(min_lat, max_lat)), 0.01),
         lon = rnorm(1, mean = mean(c(min_lon, max_lon)), 0.01),
         lat = ifelse(is.nan(lat), NA, lat),
         lon = ifelse(is.nan(lon), NA, lon))

individuals %>% names

tw_final@other$latlon$lat <- individuals$lat
tw_final@other$latlon$lon <- individuals$lon

gl.map.interactive(tw_final)


# tidy metadata -----------------------------------------------------------

tw_final@other$ind.metrics %>% head

tw_final@other$ind.metrics$pop <- tw_final@pop
tw_final@other$ind.metrics$lat <- tw_final@other$latlon$lat
tw_final@other$ind.metrics$lon <- tw_final@other$latlon$lon
tw_final@other$ind.metrics$age <- tw_final@other$ind.metrics$Age

tw_final@other$ind.metrics[,c('id', 'pop', 'lat', 'lon', 'year', 'sex', 'age')] %>% head
tw_final@other$ind.metrics <- tw_final@other$ind.metrics[,c('id', 'pop', 'lat', 'lon', 'year', 'sex', 'age')]

tw_final@other$ind.metrics

# checks ------------------------------------------------------------------
#gl.subsample.ind(tw_final, n = 20)
tw_final

## smearplot ---------------------------------------------------------------
gl.smearplot(tw_final)

## filter ------------------------------------------------------------------

tw_final_filtered <- gl.filter.callrate(tw_final, 
                                        threshold = 0.95,method = "loc")
tw_final_filtered <- gl.filter.rdepth(tw_final_filtered,lower = 5, upper=50)
tw_final_filtered <- gl.filter.callrate(tw_final_filtered, threshold = 0.95,method = "ind")
tw_final_filtered <- gl.filter.reproducibility(tw_final_filtered)

nInd(tw_final_filtered)
nLoc(tw_final_filtered)

gl.report.monomorphs(tw_final_filtered)
## pcoa --------------------------------------------------------------------
pop(tw_final_filtered) <- tw_final_filtered@other$ind.metrics$pop
pc <- gl.pcoa(tw_final_filtered)
gl.pcoa.plot(pc, tw_final_filtered, yaxis = 1, xaxis = 2)



## diversity ---------------------------------------------------------------

pop(tw_final_filtered) <- paste(tw_final_filtered@other$ind.metrics$pop,
                       tw_final_filtered@other$ind.metrics$year,
                       sep = '-')
                       



ar <- gl.report.allelerich(tw_final_filtered)
ar$`Richness per population` %>% 
  tidyr::separate(pop, into = c('pop', 'year'), sep = '-') %>% 
#  filter(popsize>1) %>% 
  ggplot(aes(year, mean_richness, colour = pop))+
  geom_point(aes(size = popsize))+
  geom_smooth(method = 'lm', aes(group=pop))+
  theme_classic()

h <- gl.report.heterozygosity(tw_final_filtered)
h %>% 
  tidyr::separate(pop, into = c('pop', 'year'), sep = '-') %>% 
  ggplot(aes(year, He, colour = pop))+
  geom_point()+
  geom_smooth(method = 'lm', aes(group = pop))+
  theme_classic()



## structure ---------------------------------------------------------------

#gl.download.binary(software = 'structure', out.dir = getwd())
struct <- gl.run.structure(tw_final, k.range = 1:5,
                           exec = "./structure/structure.exe",noadmix = F)

gl.plot.structure(struct, K = 3)
gl.evanno(struct)



# csv keep files ----------------------------------------------------------

indkeep <- data.frame(id = tympos1@ind.names, 
                      keep = tympos1@ind.names %in% tw_final@ind.names)


lockeep <- data.frame(id = tympos1@loc.names, 
                      keep = tympos1@loc.names %in% tw_final@loc.names)

lockeep <- lockeep[rep(nrow(lockeep), each = 2),]
table(duplicated(lockeep$id))

write.csv(indkeep, './data-raw/indkeep.csv', row.names = F)
write.csv(lockeep, './data-raw/lockeep.csv', row.names = F)

