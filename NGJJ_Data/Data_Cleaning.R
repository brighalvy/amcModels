## New data organization:
library(tidyverse)
rawdat <- read.csv("ngW3.csv")

# Fix class categories:
rawdat[rawdat$Year == "2022" & rawdat$Category == "99",]$Category <- "99plus"
rawdat[rawdat$Year == "2022" & rawdat$Category == "-99",]$Category <- "99"
rawdat <- rawdat[rawdat$Category != "Absolute",]

## Fix count data:
#############################
## Prep for Count Matrices ##
#############################
codes <- c("CLI1", "CLI2", "TKDA1", "BTKA1", "SUBA1", "BTK1",  "PGD1",  "SWP1",  "GPSA1", "GPS1", "SWPA1", "TKD1", "TKDA2", "BTKA2", "SUBA2", "BTK2",  "PGD2",  "SWP2",  "GPSA2", "GPS2",  "SWPA2",  "TKD2", "SUB1", "SUB2", "RST1", "RST2")
#col.codes <- c("CLI", "TKDA1", "BTKA1", "SUBA1", "BTK1",  "PGD1",  "SWP1",  "GPSA1", "GPS1", "SWPA1", "TKD1", "TKDA2", "BTKA2", "SUBA2", "BTK2",  "PGD2",  "SWP2",  "GPSA2", "GPS2",  "SWPA2",  "TKD2", "SUB1", "SUB2", "RST")
counts <- array(0, dim = c(5,26, 26), dimnames = list(unique(rawdat$Category), codes, codes))
for(i in 1:(nrow(rawdat) - 1)){
  if(rawdat$Combat[i] != rawdat$Combat[i + 1]){
    print("hello")
  } else if(rawdat$Athlete[i] == rawdat$Athlete[i+1]){
    cat <- rawdat$Category[i]
    row.code <- paste(rawdat$Code[i], "1", sep = "")
    col.code <- paste(rawdat$Code[i+1], "1", sep = "")
    counts[cat, row.code, col.code] <- counts[cat, row.code, col.code] + 1
  } else if(rawdat$Athlete[i] != rawdat$Athlete[i+1]){
    cat <- rawdat$Category[i]
    row.code <- paste(rawdat$Code[i], "1", sep = "")
    col.code <- paste(rawdat$Code[i+1], "2", sep = "")
    counts[cat, row.code, col.code] <- counts[cat, row.code, col.code] + 1
  } else{
    print(paste("trouble on row",i))
  }
}

## Get analysis ready counts matrices:

counts.fun <- function(dat){
  dat <- data.frame(dat)
  dat$RST <- dat$RST1 + dat$RST2
  dat$CLI <- dat$CLI1 + dat$CLI2
  dat[is.na(dat$RST), "RST"] = 0
  dat[is.na(dat$CLI), "CLI"] = 0 
  dat$RST = dat$RST + dat$CLI
  col.codes <- c("CLI", "TKDA1", "BTKA1", "SUBA1", "BTK1",  "PGD1",  "SWP1",  "GPSA1", "GPS1", "SWPA1", "TKD1", "TKDA2", "BTKA2", "SUBA2", "BTK2",  "PGD2",  "SWP2",  "GPSA2", "GPS2",  "SWPA2",  "TKD2", "SUB1", "SUB2", "RST")
  row.codes <- c("TKDA1", "BTKA1", "SUBA1", "BTK1",  "PGD1",  "SWP1",  "GPSA1", "GPS1", "SWPA1", "TKD1")
  dat = dat[row.codes, col.codes]
  #dat <- dat %>% relocate(SUB1, .after = SUBA2)

  return(dat)
}

# Fix this application the function works.
counts.fin <- apply(counts, 1, counts.fun)


## Code impossible moves
jj.data <- read.csv("countdata.csv", header = T, sep = ";")
rownames(jj.data) <- jj.data[,1]
jj.data = jj.data[,-1]
jj.data[is.na(jj.data$RST), "RST"] = 0
jj.data[is.na(jj.data$CLI), "CLI"] = 0 
jj.data$RST = jj.data$RST + jj.data$CLI
jj.data = jj.data[-1,-1]
head(jj.data)

jj.data <- jj.data %>% relocate(SUB1, .after = SUBA2)
q.counts1 <- jj.data[1:10,1:10] + jj.data[11:20,11:20]
q.counts2 <- jj.data[1:10, 11:20] + jj.data[11:20,1:10]
s.counts1 <- jj.data[1:10, 21] + jj.data[11:20,22]
s.counts2 <- jj.data[1:10,22] + jj.data[11:20, 21]
s.counts3 <- jj.data[1:10,23]+jj.data[11:20,23]
counts <- cbind(q.counts1, q.counts2, s.counts1, s.counts2, s.counts3)
names(counts) <- names(jj.data[1:10,])

## reorder counts.fin to match:
for(k in 1:5){
  counts.fin[[k]] <- counts.fin[[k]][rownames(counts), colnames(counts)]
}


## Manipulate to have impossible moves:
count.dat <- array(NA, dim = c(5, 10, 23), dimnames = list(unique(rawdat$Category), rownames(counts), colnames(counts)))
for(i in 1:5){
  counts.fin[[i]][is.na(counts)] <- NA
  count.dat[i,,] <- as.matrix(counts.fin[[i]])
}

## Save Data:
saveRDS(count.dat, "count.dat.3.RDS")
