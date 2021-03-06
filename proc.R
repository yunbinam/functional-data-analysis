library(Rvcg)
library(rgl)
library(shapes)
library(Matrix)

# setwd('/Users/yunbi/Library/Mobile Documents/com~apple~CloudDocs/20_6_summer/Eardi/brain-master/DataAD')
# setwd('DataOut')
# setwd('/Users/yunbi/Library/Mobile Documents/com~apple~CloudDocs/20_6_summer/Eardi/brain-master/procrustes analysis and adjacency matrix')

######################################################################
# 1. Load coordinates for each subject
## use package 'shapes' for Procrustes analysis
## option for not normalizing(scaling)

p=32492
shape_files = grep('^A.+L\\.midthickness\\.32k_fs_LR\\.ply$', list.files(), value=TRUE)
thickness_files = grep('^A.+L\\.thickness\\.32k_fs_LR\\.func\\.1D', list.files(), value=TRUE)
rid = gsub(".L.+", "", shape_files)
rid = as.numeric(gsub("A0+", "", rid))

az = read.csv('DXSUM_PDXCONV_ADNIALL.csv', header=TRUE)
az_dgn = az[az$VISCODE=="bl", c("RID", "VISCODE", "DXCURREN")]
az_dgn = az_dgn[az_dgn$RID %in% rid, ]

# table(az_dgn$DXCURREN) # 1(healthy): 36, 2(mild): 33, 3(AZ): 19 (88 in total)
# rid[!(rid %in% az_dgn$RID)] # we do not have diagnosis data for RID # 20, 24, 340

df = data.frame(rid = rid, sfiles = shape_files, tfiles = thickness_files)
df = df[df$rid %in% az_dgn$RID,]
az1_dgn = df[az_dgn$DXCURREN==1,]
az2_dgn = df[az_dgn$DXCURREN==2,]
az3_dgn = df[az_dgn$DXCURREN==3,]

# raw coordinates
raw_cd = array(sapply(shape_files,
                      function(x) t(vcgPlyRead(x, updateNormals = TRUE, clean = TRUE)$vb[1:3,])), 
               dim=c(p, 3, length(shape_files)))

# coordinates after Procrustes analysis
# https://cran.r-project.org/web/packages/shapes/shapes.pdf

# proc = procGPA(raw_cd, scale=FALSE)
# proc_cd = proc$rotated
# template = proc$mshape
saveRDS(proc_cd, 'proc.RDS')
saveRDS(template, 'template.RDS')

# triangle mesh
trg = t(vcgPlyRead('A0002.L.midthickness.32k_fs_LR.ply', updateNormals = TRUE, clean = TRUE)$it)

# f values (thickness)
fdata = array(sapply(thickness_files,
                     function(x) as.matrix(read.csv(x, header = F))),
              dim=c(p, 3, length(thickness_files)))

save.image(file='proc.RData')
load('proc.RData')
