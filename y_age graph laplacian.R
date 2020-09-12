library(Rvcg)
library(rgl)
library(shapes)
library(Matrix)
library(dplyr)
library(igraph)
library(PEIP)

setwd('/Users/yunbi/Library/Mobile Documents/com~apple~CloudDocs/20_6_summer/Eardi/brain-master/az_and_coords/brain/DataAD/DataOut')
setwd('DataOut_C')

p=32492
shape_files = grep('^A.+L\\.midthickness\\.32k_fs_LR\\.ply$', list.files(), value=TRUE)
thickness_files = grep('^A.+L\\.thickness\\.32k_fs_LR\\.func\\.1D', list.files(), value=TRUE)
rid = gsub(".L.+", "", shape_files)
rid = as.numeric(gsub("A0*", "", rid))

# regression setting with age variable
# predict age given the cortical thickness (the older, the thinner - furthermore, which specific part of brain is going to be thin or stable)
# Grp C: control, NAD: no az but with symptom, AD: az with symptom, PD: parkinson disease
age = read.csv('All_grp_ATN.csv', header=TRUE)
age = age[grep("A0*", age$Subject),]
age = age %>% 
        mutate(Subject=as.numeric(gsub("A0*","",Subject))) %>%
        arrange(age, Subject)
age = age[age$Subject %in% rid, ]
which(table(age$Subject)==2) # There are 2 rows for subject 4403 and 5040
age = age[-64,]
age = age[-88,]

# y: age - remove mean (intercept)
y_age = matrix(age$Age-mean(age$Age), nrow=nrow(age))

#####################################################################
az = read.csv('DXSUM_PDXCONV_ADNIALL.csv', header=TRUE)
az_age = az[az$VISCODE=="bl", c("RID", "VISCODE", "DXCURREN")]
az_age = az_age[az_age$RID %in% age$Subject, ]

# Why different ???
table(az_age$DXCURREN)
table(age$Grp)
#####################################################################

df = data.frame(rid = rid, sfiles = shape_files, tfiles = thickness_files)
df = df[df$rid %in% age$Subject,]

# raw coordinates
raw_cd = array(sapply(df$sfiles,
                      function(x) t(vcgPlyRead(x, updateNormals = TRUE, clean = TRUE)$vb[1:3,])), 
               dim=c(p, 3, length(df$sfiles)))

# coordinates after Procrustes analysis
# proc = procGPA(raw_cd, scale=FALSE)
proc = readRDS('proc_age_C.rds')
proc_cd = proc$rotated
template = proc$mshape # 32492*3 (coordinates of nodes)

# triangle mesh (64980 triangles * 3 nodes)
trg = t(vcgPlyRead('A0002.L.midthickness.32k_fs_LR.ply', updateNormals = TRUE, clean = TRUE)$it)

# f values (thickness) (87 subjects * 32492 nodes)
fdata = t(sapply(df$tfiles, function(x) as.matrix(read.csv(x, header = F))))

# egde set ((64980*3/2) undirected edges)
edge = matrix(0, nrow=nrow(trg)*3, ncol=2) # (64980*3) directed edges
k = 1
for (i in 1:nrow(trg)){
        tmp = trg[i,]
        edge[k,1] = tmp[1]
        edge[k,2] = tmp[2]
        k = k+1
        edge[k,1] = tmp[2]
        edge[k,2] = tmp[3]
        k = k+1
        edge[k,1] = tmp[3]
        edge[k,2] = tmp[1]
        k = k+1
}
edge_unq = t(apply(edge,1,sort))
edge_unq = unique(edge_unq) # (64980*3/2) undirected edges

# distance weight
# weight defined by Gaussian kernel? 
distance=rep(0, nrow(edge))
for(i in 1:nrow(edge)){
        v1=edge[i,1]
        v2=edge[i,2]
        distance[i]=sqrt(sum((template[v1,]-template[v2,])^2))
        # scaling? (makes not that much difference but will try it out later)
}

# adjacency matrix A (32492 nodes * 32492 nodes)
adj = sparseMatrix(i=edge[,1], j=edge[,2]) # FALSE (0) if there is no edge, TRUE (1) if there is an edge e_ij
wgt_adj = sparseMatrix(i=edge[,1], j=edge[,2], x=distance)

# (undirected) graph laplacian matrix
# from (64980*3/2) undirected edges
graph_obj = graph_from_adjacency_matrix(as.matrix(adj), mode="undirected")
lpl = laplacian_matrix(graph_obj)

#####################################################################
# takes too long time to run - we use sprase solver
# beta_ls = solve(t(fdata) %*% fdata) %*% t(fdata) %*% y
# beta_ridge_1 = solve(t(fdata) %*% fdata + lambda[1]*wgt_adj) %*% t(fdata) %*% y
# beta_ridge_100 = solve(t(fdata) %*% fdata + lambda[2]*wgt_adj) %*% t(fdata) %*% y

# lambda = c(1, 10, 100)
# beta_lpl = sapply(lambda, FUN=function(x) solve(t(fdata) %*% fdata + x*lpl) %*% t(fdata) %*% y)
# yhat = fdata %*% beta_lpl
#####################################################################

# 2/3 training data, 1/3 test data
train = sample(1:nrow(df), size=2/3*nrow(df))
test = -train
f_train = fdata[train,]
y_train = matrix(y_age[train,], nrow=length(train))

# take too long time to run
# beta_lpl_1 = solve(t(f_train) %*% f_train + 1*lpl) %*% t(f_train) %*% y_train # lambda 1
# yhat_test = fdata[test,] %*% beta_lpl

# y (Y), fdata (X), template, trg (penalty)
# solve it faster by using overdetermined linear system
chol_L = chol(lpl)
y_tilda = t(f_train) %*% y_train
A = rbind(f_train, chol_L)
y_tilda = rbind(y_tilda, matrix(rep(0, times=length(y_train)), ncol=1))
# overdetermined linear system solver
# lsei(A=A, B=y_tilda, fulloutput=TRUE)
# qr.solve(A, y_tilda)
gn = glmnet(A, y_tilda)