#"create weight adjacency matrix"

library(Rvcg)
library(rgl)
library(magrittr)
library(maptools)
library(geomorph)
library(car)
library(Matrix)
#setwd('/Users/wenbozhang/Desktop/high dimensional/DataAD')
path='\\\\fs2-vip/students/wenboz4/Desktop/high dimension project/DataAD/DataOut'
setwd(path)


######################################################################
# 1. Load the coordinates and triangular for first 100 subjects

p=32492
shape_files=grep('L.midthickness.32k_fs_LR.ply',list.files(),value=TRUE)[1:100]
fvalues=grep('L.thickness.32k_fs_LR.func.1D',list.files(),value=TRUE)[1:100]


#get coordinates
raw_coord=array(sapply(shape_files,function(x) t(vcgPlyRead(x, updateNormals = TRUE, clean = TRUE)$vb[1:3,]))
                ,dim=c(p,3,100))

#get coordinates after Procrustes analysis
#https://sabifo4.github.io/blog/Morphometrics_and_Procrustes_alignment 
coords=geomorph::gpagen(shapes)$coords

saveRDS(coords,'\\\\fs2-vip/students/wenboz4/Desktop/high dimension project/DataAD/coords.rds')



######################################################################
# 2.visualization to compare brain before and after Procrustes analysis

## two example brian
coords=readRDS('\\\\fs2-vip/students/wenboz4/Desktop/high dimension project/DataAD/coords.rds')

group=factor(c(rep(1,32492),rep(2,32492)))

raw=rbind(raw_coord[,,1],raw_coord[,,2])
scatter3d(x = raw[,1], y = raw[,2], z = raw[,3],
          groups=group, surface=FALSE)

coord=rbind(coords[,,1],coords[,,2])
scatter3d(x = coord[,1], y = coord[,2], z = coord[,3],
          groups=group, surface=FALSE)



######################################################################
# 3 get mean graph, distance and weighted adjacency matrix
#get triangle
coords=readRDS('\\\\fs2-vip/students/wenboz4/Desktop/high dimension project/DataAD/coords.rds')
trian= t(vcgPlyRead('A0002.L.midthickness.32k_fs_LR.ply', updateNormals = TRUE, clean = TRUE)$it)

get_MeanGraph<-function(x){
  #x is 32492*3*100
  z=matrix(0,nrow=32492,ncol=3)
  for(i in 1:100){
    z=z+x[,,i]
  }
  z/100
}

mean_graph=get_MeanGraph(coords)


#sort the node
sort_trian=t(apply(trian,1,sort))
dupe_edge=matrix(0,nrow=nrow(sort_trian)*3,ncol=2)

k=1
#get duplicated egde list
for (i in 1:nrow(sort_trian)){
  tri=sort_trian[i,]
  dupe_edge[k,1]=tri[1]
  dupe_edge[k,2]=tri[2]
  k=k+1
  dupe_edge[k,1]=tri[2]
  dupe_edge[k,2]=tri[3]
  k=k+1
  dupe_edge[k,1]=tri[1]
  dupe_edge[k,2]=tri[3]
  k=k+1
}

#half edge are duplicated
edge=unique(dupe_edge)

#calculate distance matrix
distance=rep(0,nrow(edge))
for(i in 1:nrow(edge)){
  coord1=edge[i,1]
  coord2=edge[i,2]
  distance[i]=sqrt(sum((mean_graph[coord1,]-mean_graph[coord2,])^2))
}

#creat weighted adjacency matrix
wei_adj=sparseMatrix(i=edge[,1],j=edge[,2],x=distance)