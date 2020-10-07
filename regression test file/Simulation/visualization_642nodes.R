library(R.matlab)
library(Rvcg)
library(rgl)
library(shapes)

# setwd('/Users/yunbi/Library/Mobile Documents/com~apple~CloudDocs/20_6_summer/Eardi/sm-fpca-master')

#######################################################
# 1. Load coordinates and f values
vertices=readMat('vertices.mat')$vertices
faces=readMat('faces.mat')$faces
X=readMat('X.mat')$X
x_1=matrix(X[1,],ncol=1)

#######################################################
# 2. Visualize brainstem images

myColorRamp <- function(colors, values) {
        v <- (values - min(values))/diff(range(values))
        x <- colorRamp(colors)(v)
        rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

for (i in 1:nrow(X)){
        cols <- myColorRamp(c("white","blue"), matrix(X[i,],ncol=1))
        
        mesh_col <- tmesh3d(
                vertices = t(vertices),
                indices = t(faces),
                homogeneous = FALSE,
                material = list(color = cols)
        )
        
        shade3d(mesh_col, meshColor="vertices")
        Sys.sleep(.2)
}