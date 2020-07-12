library(Rvcg)
library(rgl)

# Load mesh coordinates for each subject
mesh_shape2 = vcgPlyRead('DataAD/DataOut/A0002.L.midthickness.32k_fs_LR.ply', updateNormals = TRUE, clean = TRUE)
mesh2 = NULL
mesh2$p = t(mesh_shape2$vb[1:3,])

mesh_shape3 = vcgPlyRead('DataAD/DataOut/A0003.L.midthickness.32k_fs_LR.ply', updateNormals = TRUE, clean = TRUE)
mesh3 = NULL
mesh3$p = t(mesh_shape3$vb[1:3,])

mesh_shape4 = vcgPlyRead('DataAD/DataOut/A0004.L.midthickness.32k_fs_LR.ply', updateNormals = TRUE, clean = TRUE)
mesh4 = NULL
mesh4$p = t(mesh_shape4$vb[1:3,])

mesh_shape5 = vcgPlyRead('DataAD/DataOut/A0005.L.midthickness.32k_fs_LR.ply', updateNormals = TRUE, clean = TRUE)
mesh5 = NULL
mesh5$p = t(mesh_shape5$vb[1:3,])

mesh_shape7 = vcgPlyRead('DataAD/DataOut/A0007.L.midthickness.32k_fs_LR.ply', updateNormals = TRUE, clean = TRUE)
mesh7 = NULL
mesh7$p = t(mesh_shape7$vb[1:3,])

mesh2$t = t(mesh_shape2$it)
identical(t(mesh_shape2$it), t(mesh_shape7$it))

# Load functional data
mesh2$f = as.matrix(read.csv('DataAD/DataOut/A0002.L.thickness.32k_fs_LR.func.1D', header = F))
dim(mesh2$p) # Mesh nodes
dim(mesh2$t) # Mesh triangles
dim(mesh2$f) # Data on nodes

myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

cols2 <- myColorRamp(c("white", "blue"), mesh2$f) 

mesh_col2 <- tmesh3d(
  vertices = t(mesh2$p),
  indices = t(mesh2$t),
  homogeneous = FALSE,
  material = list(color = cols2)
)

shade3d(mesh_col2, meshColor = "vertices")
