library(Rvcg)
library(rgl)

setwd("DataAD/DataOut")
# Load coordinates for each subject
index = c("A0002", "A0003", "A0004", "A0005", "A0007", "A0008", "A0010", 
          "A0014", "A0016", "A0019", "A0020", "A0021", "A0022", "A0023", "A0024", "A0029", 
          "A0030", "A0031", "A0033", "A0038", "A0040", "A0042", "A0044", "A0045", "A0048", 
          "A0051", "A0053", "A0054", "A0056", "A0058", "A0059", "A0060", "A0061", "A0066", 
          "A0076", "A0077", "A0078", "A0081", "A0083", "A0084", "A0089", "A0090", "A0093", "A0096", 
          "A0097", "A0098", "A0110", "A0111", "A0125", "A0126", "A0129", "A0130", "A0139", 
          "A0156", "A0166", "A0168", "A0173", "A0176", "A0177", "A0183", "A0204", "A0213", "A0217", 
          "A0219", "A0228", "A0240", "A0241", "A0243", "A0257", "A0262", "A0282", "A0284", 
          "A0290", "A0291", "A0292", "A0303", "A0304", "A0307", "A0311", "A0312", "A0314", 
          "A0325", "A0326", "A0331", "A0332", "A0336", "A0340", "A0341", "A0352", "A0354", "A0359")
rid = c(2,3,4,5,7,8,10,14,16,19,20,21,22,23,24,29,30,31,33,38,40,42,44,45,48,51,53,54,56,58,59,60,61,66,
        76,77,78,81,83,84,89,90,93,96,97,98,110,111,125,126,129,130,139,156,166,168,173,176,177,183,
        204,213,217,219,228,240,241,243,257,262,282,284,290,291,292,303,304,307,311,312,314,325,326,
        331,332,336,340,341,352,354,359)
seq = 1:91

mesh_p <- data.frame("rid" = rep(rid, each=32492),
                     "ix" = rep(1:32492, length(rid)),
                     "x" = rep(0, length(rid)*32492),
                     "y" = rep(0, length(rid)*32492),
                     "z" = rep(0, length(rid)*32492))
for(i in seq){
  tmp_shape = vcgPlyRead(paste0(index[i],".L.midthickness.32k_fs_LR.ply"), updateNormals = TRUE, clean = TRUE)
  mesh_p[{1+32492*(i-1)}:{32492*i},3:5] = data.frame(t(tmp_shape$vb[1:3,]))
}

# probably will not use it
#for(i in  seq){
#  tmp_shape = vcgPlyRead(paste0(index[i],".L.midthickness.32k_fs_LR.ply"), updateNormals = TRUE, clean = TRUE)
  #tmp = NULL
  #tmp$p = t(tmp_shape$vb[1:3,])
  
#  assign(paste0("mesh",rid[i]), data.frame(t(tmp_shape$vb[1:3,])))
#}

# registration (mean) and procrustes analysis

# same triangles for each subject
mesh_t = t(tmp_shape$it)

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
