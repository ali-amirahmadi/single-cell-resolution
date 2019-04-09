##statistical dat analyse




options(repos='http://cran.rstudio.com/')
setwd("D:/computer lessons/arshad/thesis/proposal/code/single-cell-resolution")

# distributed analysis belong x,y,z
geo = read.delim("data/geometry.txt", header = TRUE, sep = " ", dec = ".")

#x
x_ax = geo[1][]
head(x_ax)
typeof(x_ax)
x_ax = as.numeric(unlist(x_ax))
dx = density(x_ax)
plot(dx, main="Kernel Density of cells along x axis")
dev.copy(png,'result/visual/distbuted x.png')
dev.off()
#polygon(dx, col="red", border="blue") #for coloring chart
#y
y_ax = geo[2][]
head(y_ax)
y_ax = as.numeric(unlist(y_ax))
dy = density(y_ax)
plot(dy, main="Kernel Density of cells along y axis")
dev.copy(png,'result/visual/distbuted y.png')
dev.off()
#z
z_ax = geo[3][]
head(z_ax)
z_ax = as.numeric(unlist(z_ax))
dz = density(z_ax)
plot(dz, main="Kernel Density of cells along y axis")
dev.copy(png,'result/visual/distbuted z.png')
dev.off()
##### visulize  in 3d
#install.packages("rgl")
#install.packages('rgl', dependencies=TRUE, repos='http://cran.rstudio.com/')
library("rgl")
x = as.double(unlist(geo[1][]))
y = as.double(unlist(geo[2][]))
z = as.double(unlist(geo[3][]))

##
#' @param new.device a logical value. If TRUE, creates a new device
#' @param bg the background color of the device
#' @param width the width of the device
rgl_init <- function(new.device = FALSE, bg = "white", width = 640) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}
##


rgl.open() # Open a new RGL device
rgl.points(x, y, z, color ="lightgray") # Scatter plot
rgl.snapshot("result/visual/3d1.png")
library(magick)
movie3d(spin3d(axis = c(1, 1, 1)), duration = 3,
        dir = getwd())
##
rgl.open()# Open a new RGL device
rgl.bg(color = "white") # Setup the background color
rgl.points(x, y, z, color = "blue", size = 5) # Scatter plot
rgl.snapshot("result/visual/3d2.png")
library(magick)
movie3d(spin3d(axis = c(0, 0, 1)), duration = 3,
        dir = getwd())
##
rgl.open()# Open a new RGL device
rgl.bg(color = "white") # Setup the background color
rgl.spheres(x, y, z, r = 0.2, color = "grey") 
rgl.snapshot("result/visual/3d4.png")
library(magick)
movie3d(spin3d(axis = c(0, 0, 1)), duration = 3,
        dir = getwd())
##
rgl_init()
rgl.spheres(x, y, z, r = 0.2, color = "yellow")  # Scatter plot
rgl.bbox(color = "#333377") # Add bounding box decoration
rgl.snapshot("result/visual/3d6.png")
library(magick)
movie3d(spin3d(axis = c(0, 0, 1)), duration = 3,
        dir = getwd())
##
rgl_init()
rgl.spheres(x, y, z, r = 0.2, color = "yellow")  
# Add bounding box decoration
rgl.bbox(color=c("#333377","black"), emission="#333377",
         specular="#3333FF", shininess=5, alpha=0.8 ) 
##
# Make a scatter plot
rgl.open()
rgl.spheres(x, y, z, r = 0.2, color = "yellow") 
# Add x, y, and z Axes
rgl.lines(c(min(x), max(x)), c(0, 0), c(0, 0), color = "black")
rgl.lines(c(0, 0), c(min(y),max(y)), c(0, 0), color = "red")
rgl.lines(c(0, 0), c(0, 0), c(min(z),max(z)), color = "green")
##
# Make a scatter plot
rgl.open()
rgl.spheres(x, y, z, r = 0.2, color = "yellow") 
# Add x, y, and z Axes
rgl.lines(c(min(x), max(x)), c(0, 0), c(0, 0), color = "black")
rgl.lines(c(0, 0), c(min(y),max(y)), c(0, 0), color = "red")
rgl.lines(c(0, 0), c(0, 0), c(min(z),max(z)), color = "green")

##
library(magick)
rgl.open()
rgl.spheres(x, y, z, r = 0.2, color = "#D95F02") 
rgl_add_axes(x, y, z, show.bbox = TRUE)
# Compute and draw the ellipse of concentration
ellips <- ellipse3d(cov(cbind(x,y,z)), 
                    centre=c(mean(x), mean(y), mean(z)), level = 0.95)
wire3d(ellips, col = "#D95F02",  lit = FALSE)
aspect3d(1,1,1)
# Create a movie
movie3d(spin3d(axis = c(0, 0, 1)), duration = 3,
        dir = getwd())
##
library(rgl)
library("RColorBrewer")
rgl.open()
rgl.spheres(x, y, z, r = 0.2, color = "#D95F02") 
rgl_add_axes(x, y, z, show.bbox = T)
aspect3d(1,1,1)
rg.snapshot("result/plot.png")

####################################################
####################################
#######################
# GEN DAta visualization
#gen1 is normalize gen data from 1297 cell and contains 8924 gen
gen1 = read.delim("data/dge_normalized.txt", header = TRUE, sep = "\t", dec = ".")
gen1 = read.delim("data/bdtnp.txt", header = TRUE, sep = "\t", dec = ".")
geo = read.delim("data/geometry.txt", header = TRUE, sep = " ", dec = ".")
dim(geo)
typeof(gen1)
head(gen1)
dim(gen1)

# calculate var and sd for genes
gen1_var1 = apply(gen1, 2, var)
gen1_var2 = as.data.frame(gen1_var1)
dim(gen1_var2)
colnames(gen1_var2)



gen1_var3 = as.numeric(gen1_var1)
hist(gen1_var3, breaks=100, col="#BBFFDD", main="Variance per gene", 
     xlab="Variance", ylab="Number of genes")
dev.copy(png,'result/visual/variance among genes 8000.png')
dev.off()

gen1_sd1 = apply(gen1, 2, sd)
gen1_sd2 = as.data.frame(gen1_sd1)
dim(gen1_sd2)
gen1_sd3 = as.numeric(gen1_sd1)
hist(gen1_sd3, breaks=100, col="#BBFFDD", main="sd per gene", 
     xlab="sd", ylab="Number of genes")
dev.copy(png,'result/visual/sd among genes 8000.png')
dev.off()


###
###
x = c(10,10,10)
x = as.data.frame(x)
dim(x)
colnames(x) = c("gen1_var1")
colnames(x)
gen1_var4 = rbind(gen1_var2,x)
dim(gen1_var4)


#concatenate with geo
library(dplyr)
dim(geo)
gen_geo = cbind(gen1,geo)
dim(gen_geo)
colnames(gen_geo)

###
#concatenate variances and gen_cell
gen_cell_var = cbind(as.data.frame(t(gen_geo)),gen1_var4) 
dim(gen_cell_var)

#select 20 genes with most variance
j = gen_cell_var
gen_cell_var = j
#sort gens by variance
gen1_20var = gen_cell_var[order(- gen_cell_var$gen1_var1),]
dim(gen1_20var)
typeof(gen1_20var)

gen1_20var = as.data.frame(gen1_20var)
rownames(gen1_20var)
gen1_20var = gen1_20var[1:23,]
dim(gen1_20var)
#delete variance column
gen1_20var$gen1_var1 = NULL
dim(gen1_20var)

gen1_20var = as.data.frame(t(gen1_20var))
dim(gen1_20var)
colnames(gen1_20var)

################
#####devide data set to train and test
smp_siz = floor(0.75*nrow(gen1_20var))  # creates a value for dividing the data into train and test
set.seed(123)   # set seed to ensure you always have same random numbers generated
train_ind = sample(seq_len(nrow(gen1_20var)),size = smp_siz)  # Randomly identifies therows equal to sample size
train =gen1_20var[train_ind,] #creates the training dataset with row numbers stored in train_ind
test=gen1_20var[-train_ind,]  # creates the test dataset excluding the row numbers mentioned in train_ind


###############################
#####prediction model
### decision tree
##regression tree
library(rpart)
#grow tree
fit = rpart(xcoord ~ ftz+ sna+ eve+ twi+ tsh + odd+ danr+ Ama, method = "anova", data = train)

##calculate RMSE
pred_base <- predict(object=fit,
                     newdata = test)
options(repos='http://cran.rstudio.com/')# for instaliing package

#install.packages("Metrics")
library(Metrics)
rmse_base <- rmse(actual=test$xcoord, #Actual values
                  predicted = pred_base )

rmse_base
#so rmse base says our rmse erroe
######
printcp(fit) # display the results 

write.table(printcp(fit), "result/regression tree/error regression tree by 10 high var.txt", sep="\t")

plotcp(fit) # visualize cross-validation results 
dev.copy(png,'result/regression tree/plot cpregression tree by 10 high var.png')
dev.off()


summary(fit) # detailed summary of splits
m = as.character(summary(fit))
write.table(m, "result/regression tree/summary(fit) regression tree by 10 high var.txt", sep="\n")

print(fit)
m = as.character(print(fit))
write.table(m, "result/regression tree/print(fit) regression tree by 10 high var.txt", sep="\n")

# create additional plots 
par(mfrow=c(1,2)) # two plots on one page 
rsq.rpart(fit) # visualize cross-validation results
dev.copy(png,'result/regression tree/rsquare cpregression tree by 10 high var.png')
dev.off()

# plot tree 
par(mfrow=c(1,1)) # two plots on one page 
plot(fit, uniform=TRUE, 
     main="Regression Gen")
text(fit, use.n=TRUE, all=TRUE, cex=.7)
dev.copy(png,'result/regression tree/regression tree by 10 high var.png')
dev.off()

# create attractive postcript plot of tree 
post(fit, file = "result/tree2.ps", 
     title = "Regression Tree for Mileage ")


###baging: http://uc-r.github.io/regression_trees
set.seed(123)
library(rpart)
library(ipred)       # bagging
library(caret)       # bagging
ipred::bagging
bagged_m1 <- bagging(
  formula = xcoord ~ ftz+ sna+ eve+ twi+ tsh + odd+ danr+ Ama,
  data    = gen1_20var,
  coob    = TRUE
)

bagged_m1


# assess 10-50 bagged trees
ntree <- 10:50

# create empty vector to store OOB RMSE values
rmse <- vector(mode = "numeric", length = length(ntree))

for (i in seq_along(ntree)) {
  # reproducibility
  set.seed(123)
  
  # perform bagged model
  model <- bagging(
    formula = xcoord ~ ftz+ sna+ eve+ twi+ tsh + odd+ danr+ Ama,
    data    = gen1_20var,
    coob    = TRUE,
    nbagg   = ntree[i]
  )
  # get OOB error
  rmse[i] <- model$err
}

plot(ntree, rmse, type = 'l', lwd = 2)
abline(v = 25, col = "red", lty = "dashed")
rmse

#d

























##########################################
###############################################feature selection
######first we provide data independently
gen1 = read.delim("../data/bdtnp.txt", header = TRUE, sep = "\t", dec = ".")
geo = read.delim("../data/geometry.txt", header = TRUE, sep = " ", dec = ".")
dim(geo)
typeof(gen1)
head(gen1)
dim(gen1)
feature_selectd = cbind(gen1,geo)
dim(feature_selectd)
head(feature_selectd)
#########feature selections
######by bagging with cart package
# Specify 10-fold cross validation
ctrl <- trainControl(method = "cv",  number = 10) 
#delete ycoord and zcoord 
data_x = feature_selectd
data_x$ycoord = NULL
data_x$zcoord = NULL
dim(data_x)
# CV bagged model
bagged_cv = train(
  xcoord ~ .,
  data = data_x,
  method = "treebag",
  trControl = ctrl,
  importance = TRUE
)

# assess results
bagged_cv

plot(varImp(bagged_cv), 20)  


##########feature selection by sequencatial backward forward selection
#####http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/154-stepwise-regression-essentials-in-r/
###
library(caret)
library(leaps)
library(MASS)
data_x = feature_selectd
data_x$ycoord = NULL
data_x$zcoord = NULL
dim(data_x)
# Fit the full model 
full.model <- lm(xcoord ~., data = data_x)
summary(full.model)
# Stepwise regression model
step.model <- stepAIC(full.model, direction = "both", 
                      trace = FALSE)
summary(step.model)
###
models <- regsubsets(xcoord ~., data = data_x, nvmax = 5,
                     method = "seqrep")
summary(models)

#######
# Set seed for reproducibility
set.seed(123)
# Set up repeated k-fold cross-validation
train.control <- trainControl(method = "cv", number = 10)
# Train the model
step.model <- train(xcoord ~., data = data_x,
                    method = "leapSeq", 
                    tuneGrid = data.frame(nvmax = 1:85),
                    trControl = train.control
)
step.model$results
step.model$bestTune

####Adaptive-Network-Based Fuzzy Inference System
library("frbs")
step.model <- train(xcoord ~., data = data_x,
                    method = "ANFIS", 
                    max.iter = 10 ,
                    num.lable = 
                    tuneGrid = data.frame(nvmax = 1:85),
                    trControl = train.control
)
step.model$results
step.model$bestTune
###gradian boosting tree
library("e1071")
library("plyr")
library("ipred")
library(xgboost)
param <-  data.frame(nrounds=c(1000), max_depth = c(8),eta =c(0.3),gamma=c(0),
                     colsample_bytree=c(0.8),min_child_weight=c(1),subsample=c(1))
step.model <- train(xcoord ~., data = data_x,
                    method = "xgbTree", 
                    tuneGrid = param,
                    trControl = train.control
                    
)
step.model$results
step.model$bestTune
#### method = 'xgbDART'
param <-  data.frame(nrounds=c(1000), max_depth = c(4),eta =c(0.3),gamma=c(0),
                     colsample_bytree=c(0.8),min_child_weight=c(1),subsample=c(1),rate_drop = c(0.2),
                     skip_drop = c(0.2))
step.model <- train(xcoord ~., data = out4,
                    method = "xgbDART", 
                    tuneGrid = param,
                    trControl = train.control
                    
)
step.model$results
step.model$bestTune
###

###
##my clustring
# Set seed for reproducibility
set.seed(123)
# Set up repeated k-fold cross-validation
train.control <- trainControl(method = "cv", number = 10)
# Train the model
step.model <- train(xcoord ~., data = out4,
                    method = "leapSeq", 
                    tuneGrid = data.frame(nvmax = 1:85),
                    trControl = train.control
)
step.model$results
step.model$bestTune

summary(step.model$finalModel)


##

########################################
###########################################clustring for improve feature selection (bacward forward)
cl_data = feature_selectd
xcoord = cl_data$xcoord
cl_data$xcoord = NULL
cl_data$ycoord =NULL
cl_data$zcoord = NULL
mydata = cl_data
dim(mydata)
# Determine number of clusters
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:84) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:84, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")



#plot silhouette
require(cluster)
require(factoextra)
fviz_nbclust(mydata, kmeans, method = "silhouette")
##
# K-Means Cluster Analysis
fit <- kmeans(mydata, 4) # 2 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)
mydata$xcoord = xcoord
out = split(mydata, f= mydata$fit.cluster)
out1 = out[[1]]
out2 = out[[2]]
out3 = out[[3]]
out4 = out[[4]]
dim(out4)
max(out3$xcoord)
min(out3$xcoord)
###################### kmedoeids clustring
library(cluster)
library(factoextra)
fviz_nbclust(mydata, pam, method = "silhouette")
clus_pam = pam(mydata, 3)
print(clus_pam)
# append cluster assignment
mydata <- data.frame(mydata, clus_pam$cluster)
mydata$xcoord = xcoord
out = split(mydata, f= mydata$clus_pam.cluster)
out1 = out[[1]]
out2 = out[[2]]
out3 = out[[3]]

######## fuzzi clustring
library(cluster)
library(factoextra)
fanny(mydata, 4, metric = "euclidean", stand = FALSE)
fuzzy_clus = fanny(mydata, 4, metric = "euclidean", stand = FALSE)
print(fuzzy_clus)
# append cluster assignment
fuzzy_clus$membership
mydata <- data.frame(mydata, fuzzy_clus$membership)
mydata$xcoord = xcoord
dim(mydata)
head(mydata)
mydata$fit.cluster = NULL
##find best teshold for fuzzi clustring by ga
library(GA)
cell_fitness = function(b1,b2){
  mydata1 = mydata
  mydata1 = NULL
  mydata2 = mydata
  mydata2 = NULL
  
  for (cell in 1:nrow(mydata)) {
    if(mydata[cell, 87] > b1){
      mydata1 = rbind(mydata1,mydata[cell,])
    }
  }
  for (cell in 1:nrow(mydata)) {
    if(mydata[cell, 88] > b2){
      mydata2 = rbind(mydata1,mydata[cell,])
    }
  }
  mydata1$fuzzy_clus.cluster = NULL
  mydata1$X1 = NULL
  mydata1$X2 = NULL
  mydata1$x3 = NULL
  mydata1$x4 = NULL
  mydata2$fuzzy_clus.cluster = NULL
  mydata2$X1 = NULL
  mydata2$X2 = NULL
  mydata2$x3 = NULL
  mydata2$x4 = NULL
  
  
  train.control <- trainControl(method = "cv", number = 10)
  # Train the model
  step.model <- train(xcoord ~., data = mydata1,
                      method = "leapSeq", 
                      tuneGrid = data.frame(nvmax = 1:86),
                      trControl = train.control
  )
  res_df = step.model$results
  res_df = as.data.frame(res_df)
  best = as.data.frame(step.model$bestTune)
  
  train.control <- trainControl(method = "cv", number = 10)
  # Train the model
  step.model <- train(xcoord ~., data = mydata2,
                      method = "leapSeq", 
                      tuneGrid = data.frame(nvmax = 1:86),
                      trControl = train.control
  )
  res_df2 = step.model$results
  res_df2 = as.data.frame(res_df2)
  best2 = as.data.frame(step.model$bestTune)
  return (res_df[best[1,1],2] + res_df2[best2[1,1],2])
  
  
}



GA <- ga(type = "real-valued", 
         fitness =  function(x) -cell_fitness(x[1], x[2]),
         lower = c(.5, .5), upper = c(.8, .8), 
         popSize = 2, maxiter = 2
         )


plot(GA)
myresult = cell_fitness(.8)
mydata1 = mydata
mydata1 = NULL
##test for
#mydata1 = NULL
#mydata1 = as.data.frame(mydata1)
rm(mydata1)
count = 0
for (cell in 1:nrow(mydata)) {
  if(mydata[cell, 87] > .80){
    count = count + 1
    mydata1 = rbind(mydata1,mydata[cell,])
    #mydata1[cell,1:90] = mydata[cell,1:90]
    #print((mydata[cell,]))
  }
}
dim(mydata1)
count
mydata1
head(mydata1)
###
#gen2 is binarized gen data that points to 84 most important gene
gen2 = read.csv("data/dge_binarized_distMap.csv " ,check.names = FALSE)
typeof(gen2)
head(gen2)
dim(gen2)
typeof(gen1_var)
gen1_var2 = as.data.frame(gen1_var)
dim(gen1_var2)
var1 = var(gen1)
nrow(var1)
ncol(var1)

#//////data extraction review
gen1 = read.delim("data/dge_normalized.txt", header = TRUE, sep = "\t", dec = ".")
dim(gen1)
geo = read.delim("data/geometry.txt", header = TRUE, sep = " ", dec = ".")
dim(geo)
gen2 = read.csv("data/dge_binarized_distMap.csv " ,check.names = FALSE)
dim(gen2)
###########
####bdntp is our gen that we work on it
bdtnp1 = read.delim("data/bdtnp.txt", header = TRUE, sep = "\t", dec = ".")
dim(bdtnp1)
#####
dge_raw1 = read.delim("data/dge_raw.txt", header = TRUE, sep = "\t", dec = ".")
dim(dge_raw1)
dge_binarized_distMap1 = read.csv("data/dge_binarized_distMap.csv" ,check.names = FALSE)
dim(dge_binarized_distMap1)  






#///////////////
library(DistMap)
raw.data = read.csv("data/dge_raw.txt",sep = "\t",header = F)
rownames(raw.data) = raw.data$V1
raw.data$V1 = NULL

normalized.data = read.csv("data/dge_normalized.txt", sep = "\t")
insitu.matrix = read.csv("data/binarized_bdtnp.csv",check.names=F)


geometry = read.csv("data/geometry.txt",sep = " ")
dm = new("DistMap",
         raw.data=as.matrix(raw.data),
         data=as.matrix(normalized.data),
         insitu.matrix=as.matrix(insitu.matrix),
         geometry=as.matrix(geometry))
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dm <- mapCells(dm)
dm

computeVISH(dm, 'sna', threshold=0.75)
computeGeneGradient(dm, 'sna')



