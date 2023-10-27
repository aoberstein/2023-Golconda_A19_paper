library(MLSeq)
library(S4Vectors)
library(DESeq2)
library(caret)
library(pROC)
library(doParallel)

## adjust this to the number of desired cpu threads:
cl <- makePSOCKcluster(32)
registerDoParallel(cl)

# User defined variables --------------------------------------------------
ttgFile = "../ttg.tab"
metadataFile = "../metadata_merged.tab"
countFile = paste0("../2-limmaGene_GeTMM/", "GeTMM_allData_TMM_expValues_cpm.tab")
ttg = read.csv(ttgFile, header = T, sep = '\t', stringsAsFactors = F)
counts = read.csv(countFile, sep='\t', header= TRUE, stringsAsFactors = FALSE)
metadata = read.csv(metadataFile, sep='\t', header= TRUE, stringsAsFactors = FALSE)

## No rounding of TMM normalized counts log2 transform
counts2 = data.frame(row.names = counts[,1], counts[,2:length(colnames(counts))])
counts3 = log2(counts2+1)
# 
# ## filter HCMV, GFP, puro, and neo genes
gidsToRemove = c("EGFP_A206K", "neomycin", "puromycin")
counts4 = subset(counts3, !(row.names(counts3) %in% gidsToRemove))
cmvIds = ttg$ens_gene[grepl("viral", ttg$target_id)]
counts5 = subset(counts4, !row.names(counts4) %in% cmvIds)

## create class conversion table
classValuesEM = ifelse(metadata$cell == "RWPE-1"|metadata$cell == "MCF10A", "E", "M")
classConversionTable = data.frame(metadata$run, metadata$cell, classValuesEM)
names(classConversionTable)[1:3] = c("run", "cell", "class") 

###### Train classfier(s)
# Remove ARPE-19 samples, whose "E" vs "M" phenotype will be predicted
# using the resulting classifier.
A19_ids = subset(metadata$run, metadata$cell == "ARPE-19")

## select top x features having the highest
## gene-wise variances in order to decrease computational cost.
set.seed(2128)
vars <- sort(apply(counts5, 1, var, na.rm = TRUE), decreasing = TRUE)
counts6 <- counts5[names(vars)[1:5000], ]

# subset training and testing data
data.train <- as.matrix(counts6[, !colnames(counts6) %in% A19_ids])
data.test <- as.matrix(counts6[, colnames(counts6) %in% A19_ids])
samplesTrain = colnames(data.train) 
samplesTest = colnames(data.test)
class.train = DataFrame(condition = classConversionTable$class[match(samplesTrain, classConversionTable$run)])
class.test = DataFrame(condition = classConversionTable$class[match(samplesTest, classConversionTable$run)])
class.train = factor(classConversionTable$class[match(samplesTrain, classConversionTable$run)])
class.test = factor(classConversionTable$class[match(samplesTest, classConversionTable$run)])

## for QC (make sure classes are assigned to test and training data)
train.df = data.frame(samplesTrain, 
                      class.train, 
                      classConversionTable$run[match(samplesTrain, classConversionTable$run)],
                      classConversionTable$cell[match(samplesTrain, classConversionTable$run)],
                      classConversionTable$class[match(samplesTrain, classConversionTable$run)]
                      ) 
names(train.df)[3:5] = c("orig.run", "orig.cell", "orig.clas")
test.df = data.frame(samplesTest, 
                      class.test, 
                      classConversionTable$run[match(samplesTest, classConversionTable$run)],
                      classConversionTable$cell[match(samplesTest, classConversionTable$run)],
                      classConversionTable$class[match(samplesTest, classConversionTable$run)]
) 
names(test.df)[3:5] = c("orig.run", "orig.cell", "orig.class")

# final DESeqDataSet for MLSeq (not needed for caret)
# data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = class.train,
#                                       design = formula(~condition))
# data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = class.test,
#                                      design = formula(~1))
#####################################################################################
## Caret svmLinear
## see: https://cafernandezlo.github.io/es_fic_mubics_caret/caret.html#Machine_Learning_models:_Random_Forest_and_Linear_Discriminant_Analysis
## ranger is a fast version of Random Forest
cv.k  <- 10
reps <- 5

train.control <- trainControl(method = "repeatedcv", number = cv.k,
                              repeats = reps, verboseIter = TRUE,
                              classProbs = TRUE, savePredictions = TRUE,
                              allowParallel = TRUE)
# linear SVM model:
caretSVM.model <- train(x = t(data.train),
                  y = class.train,
                  method = "svmLinear",
                  metric = "Accuracy",
                  trControl = train.control)

show(caretSVM.model)
confusionMatrix(caretSVM.model)

# for svmLinear
selectedIndices <- caretSVM.model$pred$C == 1 

plot.roc(caretSVM.model$pred$obs[selectedIndices],
         caretSVM.model$pred$M[selectedIndices])
pred.caretSVM = predict(caretSVM.model, t(data.test))
pred.caretSVM
data.frame(test.df, pred.caretSVM)
correct.caretSVM = length(pred.caretSVM[pred.caretSVM == test.df$orig.class])/length(pred.caretSVM)
correct.caretSVM

####################################################################
## Caret svmPoly
cv.k  <- 10
reps <- 5
# hyperparameters <- expand.grid(mtry = c(3, 5, 7),
#                                min.node.size = c(2, 4, 6),
#                                splitrule = "gini")

train.control <- trainControl(method = "repeatedcv", number = cv.k,
                              repeats = reps, verboseIter = TRUE,
                              classProbs = TRUE, savePredictions = TRUE,
                              allowParallel = TRUE)
# polynomial basis kernel:
caretSVMpoly.model <- train(x = t(data.train),
                        y = class.train,
                        method = "svmPoly",
                        # tuneGrid = hyperparameters,
                        metric = "Accuracy",
                        # importance= 'impurity',
                        trControl = train.control)

show(caretSVMpoly.model)
confusionMatrix(caretSVMpoly.model)
# for SVMpoly
selectedIndices <- caretSVMpoly.model$pred$degree == 1 & caretSVMpoly.model$pred$scale == 0.001 & caretSVMpoly.model$pred$C == 0.25
plot.roc(caretSVMpoly.model$pred$obs[selectedIndices],
         caretSVMpoly.model$pred$M[selectedIndices])
pred.caretSVMpoly = predict(caretSVMpoly.model, t(data.test))
pred.caretSVMpoly
data.frame(test.df, pred.caretSVMpoly)
correct.caretSVMpoly = length(pred.caretSVMpoly[pred.caretSVMpoly == test.df$orig.class])/length(pred.caretSVMpoly)
correct.caretSVMpoly
# print(caretSVMpoly.model)
# predictors(caretSVMpoly.model)
# selectedGenes(caretSVMpoly.model)
# #####################################################################################
# ## Caret Random forest
# ## see: https://cafernandezlo.github.io/es_fic_mubics_caret/caret.html#Machine_Learning_models:_Random_Forest_and_Linear_Discriminant_Analysis
# ## ranger is a fast version of Random Forest
cv.k  <- 10
reps <- 5
hyperparameters <- expand.grid(mtry = c(3, 5, 7),
                               min.node.size = c(2, 4, 6),
                               splitrule = "gini")

train.control <- trainControl(method = "repeatedcv", number = cv.k,
                              repeats = reps,
                              returnResamp = "final", verboseIter = TRUE,
                              allowParallel = TRUE)

caretRF.model <- train(x = t(data.train),
                  y = class.train,
                  method = "ranger",
                  tuneGrid = hyperparameters,
                  metric = "Accuracy",
                  importance= 'impurity',
                  trControl = train.control,
                  num.trees = 500)

show(caretRF.model)
confusionMatrix(caretRF.model)

pred.caretRF = predict(caretRF.model, t(data.test))
pred.caretRF
data.frame(test.df, pred.caretRF)
correct.caretRF = length(pred.caretRF[pred.caretRF == test.df$orig.class])/length(pred.caretRF)
correct.caretRF


# #####################################################################################
# ## Caret Binomial (two classes) Logistic Regression
# ## see: https://cafernandezlo.github.io/es_fic_mubics_caret/caret.html#Machine_Learning_models:_Random_Forest_and_Linear_Discriminant_Analysis
# ## ranger is a fast version of Random Forest
cv.k  <- 10
reps <- 5
# hyperparameters <- expand.grid(mtry = c(3, 5, 7),
#                                min.node.size = c(2, 4, 6),
#                                splitrule = "gini")

train.control <- trainControl(method = "repeatedcv", number = cv.k,
                              repeats = reps,
                              returnResamp = "final", verboseIter = TRUE,
                              allowParallel = TRUE)

caretLogistic.model <- train(x = t(data.train),
                       y = class.train,
                       method = "glmnet",
                       family = "binomial",
                       # tuneGrid = hyperparameters,
                       metric = "Accuracy",
                       importance = 'impurity',
                       trControl = train.control)

show(caretLogistic.model)
confusionMatrix(caretLogistic.model)

pred.caretLogistic = predict(caretLogistic.model, t(data.test))
pred.caretLogistic
data.frame(test.df, pred.caretLogistic)
correct.caretLogistic = length(pred.caretLogistic[pred.caretLogistic == test.df$orig.class])/length(pred.caretLogistic)
correct.caretLogistic

# #####################################################################################
# ## KNN Regression
# ## see: https://cafernandezlo.github.io/es_fic_mubics_caret/caret.html#Machine_Learning_models:_Random_Forest_and_Linear_Discriminant_Analysis
# ## ranger is a fast version of Random Forest
cv.k  <- 10
reps <- 5
# hyperparameters <- expand.grid(mtry = c(3, 5, 7),
#                                min.node.size = c(2, 4, 6),
#                                splitrule = "gini")

train.control <- trainControl(method = "repeatedcv", number = cv.k,
                              repeats = reps, classProbs = TRUE,
                              returnResamp = "final", verboseIter = TRUE,
                              allowParallel = TRUE)

caretKNN.model <- train(x = t(data.train),
                             y = class.train,
                             method = "knn",
                             # tuneGrid = hyperparameters,
                             metric = "Accuracy",
                             trControl = train.control)

show(caretKNN.model)
confusionMatrix(caretKNN.model)

pred.caretKNN = predict(caretKNN.model, t(data.test))
pred.caretKNN
data.frame(test.df, pred.caretKNN)
correct.caretKNN = length(pred.caretKNN[pred.caretKNN == test.df$orig.class])/length(pred.caretKNN)
correct.caretKNN

# #####################################################################################
# ## Naive Bayes Classifier
# ## see: https://cafernandezlo.github.io/es_fic_mubics_caret/caret.html#Machine_Learning_models:_Random_Forest_and_Linear_Discriminant_Analysis
# ## ranger is a fast version of Random Forest
cv.k  <- 10
reps <- 5
# hyperparameters <- expand.grid(mtry = c(3, 5, 7),
#                                min.node.size = c(2, 4, 6),
#                                splitrule = "gini")

train.control <- trainControl(method = "repeatedcv", number = cv.k,
                              repeats = reps,
                              returnResamp = "final", verboseIter = TRUE,
                              allowParallel = TRUE)

caretNB.model <- train(x = t(data.train),
                        y = class.train,
                        method = "nb",
                        # usepoisson = TRUE,
                        # tuneGrid = hyperparameters,
                        metric = "Accuracy",
                        trControl = train.control)

show(caretNB.model)
confusionMatrix(caretNB.model)

pred.caretNB = predict(caretNB.model, t(data.test))
pred.caretNB
data.frame(test.df, pred.caretNB)
correct.caretNB = length(pred.caretNB[pred.caretNB == test.df$orig.class])/length(pred.caretNB)
correct.caretNB

################################# summarize data
summaryDat = data.frame(test.df, pred.caretKNN, pred.caretLogistic, pred.caretNB, pred.caretSVM, pred.caretSVMpoly, pred.caretRF)
write.table(file = "ARPE19_caretML_predictions.tab", summaryDat, row.names = F, sep = '\t', quote = F)

