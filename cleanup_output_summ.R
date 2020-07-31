SE_output <- read.csv("output_SE.txt", sep="", header=FALSE)
SE_output_dim <- SE_output[c(2,3,5,6,8,9,11,12,14,15,17,18),]
SE_output_dim$V1 <- NULL
SE_output_dim$V3 <- NULL
rownames(SE_output_dim) <- c("All MATS events", "Significant MATS events", "All MISO Events", "significant MISO events", "All SUPPA events", "Significant SUPPA events", "MISO MATS overlap", "MISO MATS sig", "MISO SUPPA overlap", "MISO SUPPA sig", "MATS SUPPA overlap", "MATS SUPPA sig")
colnames(SE_output_dim) <- c("Events")
SE_output_Rsquare <- SE_output[c(77,141,205),]                    
SE_output_Rsquare$V1 <- NULL
SE_output_Rsquare$V2 <- NULL 
rownames(SE_output_Rsquare) <- c("Adj. Rsq. MATS MISO", "Adj. Rsq MISO SUPPA", "Adj. Rsq. MATS SUPPA")

RI_output <- read.csv("output_RI.txt", sep="", header=FALSE)
RI_output_dim <- RI_output[c(2,3,5,6,8,9,11,12,14,15,17,18),]
RI_output_dim$V1 <- NULL
RI_output_dim$V3 <- NULL
rownames(RI_output_dim) <- c("All MATS events", "Significant MATS events", "All MISO Events", "significant MISO events", "All SUPPA events", "Significant SUPPA events", "MISO MATS overlap", "MISO MATS sig", "MISO SUPPA overlap", "MISO SUPPA sig", "MATS SUPPA overlap", "MATS SUPPA sig")
colnames(RI_output_dim) <- c("Events")
RI_output_Rsquare <- RI_output[c(77,141,204),]                    
RI_output_Rsquare$V1 <- NULL
RI_output_Rsquare$V2 <- NULL 
rownames(RI_output_Rsquare) <- c("Adj. Rsq. MATS MISO", "Adj. Rsq MISO SUPPA", "Adj. Rsq. MATS SUPPA")

A5_output <- read.csv("output_A5.txt", sep ="", header=FALSE)
A5_output_dim <- A5_output[c(2,3,5,6,8,9,11,12,14,15,17,18),]
A5_output_dim$V1 <- NULL
A5_output_dim$V3 <- NULL
rownames(A5_output_dim) <- c("All MATS events", "Significant MATS events", "All MISO Events", "significant MISO events", "All SUPPA events", "Significant SUPPA events", "MISO MATS overlap", "MISO MATS sig", "MISO SUPPA overlap", "MISO SUPPA sig", "MATS SUPPA overlap", "MATS SUPPA sig")
colnames(A5_output_dim) <- c("Events")
A5_output_Rsquare <- A5_output[c(76,140,203),]                    
A5_output_Rsquare$V1 <- NULL
A5_output_Rsquare$V2 <- NULL 
rownames(A5_output_Rsquare) <- c("Adj. Rsq. MATS MISO", "Adj. Rsq MISO SUPPA", "Adj. Rsq. MATS SUPPA")

A3_output <- read.csv("output_A3.txt", sep="", header=FALSE)
A3_output_dim <- A3_output[c(2,3,5,6,8,9,11,12,14,15,17,18),]
A3_output_dim$V1 <- NULL
A3_output_dim$V3 <- NULL
rownames(A3_output_dim) <- c("All MATS events", "Significant MATS events", "All MISO Events", "significant MISO events", "All SUPPA events", "Significant SUPPA events", "MISO MATS overlap", "MISO MATS sig", "MISO SUPPA overlap", "MISO SUPPA sig", "MATS SUPPA overlap", "MATS SUPPA sig")
colnames(A3_output_dim) <- c("Events")
A3_output_Rsquare <- A3_output[c(74,137,200),]                    
A3_output_Rsquare$V1 <- NULL
A3_output_Rsquare$V2 <- NULL 
rownames(A3_output_Rsquare) <- c("Adj. Rsq. MATS MISO", "Adj. Rsq MISO SUPPA", "Adj. Rsq. MATS SUPPA")






sink(file= "clean_output")

print("SE")
SE_output_dim
SE_output_Rsquare

print("RI")
RI_output_dim
RI_output_Rsquare

print("A5")
A5_output_dim
A5_output_Rsquare

print("A3")
A3_output_dim
A3_output_Rsquare



sink()