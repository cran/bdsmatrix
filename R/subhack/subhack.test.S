tmat <- bdsmatrix(c(3,2,2,4), 
              c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
              matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dimnames(tmat) <- list(NULL, letters[1:13])

subhack(tmat, c(1, 14), c(2,1), 21:22)  #should fail
subhack(tmat, c(1:3), c(2,1), 21:22)  #should fail
subhack(tmat, 7,1,2)      #should fail, sparse

temp <- subhack(tmat, c(1, 11,11), c(3, 9,10), -(1:3))
temp[1,3] == -1
temp[11,9] == -2
temp[9,11] == -2
temp[11,10] ==-3
