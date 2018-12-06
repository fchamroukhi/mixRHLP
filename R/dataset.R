source("R/utils.R")


MyData <- setRefClass(
  # Set the name for the class
  "MyData",

  # Define the fields
  fields = list(
    X="matrix",  # the data set
    XR="matrix", # reshaped data
    n="numeric",
    m="numeric",
    t="matrix"),

  # Set the methods
  methods=list(
    # set Data from a file
    setData = function(fileName){
      X <<- as.matrix(read.table(fileName,
                                 header = FALSE))

      setDataProperties(X)
    },
    # Generate example data
    generateExampleData = function(){
      nn<-c(10,10,10)
      n1<-nn[1]
      n2<-nn[2]
      n3<-nn[3]
      mean_flou <- as.matrix(read.table("data/mean_1_flou.txt",
                                 header = FALSE))
      y1 <- ones(n1,1)%*%t(mean_flou) + t(drnorm(length(mean_flou),n1,5,1))+1
      a <- drnorm(80,n2,7,1)
      b <- drnorm(130,n2,5,1)
      c <- drnorm(140,n2,4,1)
      y3 = t(rbind(a, b, c))

      a <- drnorm(120,n3,5,1)
      b <- drnorm(70,n3,7,1)
      c <- drnorm(160,n3,5,1)
      y2 <- t(rbind(a, b, c))
      X <<- rbind(y1, y2, y3)

      setDataProperties(X)
    },

    setDataProperties = function(X){
      XR <<- matrix(t(X), ncol = 1)
      n <<- nrow(X)
      m <<- ncol(X)
      #construction des matrices de regression
      t <<- t(matrix(seq(0, 1, length.out = m)));# ou rentrer le vecteur de covariables des courbes
    }
  )
)
