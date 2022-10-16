setwd(".\\medium\\NMF\\lfw")
set.seed(123)
main.dir<-".\\medium\\NMF\\lfw"
my.list<-list(NULL)
index<-1
for(j in seq_along(dir())){ #reading all the image and storing into a list
  
  paste(main.dir,"\\",dir()[j], sep="")%>%setwd()
  
  for(i in seq_along(dir())){
    
    my.list[[index]]<-readImage(dir()[i])
    index<-index+1
    
  }
  setwd(main.dir)
  
}
#-------------------------Functions------------------------#
zscale<-function(x, u, sigma){ (x-mean(x, na.rm=T))*sigma/sd(x, na.rm = T) + u} #standardization function
minmax<-function(x){ (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)) } #min max function

NMF<-function(X, K, ite = 500, verbose){
  
  #..........Initialize.................
  N=nrow(X); M=ncol(X)
  W<-matrix(runif(N*K), nrow = N)
  H<-matrix(runif(K*M), nrow = K)
  
  loss<-NULL
  
  #................Update..............#
  for(i in 1:ite){
    
    H_num<-t(W)%*%X
    H_deno<-t(W)%*%W%*%H + 1e-10
    H<-H*(H_num/H_deno)
    
    W_num<-X%*%t(H)
    W_deno<-W%*%H%*%t(H) + 1e-10
    W<-W*(W_num/W_deno)
    
    loss[i]<-sum((X-W%*%H)^2)
    
    if(i%%500==0 & verbose == T){print(paste("iteration----------",i, sep = " "))
      print(sum((X-W%*%H)^2))}
    
  }
  
  return(list(Basis = W, Encoding = H, loss = loss))
}


#---------------------------------#
images=500 #number of images to be selected from the data base
pixel=150
sample.image<-sample(size = images, c(1:length(my.list)), replace = F)
my.new.list<-my.list[c(sample.image)]
#-------Step 1-------------------#
for(i in 1:images){my.new.list[[i]]<-c(channel(my.new.list[[i]], "grey")%>%resize(pixel, pixel))%>%zscale(u=0.25, sigma = 0.25)%>%minmax()}
#---------------Step 2--------------------------#
V<-unlist(my.new.list)%>%matrix(nrow = pixel*pixel, byrow = F)
#---------------Step 3------------------#
a<-NMF(X = V, K = 100, ite=1000, verbose = F)