
functional_dist<-function(r,prob=0.4){
c=r
ver_nbd=matrix(nrow=r,ncol=r)
ver_nbd=apply(ver_nbd, c(1,2), function(x) sample(c(0,1),1,prob = c(1-prob,prob))) 
ver_nbd[upper.tri(ver_nbd)]<-0
ver_nbd<-ver_nbd+t(ver_nbd)
diag(ver_nbd)<-1

## finished making vertex neighbourhood matrix

# Create a label matrix to keep track of the vertex indices
loc_label=matrix(nrow = r,ncol=r)

lab=1;

for (i in 1:r){
  for (j in i:c){
    loc_label[i,j]=lab;
    lab=lab+1;
  }
}

# symmetrizing the label matrix
loc_label=t(loc_label)
loc_label[upper.tri(loc_label)]=t(loc_label)[upper.tri(loc_label)] 


# Create the functional overlap matrix and metrics 
n_pairs=r*(r+1)/2
pair_nbd=matrix(nrow=n_pairs,ncol=n_pairs)

for (i1 in 1:r){
  for(j1 in 1:c){
    for(i2 in 1:r){
      for(j2 in 1:c){
        if (i1==i2 & j1==j2){
          next;
        } else {pair_nbd[loc_label[i1,j1],loc_label[i2,j2]]=(ver_nbd[i1,i2]+ver_nbd[i1,j2]+ver_nbd[i2,j1]+ver_nbd[j1,j2])/4;
        }
      }
    }
  }
}

func_overlp=pair_nbd #this is the functional overlap matrix
func_dist=1/func_overlp #functional distance is inverse of functional overlap
func_dist=round(func_dist,2)

diag(func_dist)=0 #setting diagonal =0, since distance of each pair from itself is 0

func_dist[!is.finite(func_dist)] <- 8
write.csv(func_dist,file="func_dist.csv",row.names = FALSE)
#test_mat=read.csv("func_dist.csv")
}