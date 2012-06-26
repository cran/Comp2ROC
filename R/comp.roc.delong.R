comp.roc.delong <-
function(sim1.ind,sim1.sta,sim2.ind,sim2.sta) {
imax=length(sim1.sta[sim1.sta==0])
jmax=length(sim2.sta[sim2.sta==1])
total_cases=imax+jmax
dados=rbind(cbind(sim1.ind,sim1.sta),cbind(sim2.ind,sim2.sta))

mod=2

# norm defines the negative cases
norm=array(data=NA,c(imax,mod))
# abnorm defines the positive cases
abnorm=array(data=NA,c(jmax,mod))
# data are atributed to the arrays norm and abnorm
for (i in 1:imax) {
    for (k in 1:mod){
        norm[i,k]=dados[i+(total_cases*(k-1)),1]
    }
}
for (j in 1:jmax) {
    for (k in 1:mod){
     abnorm[j,k]=dados[j+imax+(total_cases*(k-1)),1]
    }
}
# T1 defines the Wilcoxon Mann Whitney matrix for each modality
T1=array(data=NA,c(jmax,imax,mod))
for (k in 1:mod) {
    for (i in 1:imax){
	      for (j in 1:jmax){
	                  dif=abnorm[j,k]-norm[i,k]
	                  if(dif>0) T1[j,i,k]=1 else (if(dif==0) T1[j,i,k]=0.5 else T1[j,i,k]=0)
        }
    }
}
# U and V are the arrays obtained by summing the columns and the rouws of T1
U=array(data=NA,c(jmax,mod))
V=array(data=NA,c(imax,mod))
for (k in 1:mod) {
    for (j in 1:jmax){
    U[j,k]=sum(T1[j,,k])/imax
    }
    for (i in 1:imax){
    V[i,k]=sum(T1[,i,k])/jmax
    }
}
# Areas under the ROC curves for each modality
AUC=array(data=NA,c(mod))
for (k in 1:mod){
    AUC[k]=1/(imax*jmax)*sum(T1[,,k])
}
# Standard deviations for eachj modality
S=array(data=NA,c(mod))
for (k in 1:mod) {
     S[k]=(1/imax*sum((V[,k]-AUC[k])^2)/(imax-1)+1/jmax*sum((U[,k]-AUC[k])^2)/(jmax-1))^(0.5)
}
# Vari?nces for each modality, positive and negative cases and global
S2V=array(data=NA,c(mod,mod))
S2U=array(data=NA,c(mod,mod))
S2=array(data=NA,c(mod,mod))
for (k in 1:mod){
     for (z in 1:mod){
          S2V[k,z]= sum((V[,k]-AUC[k])*(V[,z]-AUC[z]))/(imax-1)
          S2U[k,z]= sum((U[,k]-AUC[k])*(U[,z]-AUC[z]))/(jmax-1)
          S2[k,z]=1/imax*S2V[k,z]+1/jmax*S2U[k,z]
     }
}
# Global correlations
R=array(data=NA,c(mod,mod))
for (k in 1:mod){
     for (z in 1:mod){
     if (z==k) {
      if (S2[k,k]*S2[z,z]==0)
        R[k,z]=1
      else
        R[k,z]=S2[k,z]/((S2[k,k]*S2[z,z])^(0.5))
        }
      else
      {
        if (S2[k,k]*S2[z,z]==0)
          R[k,z]=0
        else
          R[k,z]=S2[k,z]/((S2[k,k]*S2[z,z])^(0.5))
      }
     }
}
Z=(AUC[1]-AUC[2])/sqrt(S2[1,1]+S2[2,2]-2*R[1,2]*sqrt(S2[1,1]*S2[2,2]))
if (Z>0) p.value=2*(1-pnorm(Z,0,1)) else p.value=2*pnorm(Z,0,1)
answer=list(Z=Z,pvalue=p.value,AUC=AUC,SE=S,S=S2,R=R)
return(answer)
}
