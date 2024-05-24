# load the ncdf4 package
library(ncdf4)


# set path and filename
ncname <- "SOCAL_5PC_PERTURBED_ELEMGFS_updated4"  
ncfname <- paste(ncname, ".nc4", sep="")
dname <- "SoCal"  
dname2 <- "obsdata"
dname3 <- "Gtensors"

# open a netCDF file
ncfname <- "C:/Users/johns/Documents/Project Files/SOCAL_5PC_PERTURBED_ELEMGFS_updated4.nc4"
ncin <- nc_open(ncfname)
print(ncin)


###################################notes
#	15 global attributes:
#        nm: 300
#        nm_note: Number of perturbed 1D Earth models
#        ns: 6
#        ns_note: Number of stations
#        nc: 3
#        nc_note: Number of components (Vertical up, North, East)
#        ne: 6
#        ne_note: Number of elementary Green's function: (m_nn, m_ee, m_dd, m_ne, m_nd, m_ed)
#        nt: 180
#        nt_note: Number of point in time; seismograms are samped at 1 Hz
#        note1: You might notice there is different in the convention of channel (up, north, east)and the moment tensor components (north, east, down), but rest assured that they are checked.The MT solution you obtain will be in the (north, east, down) coordinates.
#        note2: If your moment tensor, m, is (1xne) vector, then pred = m . Gtensors, where . is a matrix multiplication.
#        SoCal_note: Green function corresponding to the original SoCal model
#        obsdata_note: Observed data
#        Gtensors_note: Green function of 300 perturbed Earth models








#########################Greens functions covariates


# get SoCal
elem_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(elem_array)
elem_array[[1,1,1,2]]



dimnames(elem_array)[[1]] <- 1:180
dimnames(elem_array)[[2]] <- c("m_nn","m_ee","m_dd","m_ne","m_nd","m_ed")
dimnames(elem_array)[[3]] <- c("Vertical_up","North","East")
dimnames(elem_array)[[4]] <- c("BK.BKS","BK.CMB","BK.KCC","BK.ORV","BK.PKD","US.MNV")


my_data <- as.data.frame.table(elem_array) 
#my_data2 <- as.data.frame(elem_array) 
names(my_data)=c("time","element","channel","station","value")
my_data[1:10,]
#my_data[1:370,]


# get a single slice or layer 
m <- 1
g_slice1 <- elem_array[,m,,]
dim(g_slice1)

# get a single slice or layer 
m <- 2
g_slice2 <- elem_array[,m,,]
dim(g_slice2)

# get a single slice or layer 
m <- 3
g_slice3 <- elem_array[,m,,]
dim(g_slice3)

# get a single slice or layer 
m <- 4
g_slice4 <- elem_array[,m,,]
dim(g_slice4)


# get a single slice or layer 
m <- 5
g_slice5 <- elem_array[,m,,]
dim(g_slice5)

# get a single slice or layer 
m <- 6
g_slice6 <- elem_array[,m,,]
dim(g_slice6)




my_data_slice1 <- as.data.frame.table(g_slice1) 
names(my_data_slice1)=c("time","channel","station","value")
my_data_slice1[1:10,]

my_data_slice2 <- as.data.frame.table(g_slice2) 
names(my_data_slice2)=c("time","channel","station","value")
my_data_slice2[1:10,]

my_data_slice3 <- as.data.frame.table(g_slice3) 
names(my_data_slice3)=c("time","channel","station","value")
my_data_slice3[1:10,]

my_data_slice4 <- as.data.frame.table(g_slice4) 
names(my_data_slice4)=c("time","channel","station","value")
my_data_slice4[1:10,]


my_data_slice5 <- as.data.frame.table(g_slice5) 
names(my_data_slice5)=c("time","channel","station","value")
my_data_slice5[1:10,]

my_data_slice6 <- as.data.frame.table(g_slice6) 
names(my_data_slice6)=c("time","channel","station","value")
my_data_slice6[1:10,]


#########################observed data

# get obsdata
elem_array <- ncvar_get(ncin,dname2)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(elem_array)
elem_array[[1,1,1]]


dimnames(elem_array)[[1]] <- 1:180
dimnames(elem_array)[[2]] <- c("Vertical_up","North","East")
dimnames(elem_array)[[3]] <- c("BK.BKS","BK.CMB","BK.KCC","BK.ORV","BK.PKD","US.MNV")




my_data2 <- as.data.frame.table(elem_array) 
names(my_data2)=c("time","channel","station","value")
my_data2[1:10,]


#merge slice with observed data
g1=my_data_slice1[,"value"]
g2=my_data_slice2[,"value"]
g3=my_data_slice3[,"value"]
g4=my_data_slice4[,"value"]
g5=my_data_slice5[,"value"]
g6=my_data_slice6[,"value"]
merge1=cbind(my_data2,g1,g2,g3,g4,g5,g6)
merge1[1:10,]


write.table(merge1,"merge1.csv",sep=",",row.names=FALSE,col.names=TRUE)




#crude least squares fit
mm=lm(merge1[,"value"]~-1 + merge1[,"g1"]+merge1[,"g2"]+merge1[,"g3"]+merge1[,"g4"]+merge1[,"g5"]+merge1[,"g6"])
mm0=mm
beta=matrix(0,6,1)
beta[,1]=mm0$coefficients
sigma=summary(mm)$sigma
summary(mm)
signif(beta,5)
signif(sigma^2,5)

#########shift covariates 

m_group=18
n_time=180
n_tot=m_group*n_time
x_matrix=as.matrix(merge1[,5:10])


#start with one single a
seta=matrix(0,1,11)
seta[1]=-5
seta[2]=-4
seta[3]=-3
seta[4]=-2
seta[5]=-1
seta[6]=0
seta[7]=1
seta[8]=2
seta[9]=3
seta[10]=4
seta[11]=5

#initial value for tau
tau=0

#needed for weights calculation
LL=3
KL=180/LL

#needed for lme, the groups for random effects
N_big=n_tot*11
x_comb=matrix(0,N_big,1)
gro=1
ind=0
for (j in 1:N_big)
{	
	ind=ind+1
	if(ind==sum(LL,1)){ind=1;gro=gro+1}
	x_comb[j]=gro


}




for (ii in 1:30)
{
#####################repeat across all a's

	#needed for weights calculation
	J_L=matrix(tau^2,LL,LL)
	V_i=J_L+diag(sigma^2,LL)
	V_inverse=solve(V_i)
	V_big=V_inverse
	library(Matrix)
	for (j in 1:sum(KL,-1))
	{
		V_big=bdiag(V_big,V_inverse)
	}
	#on log scale
	detV=KL*determinant(V_i, logarithm = TRUE)$modulus[1]

	p_group=matrix(0,m_group,11)

	x_big=matrix(0,1,6)
	y_big=matrix(0,1,1)

	for (kk in 1:11)
	{

		aa=seta[kk]

		# 0 or 1 or 2 (0=nothing, 1=shiftforward, 2=shiftback)
		if (aa < 0){backforward=matrix(2,m_group,1)}
		if (aa > 0){backforward=matrix(1,m_group,1)}
		if (aa == 0) {backforward=matrix(0,m_group,1)}

		#shift amount
		shift=matrix(abs(aa),m_group,1)


		groups=matrix(0,m_group,2)
		counter=0
		for (i in 1:m_group)
		{

			groups[i,1]=counter+1
			groups[i,2]=i*n_time
			counter=counter+n_time
		}


		x_switch=matrix(0,m_group*n_time,6)
		for (i in 1:m_group)
		{

			if (backforward[i]==0)
			{
				x_switch[groups[i,1]:groups[i,2],]=x_matrix[groups[i,1]:groups[i,2],]
			}
	
			if (backforward[i]==1)
			{
				position1=groups[i,1]+shift[i]
				position2=groups[i,2]
				piece1=x_matrix[position1:position2,]
				position1=groups[i,1]
				position2=groups[i,1]+shift[i]-1
				piece2=x_matrix[position1:position2,]

				x_switch[groups[i,1]:groups[i,2],]=rbind(piece1,piece2)
			}

			if (backforward[i]==2)
			{
				position1=groups[i,2]-shift[i]
				position2=groups[i,2]
				piece1=x_matrix[position1:position2,]
				position1=groups[i,1]
				position2=groups[i,2]-shift[i]-1
				piece2=x_matrix[position1:position2,]

				x_switch[groups[i,1]:groups[i,2],]=rbind(piece1,piece2)
			}
		}

		resp=merge1[,"value"]

	

	
		for (i in 1:m_group)
		{
			p_group[i,kk]=(-0.5*n_time*log(2*pi)-0.5*detV-0.5*(t(resp[groups[i,1]:groups[i,2]])-t(x_switch[groups[i,1]:groups[i,2],1:6]%*%beta))%*%V_big%*%t(t(resp[groups[i,1]:groups[i,2]])-t(x_switch[groups[i,1]:groups[i,2],1:6]%*%beta)))[1,1]
			#sum(prob[groups[i,1]:groups[i,2]])
		}

		# the above is on the log scale

		x_big=rbind(x_big,x_switch)
		y_big=rbind(y_big,t(t(resp)))

	}
	x_big=x_big[2:length(x_big[,1]),]
	y_big=y_big[2:length(y_big[,1]),]

	for (i in 1:m_group)
	{
		p_group[i,]=exp(p_group[i,]-max(p_group[i,]))/sum(exp(p_group[i,]-max(p_group[i,])))
	}

	weights1=matrix(0,n_tot,11)
	for (kk in 1:11)
	{
		for (i in 1:m_group)
		{
			weights1[groups[i,1]:groups[i,2],kk]=p_group[i,kk]
		}
	}

	#stack weights but we need to repeat within. 
	weights2=as.matrix(weights1[,1])
	for(kk in 2:11)
	{
	
		weights2=rbind(weights2,as.matrix(weights1[,kk]))
	}

	#multiply weights by 11
	weights2=weights2*11

	library(nlme)
	data_big=cbind(y_big,x_big,weights2,x_comb)
	data_big=as.data.frame(data_big)
	names(data_big)=c("resp","x1","x2","x3","x4","x5","x6","weights2","xcomb")
	data_big[,"groups"]=1:length(data_big[,1])
	data_big <- groupedData(resp ~ 1 | xcomb, data_big)
	data_big[,"weights2"]=data_big[,"weights2"]+0.0000000000000001

	#weighted least squares fit
	#mm=gls(resp ~ -1+x1+x2+x3+x4+x5+x6,data=data_big,weights=~1/weights2,method="ML")
	#mm1=mm
	#beta=matrix(0,6,1)
	#beta[,1]=mm1$coefficients
	#sigma=mm1$sigma
	#print(beta)
	#print(sigma^2)


	#generalised least squares fit
	mm=gls(resp ~ -1+x1+x2+x3+x4+x5+x6,data=data_big,correlation=corCompSymm(value=0.3,form=~1|xcomb),weights=~1/weights2,method="ML")
	mm3=mm
	beta=matrix(0,6,1)
	beta[,1]=mm3$coefficients
	rho_star=intervals(mm3)$corStruct[2]
	sig_star=mm3$sigma^2
	tau=sqrt(rho_star*sig_star)
	sigma=sqrt(sig_star*(1-rho_star))
	print(beta)
	print(sigma^2)
	print(tau^2)



}







round(p_group,2)


#the posteriors are pretty close to binary within rows, so we can just choose 
#one timeshift and refit to get SE's (uncertainties).

round(beta,4)
sigma^2
tau^2




##############################################
####identify shifted green's functions
##############################################

tot_switch=11*m_group
groups_switch=matrix(0,tot_switch,4)
counter=0
groups_switch[,2]=rep(1:18,11)
counter_2=0
counter_3=0
for (i in 1:tot_switch)
{
	groups_switch[i,3]=counter+1
	groups_switch[i,4]=i*n_time
	counter=counter+n_time
	groups_switch[i,1]=counter_2+1
	counter_3=counter_3+1
	if (counter_3==18){counter_2=counter_2+1;counter_3=0}
}


x_back=matrix(0,1,6)
p_group=round(p_group)


kk2=matrix(6,m_group,1)
for (i in 1:m_group)
{

	for (j in 1:11)
	{
		if (p_group[i,j]==1){kk2[i]=j}
	}

}

for (j in 1:m_group)
{
	for (i in 1:prod(m_group*11))
	{

		if (groups_switch[i,1]==kk2[j] & groups_switch[i,2]==j){x_back=rbind(x_back,x_big[groups_switch[i,3]:groups_switch[i,4],])}

	}
}
x_back=x_back[2:length(x_back[,1]),]



##############################################
####fit lm model with no weights
##############################################

mm=lm(merge1[,"value"]~-1 + x_back[,1]+x_back[,2]+x_back[,3]+x_back[,4]+x_back[,5]+x_back[,6])
mm1=mm
summary(mm1)$sigma^2
sum((merge1[,"value"]-mm1$fitted.values)^2)/3240

##############################################
####fit lme model with no weights
##############################################

 
#lme fit
library(nlme)
data_big=cbind(merge1[,"value"],x_back,x_comb[1:n_tot])
data_big=as.data.frame(data_big)
names(data_big)=c("resp","x1","x2","x3","x4","x5","x6","xcomb")
data_big <- groupedData(resp ~ 1 | xcomb, data_big)
mm3=lme(resp ~ -1+x1+x2+x3+x4+x5+x6,random=pdIdent(~ 1), data=data_big,control=list(maxIter=1000,msMaxIter=1000),method="ML")
beta[,1]=fixef(mm3)
tau=as.numeric(VarCorr(mm3)[1,2])
sigma=as.numeric(VarCorr(mm3)[2,2])
print(beta)
print(sigma^2)
print(tau^2)

#SE calculation
J_L=matrix(tau^2,LL,LL)
V_i=J_L+diag(sigma^2,LL)
V_inverse=solve(V_i)
V_big=V_inverse
library(Matrix)
for (j in 1:sum(KL,-1))
{
	V_big=bdiag(V_big,V_inverse)
}

V_big2=V_big
for (j in 1:sum(m_group,-1))
{
	V_big2=bdiag(V_big2,V_big)
}

summary(mm3)
#SE's
SE=solve(t(x_back)%*%V_big2%*%x_back)

round(SE,5)

sqrt(diag(SE))

2.1735 + 2*sqrt(0.04462)
2.1735 - 2*sqrt(0.04462)
2*sqrt(0.04462)


############################
######plots
#########################

pred1=predict(mm3)

#i=10
#channel: Vertical_up
#station: BK.ORV

i=10
plot(merge1[groups[i,1]:groups[i,2],"time"],merge1[groups[i,1]:groups[i,2],"value"],ylab="displacment",xlab="time",pch=16,cex=2,col="black")
lines(merge1[groups[i,1]:groups[i,2],"time"],as.matrix(mm0$fitted.values[groups[i,1]:groups[i,2]]),col="blue",pch=16,cex=2)
lines(merge1[groups[i,1]:groups[i,2],"time"],as.matrix(mm1$fitted.values[groups[i,1]:groups[i,2]]),col="red",pch=16,cex=2)
lines(merge1[groups[i,1]:groups[i,2],"time"],as.matrix(pred1[groups[i,1]:groups[i,2]]),col="green",pch=16,cex=2)

legend("topleft", legend=c("observed", "least squares","timeshift","timeshift+amplitude"),  col = c("black","blue","red","green"),lty = c(NA,1,1,1),cex=1,pch = c(16, NA, NA, NA) )

#plot(merge1[groups[i,1]:groups[i,2],"time"],merge1[groups[i,1]:groups[i,2],"value"]-as.matrix(pred1[groups[i,1]:groups[i,2]]))

