library(glarma)
library(lmtest)
library(MASS)
library(forecast)

# Get patient names
all_names <- c("HUP064","HUP070","HUP074","HUP075","HUP078","HUP080","HUP086","HUP088","HUP094","HUP105","HUP106","HUP107","HUP111A","HUP116","Study016","Study019","Study020","Study022","Study028","Study029")

# Initialize list
m <- matrix(nrow = length(all_names),ncol = 3)
m_lin <- matrix(nrow = length(all_names),ncol = 3)
m_g <- matrix(nrow = length(all_names),ncol = 3)

out_folder <- "/Users/erinconrad/Desktop/residency stuff/R25/actual work/results/for_r/"

count <- 0
# Loop through all patients
for (name in all_names){
	
	print(name)
	
	count <- count + 1

	MyData <- read.csv(file=paste(out_folder,name,"_SL.csv",sep = ""), header=FALSE, sep=",")

	# The independent and dependent variables that will be used for the GLARMA model. X (independent variable) is a column of ones (for the intercept) and a column with the alpha/delta ratio. Y is a column with spike counts
	x<-matrix(unlist(MyData[,c("V3","V1")]),nc=2)
	colnames(x)<-c("intersect","mean_ad")
	y<-matrix(unlist(MyData[,c("V2")]),nc=1)
	colnames(y)<-c("counts")

	# The variables used for the standard ARIMA model. Y2 is counts. X2 is alpha/delta ratio (no need for column of ones).
	y2<-matrix(unlist(MyData[,c("V2")]),nc=1)
	colnames(y2)<-c("prop")
	x2<-matrix(unlist(MyData[,c("V1")]),nc=1)
	colnames(x2)<-c("b")


	# Do the AR(2) linear regression 
	ar_lin <- Arima(y2,order = c(2,0,0),xreg=x2)
	test_ar <- coeftest(ar_lin)
	b_lin <- test_ar[4,1]
	z_lin <- test_ar[4,3]
	p_lin <- test_ar[4,4]
	
	# add these to the list
	m_lin[count,1] <- b_lin
	m_lin[count,2] <- z_lin
	m_lin[count,3] <- p_lin
	
	# Do the non-autoregressive GLM
	gl <- glm(y~x2,family=poisson)
	gls <- summary(gl)
	b_g <- gls$coefficients[2,1]
	z_g <- gls$coefficients[2,3]
	p_g <- gls$coefficients[2,4]
	
	# Add these to the list
	m_g[count,1] <- b_g
	m_g[count,2] <- z_g
	m_g[count,3] <- p_g


	# Do the GLARMA model. For binomial distribution, logit link function is used.
	model <- try(ar <- glarma(y,x,offset = NULL, type = "Poi",phiLags = c(1,2),method="NR"),silent=TRUE)
	
	if(identical(typeof(model),"character")){
		print("Cannot run NR")
		model <- try(glarma(y,x,offset = NULL, type = "Poi",phiLags = c(1,2),method="FS"),silent=TRUE)
	}
		
	
	if(identical(typeof(model),"character")){
		print("GLM model failed, using poisson model")
		b <- b_g
		z <- z_g
		p <- p_g
	}else{
		ar <- model
		ars <- summary(ar,tests=FALSE)
		b <- ars$coefficients1$Estimate[2]
		z <- ars$coefficients1$z[2]
		p <- ars$coefficients1$Pr[2]
		
		
		ar <- model
		plot(ar,which=1,fits = c(1,3))
		lines(fitted(ar_lin),col="green")
		invisible(readline(prompt="Press [enter] to continue"))
		dev.copy(png,paste(out_folder,name,"SL.png"))
		dev.off()
		graphics.off()
		
	}
	
	

	# Add these to the list
	m[count,1] <- b
	m[count,2] <- z
	m[count,3] <- p

}

# get them back into lists
b_lin <- c(m_lin[,1])
z_lin <- c(m_lin[,2])
p_lin <- c(m_lin[,3])

b <- c(m[,1])
z <- c(m[,2])
p <- c(m[,3])

b_g <- c(m_g[,1])
z_g <- c(m_g[,2])
p_g <- c(m_g[,3])

# Add names
lin_model <- data.frame(all_names,b_lin,z_lin,p_lin)
glarma_model <- data.frame(all_names,b,z,p)
glm_model <- data.frame(all_names,b_g,z_g,p_g)

# output
write.csv(lin_model,paste(out_folder,"lin_SL.csv",sep = ""))
write.csv(glarma_model,paste(out_folder,"glarma_SL.csv",sep = ""))
write.csv(glm_model,paste(out_folder,"glm_SL.csv",sep = ""))

# Find those with p < 0.0025
all_names[which(p < 0.0025)]
