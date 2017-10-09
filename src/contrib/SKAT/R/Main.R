SKAT_OLD = function(Z,y, X=NULL, kernel = "linear.weighted", out_type="C",method="davies", weights.beta=c(1,25) , weights = NULL, is_intercept = TRUE){

	re<-SKAT_MAIN(Z,y, X, kernel, out_type,method, weights.beta, weights,is_intercept=is_intercept)	
	return(re)
}


SKAT = function(Z,obj,...){

	ml<-match.call()
	ml.name<-names(ml)
	IDX1<-which(ml.name == "y")
	if(length(IDX1) > 0){
		re<-SKAT_MAIN(Z,...)
	} else {

		if(class(obj) == "SKAT_NULL_Model_ADJ"){
			re<-SKAT_With_NullModel_ADJ(Z,obj, ...)
		} else if(class(obj) == "SKAT_NULL_Model"){
			re<-SKAT_With_NullModel(Z,obj, ...)
		} else {
			re<-SKAT_MAIN(Z,obj, ...)
		}

	}	
	class(re)<-"SKAT_OUT"
	return(re)
}

#
#	Check the parameter y and X
#
#
SKAT_MAIN_Check_XY<-function(y,X,is_intercept){

	n = length(y)
	if (!is.null(X)) {
    		if (class(X)!= "matrix") stop("X is not a matrix")
		if(sum((X[,1] - 1)^2,na.rm = TRUE) == 0){
			X1 = X
		} else if(is_intercept){
    			X1 = cbind(1,X)
		} else {
    			X1 = X
		}

  	} else if(is_intercept == TRUE){
    		X1 = as.matrix(rep(1,n))
		
  	} else {
    		stop("Error: Additional covariates (X) are needed when is_intercept = FALSE")
  	}

	if (nrow(X1)!=n) stop("Dimensions of y and X do not match")

	# Check missing of y and X
	id_missy<-which(is.na(y))
	id_missx<-NULL
	for(k in 1:dim(X1)[2]){
		temp1<-which(is.na(X1[,k]))
		id_missx<-union(id_missx,temp1)
	}
	id_missxy<-sort(union(id_missx,id_missy))
	
	id_include<-1:n
	if(length(id_missxy) > 0){
		MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",length(id_missxy))
		warning(MSG,call.=FALSE)
		id_include<-id_include[-id_missxy]
		
	}

	y.test<-y[id_include]
	X1.test<-as.matrix(X1[id_include,])

	return(list(n=n,X1=X1,id_missxy=id_missxy,id_include=id_include,
	y.test=y.test,X1.test=X1.test))

}

#
#	Check the out_type
#
SKAT_MAIN_Check_OutType<-function(out_type){
 	
	if(out_type != "C" && out_type != "D"){
		stop("Invalid out_type!. Please use either \"C\" for the continous outcome or \"D\" for the dichotomous outcome.")
	}

}

#
#	Check the Z, and do imputation
#
#
SKAT_MAIN_Check_Z<-function(Z, n, id_missxy, id_include, SetID, weights, weights.beta, impute.method, is_check_genotype=TRUE){

	#############################################
	# Check parameters

	if (class(Z)!= "matrix") stop("Z is not a matrix")
	if (nrow(Z)!=n) stop("Dimensions of y and Z do not match")
 
	#####################################################
	# Check Z

	if(!is_check_genotype){
		Z.test<-Z[id_include,]
		if(!is.matrix(Z.test)){
			Z.test<-as.matrix(Z.test)
		}
		return(list(Z.test=Z.test,weights=weights, return=0) )
	}

	##############################################
	# Check Missing and doing imputation

	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS]<-NA
	} 

	MAF<-colMeans(Z, na.rm = TRUE)/2
	IDX.Err<-which(MAF > 0.5)	
	if(length(IDX.Err) > 0){
		#msg<-sprintf("Genotypes of some variants are not the number of minor allele! It is fixed!")
		msg<-sprintf("Genotypes of some variants are not the number of minor alleles!")
		warning(msg,call.=FALSE)

		# Fixed by SLEE
		#Z[,IDX.Err]<-2 - Z[,IDX.Err]
		#MAF[IDX.Err]<-1- MAF[IDX.Err]
	}

	###########################################
	# Check non-polymorphic

	if(length(which(MAF > 0)) == 0){
		
		if(is.null(SetID)){
			msg<-sprintf("No polymorphic SNP, so P-value = 1" )
		} else {
			msg<-sprintf("In %s, No polymorphic SNP, so P-value = 1",SetID )
		}
		warning(msg,call.=FALSE)
		re<-list(p.value = 1, p.value.resampling =NA, Test.Type = NA, Q = NA, param=NA, return=1 )   
		return(re)
	}

	##########################################
	# Missing Imputation
	if(length(IDX_MISS) > 0){
		if(is.null(SetID)){
			msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )
		} else {
			msg<-sprintf("In %s, the missing genotype rate is %f. Imputation is applied.", SetID, (length(IDX_MISS))/length(Z) )
		}

		warning(msg,call.=FALSE)
		Z<-Impute(Z,impute.method)
	} 
	
	##########################################
	# Get Weights

	if(is.null(weights)){
		weights<-Beta.Weights(MAF,weights.beta)
	}

	###########################################
	# Check missing of y and X

	if(length(id_missxy) > 0){
		
		if(dim(Z)[2] == 1){
			MAF<-colMeans(cbind(Z[id_include,]))/2
		} else{
			MAF<-colMeans(Z[id_include,])/2
		}
		id_Z<-which(MAF > 0)

		if(length(id_Z) == 0){
			stop("No polymorphic markers!!")
		} else if (length(id_Z) == 1){
			Z<-cbind(Z[,id_Z])
		} else {
			Z<-Z[,id_Z]
		}
		if(!is.null(weights)){
			weights<-weights[id_Z]
		}

	}	
	
	if( dim(Z)[2] == 1){

		if(is.null(SetID)){
			msg<-sprintf("Only one SNP in the SNP set!" )
		} else {
			msg<-sprintf("In %s, Only one SNP in the SNP set!"
			,SetID )
		}
		warning(msg,call.=FALSE)

		Z.test<-as.matrix(Z[id_include,])

	} else {

		Z.test<-Z[id_include,]

	}

	return(list(Z.test=Z.test,weights=weights, return=0) )

}

SKAT_Check_RCorr<-function(kernel, r.corr){

	if(length(r.corr) == 1 && r.corr[1] == 0){
		return(1)
	}
	if(kernel != "linear" && kernel != "linear.weighted"){
		stop("Error: non-zero r.corr only can be used for linear or linear.weighted kernels")
	}

	for(i in 1:length(r.corr)){
		if(r.corr[i] < 0 || r.corr[i] > 1){
			stop("Error: r.corr should be either >= 0 or <= 1")
		}
	}



}

SKAT_Check_Method<-function(method,r.corr){


	if(method != "liu"  && method != "davies" && method != "liu.mod" && method != "optimal" && method != "optimal.moment" && method != "adjust"  ){
		stop("Invalid method!")
	}
	
	if((method == "optimal" || method =="optimal.moment") && length(r.corr) == 1){
		r.corr = (0:10)/10
	}	
	if(method =="optimal"){
		method="davies"
	} else if (method =="optimal.moment") {
		method="liu.mod"
	}

	re<-list(method=method,r.corr=r.corr)
	return(re)

}

#
#	4 methods (keep in mind)
#
SKAT_MAIN = function(Z,y, X=NULL, kernel = "linear.weighted", out_type="C",method="davies", weights.beta=c(1,25) , weights = NULL, 
impute.method = "random", SetID = NULL, is_intercept = TRUE, r.corr=0,
is_check_genotype = TRUE){

	warning("It is old interface!, please run SKAT with SKAT_Null_Model or SKAT_Null_Model_MomentAdjust !")

	out.method<-SKAT_Check_Method(method,r.corr)
	method=out.method$method
	r.corr=out.method$r.corr

	SKAT_Check_RCorr(kernel, r.corr)
	SKAT_MAIN_Check_OutType(out_type)
	out.xy<-SKAT_MAIN_Check_XY(y,X,is_intercept)
	out.z<-SKAT_MAIN_Check_Z(Z,out.xy$n,out.xy$id_missxy
	,out.xy$id_include, SetID,weights,weights.beta, impute.method,is_check_genotype)

	if(out.z$return ==1){
		return(out.z)
	}

	if(out_type == "C"){
		re<-SKAT.linear(out.z$Z.test,out.xy$y.test,out.xy$X1.test,kernel = kernel, weights = out.z$weights, method=method, r.corr=r.corr)
	} else if (out_type == "D"){
		re<-SKAT.logistic(out.z$Z.test,out.xy$y.test,out.xy$X1.test, kernel = kernel, weights = out.z$weights, method=method, r.corr=r.corr)
	}

	return(re)

}



SKAT_With_NullModel = function(Z, obj.res, kernel = "linear.weighted", method="davies", weights.beta=c(1,25), weights = NULL,
impute.method = "random", SetID = NULL, r.corr=0, is_check_genotype=TRUE){

	
	n<-dim(Z)[1]
	m<-dim(Z)[2]

	out.method<-SKAT_Check_Method(method,r.corr)
	method=out.method$method
	r.corr=out.method$r.corr


	SKAT_Check_RCorr(kernel, r.corr)

	out.z<-SKAT_MAIN_Check_Z(Z, n, obj.res$id_missxy
	, obj.res$id_include, SetID, weights, weights.beta, impute.method,
	is_check_genotype)
	if(out.z$return ==1){
		return(out.z)
	}

	if(obj.res$out_type == "C"){
		  if( (kernel =="linear" || kernel == "linear.weighted") && n > m){
		    re = SKAT.linear.Linear(obj.res$res,out.z$Z.test
			,obj.res$X1, kernel, out.z$weights,obj.res$s2,method
			,obj.res$res.out, obj.res$n.Resampling,r.corr=r.corr)
		  } else {  
		    re = SKAT.linear.Other(obj.res$res,out.z$Z.test
			,obj.res$X1, kernel, out.z$weights,obj.res$s2,method
			,obj.res$res.out, obj.res$n.Resampling)  
		  }
	} else if (obj.res$out_type == "D"){

		if( (kernel =="linear" || kernel == "linear.weighted") && n > m){
			re = SKAT.logistic.Linear(obj.res$res, out.z$Z.test
			,obj.res$X1, kernel, out.z$weights, obj.res$pi_1,method
			,obj.res$res.out, obj.res$n.Resampling,r.corr=r.corr)
		} else {  
			re = SKAT.logistic.Other(obj.res$res,out.z$Z.test
			,obj.res$X1, kernel, out.z$weights, obj.res$pi_1, method
			,obj.res$res.out, obj.res$n.Resampling)  
		}
	}

	return(re)

}
 
#
#	Adjustment methods only use liu.mod, so it doesn't need method the "method" parameter
#	I use this field for outcome.type for subfunctions
#	
SKAT_With_NullModel_ADJ = function(Z, obj.res.a, kernel = "linear.weighted", method="adjust", weights.beta=c(1,25), weights = NULL,
impute.method = "random", SetID = NULL, r.corr=0, is_check_genotype=TRUE){

	
	n<-dim(Z)[1]
	m<-dim(Z)[2]
	obj.res<-obj.res.a$re1

	out.method<-SKAT_Check_Method(method,r.corr)
	method=out.method$method
	r.corr=out.method$r.corr

	SKAT_Check_RCorr(kernel, r.corr)
	# Use method field for the type of outcome
	method = obj.res.a$type

	out.z<-SKAT_MAIN_Check_Z(Z, n, obj.res$id_missxy
	, obj.res$id_include, SetID, weights, weights.beta, impute.method, is_check_genotype)
	if(out.z$return ==1){
		return(out.z)
	}

	res2<-NULL
	if(obj.res.a$is_kurtosis_adj){
		res2<-obj.res.a$re2$res.out
	}	

	if(length(r.corr) > 1 && dim(out.z$Z.test)[2] == 1){
		r.corr=0
	}

	if(length(r.corr) == 1){

		re = KMTest.logistic.Linear.VarMatching (obj.res$res,out.z$Z.test
			, obj.res$X1, kernel, out.z$weights, obj.res$pi_1, method
			, obj.res$res.out, obj.res$n.Resampling,r.corr=r.corr
			, obj.res$mu, res.moments=res2)

	} else {

		re = SKAT_Optimal_Logistic_VarMatching(obj.res$res, out.z$Z.test
			, obj.res$X1, kernel, out.z$weights, obj.res$pi_1, method
			, obj.res$res.out, obj.res$n.Resampling, r.corr, obj.res$mu
			, res.moments=res2)

	}
	return(re)


}
 
 
