#' Calculate/Simulate model counterfactual in matrix one sector form
#'

#' \code{calc_cf_log_m_1} is a function that calculates the
#' model for one sector assumption for a given set of data and a
#' shock. This function gets sample data in matrix form 
#' ( either arbitrary matrix or real world data)
#' , shocks, and parameters and returns changes of welfare C in two forms:

#' 1. C_hat_1 returns results based on X' computed based on psuedo inverse
#' 2. C_hat_2 returns result based on X' computed based on normalized Q.
#' At the end, C_hat_1 and C_hat_2 should convert to a same number with 
#' with an infinitesimal difference


#' @details This script calculates welfare changes as a result of shock to one sector model:                                   #
#' Here we are using matrix representation of the mode. 

#' using options() function we set digits=22 to make sure about the accuracy of the model
#' of  convergence.

#' The objective of using option function is to make sure
#' about the convergence of X' using psuedo inverse in comparison 
#' with inverse of normalized one using Q matrix which is 
#' Q=gamma*R as it is mentioned in footnote-page 10-task 1. 
#'
#' @param data a list of model variables given as examples
#'   \describe{
#'   \
#'   \item{R}{matrix of country revenues with
#'            \code{J = length(country) columns}}
#'   \item{pi}{vector of import shares across locations for each
#'            destination combination, where the import share of
#'             origin \code{i} products used in destination
#'             \code{j}
#'   \item{gamma}{gamma with items gamma_j vector of intermediate cost shares of each country countries for each
#'                    . after receiving this vector as input
#'                    it is going to be turned into diagonal matrix J*J
#' @param shock a list containing a shock with values
#'   \describe{
#'     \item{T_hat}{technological shock}}
#'     \item{tau_hat}{tariff shocks that may affect both importer and exporter}
#' @return after computing all X,X',A,Ap,Y1,Y1p,Y2,Y2p,B,Bp,Z1,Z1p,Z2,Z2p
#' Z3,Z3p,Z4p,Y3, the function returns C_hat_1 and C_hat_2 based on psuedo and normalized X'



#' calc_cf_log_m -> function(simple_data,shocks,parameters)

#' three blocks of input: 

#' 1. simple_data <-
#' R revenue or GDP vector of size J
#' J number of countries 
#' pi=matrix(c(),nrow=J) , is a vector of size J
#' needed to be divided into intermediate and final pi, pi_I:J*J for intermediate
#' pi_F:J*J for final
#' pi_I=matrix(pi[1:J*J],nrow=J,ncol=J)
#' pi_F=matrix(pi[(J*J)+1:J*J+J*J],nrow=J,ncol=J)
#' gamma=as.vector of size J 



#'2. parameters <-
#'  list(
#'   epsilon = 3,
#'    mobility = "immobile"


#' 3. shock <-
#'    list(
#'      T_hat = matrix(vector(c(), nrow = J,ncol=1) 
#'      tau_hat=matrix(c(),nrow=4J=N,ncol=1) should be divided into two types of shocks: intermediates and final
#'      tau_hat_I=matrix(tau_hat[1:(J*J)],nrow=J,ncol=J)
#'      tau_hat_F=matrix(tau_hat[(J*J)+1:(J*J)+(J*J)],nrow=J,ncol=J)








calc_cf_log_m_1=function(J,R,pi_I,pi_F,gamma,T_hat,tau_hat_I,tau_hat_F,epsilon){



  
  
  ############################################################
  ############### Libraries ##################################
  ############################################################
  
  library(MASS)
  #library(matlib)
  library(pracma)
  library(dplyr)
  library(devtools)
  
  
  options(digits=22) 

  diaggamma=diag(c(gamma));  ###### diagonal matrix of gamma

  IJ=diag(J)                 ###### Identity matrix

  Igammadiag=IJ-diaggamma    ###### Identity-diaggama

  pi_Igamma=crossprod(pi_I,Igammadiag)  ###### pi_I'(identity-diaggama)
  pi_Fgamma=crossprod(pi_F,Igammadiag)  ###### pi_F'(identity-diaggama)
  col=matrix(colSums(pi_I),ncol=2)     ###### vector w elements from pi_I pi11+pi21, pi12+pi22
  colF=matrix(colSums(pi_F),ncol=2)    ###### vector w elements from pi_F pi11+pi21, pi12+pi22

# step 1 Matrix X calc
  X=solve(IJ-pi_Igamma)   ##### matrix X def checked
  ##############[,1]                [,2]
  #[1,] 1.50000000000000022 0.50000000000000011
  #[2,] 0.25000000000000006 1.75000000000000000
  
  
  # step 2 Matrix A calc

  Rrecip=c(sapply(R, function(x) 1/x))  #getting the reciprocal of R vector (1/Ri)


  RR=crossprod(t(Rrecip),R)

  A=(crossprod(t(Rrecip),R)*(pi_I%*%Igammadiag))   #### matrix A def

  ################[,1]                [,2]
  #[1,] 0.29999999999999999 0.15000000000000002
  #[2,] 0.13333333333333333 0.40000000000000002

  # step 3 Matrix diagZ calc
  diagz=diag(c(epsilon*Rrecip*(R%*%(diag(c(col))%*%t(pi_I%*%Igammadiag))))) ###matrix(diag(z))
  ####             [,1]               [,2]
  ##[1,] 1.3500000000000001 0.0000000000000000
  ##[2,] 0.0000000000000000 1.6000000000000001
  
  # step 4 Matrix A' calc
  Ap=(crossprod(t(Rrecip),R)*(pi_F%*%diaggamma))  ##### matrix A'
  #             [,1]                [,2]
  #[1,] 0.25000000000000000 0.30000000000000004
  #[2,] 0.16666666666666666 0.29999999999999999
  
  # step 5 Matrix Y2 calc

  ###### matrix Y2 checked. correct
  Y2=crossprod(pi_F,(diaggamma+crossprod(Igammadiag,crossprod(t(X),crossprod(pi_I,diaggamma)))))
  #                   [,1]   [,2]
  #[1,] 0.43750000000000000 0.5625
  #[2,] 0.37500000000000006 0.6250
  
  
  # step 6 Matrix Z' calc
  #####matrix diag(z')
  diagzp=diag(c(epsilon*Rrecip*(R%*%(diag(c(colF))%*%t(pi_F%*%diaggamma)))))
  ###############[,1]               [,2]
  #[1,] 1.6500000000000001 0.0000000000000000
  #[2,] 0.0000000000000000 1.3999999999999999



  # step 7 Matrix Z2 calc
  ddzx=(diagz*Igammadiag)*(X%*%t(pi_I))
  ddzxp=(diagzp*Igammadiag)*(X%*%t(pi_I))

  zg=diagz*diaggamma
  zIg=diagz*Igammadiag
  zpg=diagzp%*%diaggamma
  zpIg=diagzp%*%Igammadiag

  #######matrix Z2 correct
  Z2=(A-zg)+((epsilon*A-zIg)%*%X%*%t(pi_I)%*%diaggamma) ##checked Matrix Z2

  ####            [,1]                [,2]
  ##[1,] -0.15000000000000008 0.60000000000000009
  ##[2,]  0.43333333333333335 0.10000000000000009


  # step 8 Matrix Z2' calc
  ########matrix Z2' checked correct
  Z2P=(Ap+(epsilon*(Ap%*%Y2))-zpg-(zpIg%*%X%*%t(pi_I)%*%diaggamma))
  ##############[,1]                  [,2]
  #[1,] -0.32187500000000002  0.871874999999999956
  #[2,]  0.54791666666666661 -0.081249999999999933

  # step 9 Matrix Z1 calc
  ########matrix Z1 checked
  Z1=(1/epsilon)*(diagz-((epsilon*A-zIg)%*%X%*%t(pi_I)))
  ############### [,1]                 [,2]
  ##[1,]  0.30000000000000004 -0.29999999999999999
  ##[2,] -0.20000000000000001  0.19999999999999996
  
  
  # step 10 Matrix Y1 calc

  ######## matrix Y1 checked correct
  Y1=((-1)/epsilon)*(t(pi_F)+crossprod(pi_F,crossprod(Igammadiag,crossprod(t(X),t(pi_I)))))
  ###############[,1]                 [,2]
  ## [1,] -0.29166666666666663 -0.37500000000000000
  ## [2,] -0.25000000000000000 -0.41666666666666663
  
  
  # step 11 Matrix Z1' calc
  ######## matrix Z1' Checked
  Z1p=((1/epsilon)*(diagzp+(zpIg%*%X%*%t(pi_I))))+(epsilon*Ap%*%Y1)
  ################[,1]                 [,2]
  ##[1,]  0.38125000000000020 -0.38124999999999998
  ##[2,] -0.25416666666666665  0.25416666666666665

  # step 12 Matrix B calc
  ######### matrix B

  #Check

  B=epsilon*(crossprod(t(Rrecip),R))*((pi_I)%*%(Igammadiag)%*%diag(c(colSums(pi_I))))
  ###
  ##             [,1]                [,2]
  #[1,] 0.89999999999999991 0.45000000000000001
  #[2,] 0.40000000000000002 1.20000000000000018


  # step 13 Matrix B' calc
  ###check Bp

  Bp=epsilon*(crossprod(t(Rrecip),R))*(pi_F%*%diaggamma%*%diag(c(colSums(pi_F))))
  ######[,1]                [,2]
  #[1,] 0.75 0.90000000000000002
  #[2,] 0.50 0.89999999999999991
  
  
  
  # step 14 Matrix Z3 calc
  ######### matrix Z3 checked correct
  Z3=(epsilon*A-zIg)%*%X

  #################[,1]                [,2]
  #[1,] 0.44999999999999990 0.90000000000000002
  #[2,] 0.70000000000000018 0.90000000000000024


  # step 15 Matrix Y3 calc
  ######### matrix Y3 checked correct
  Y3=crossprod(t(pi_Fgamma),X)
  ##################[,1]   [,2]
  #[1,] 0.43750000000000006 0.5625
  #[2,] 0.37500000000000006 0.6250
  
  
  # step 16 Matrix Z3' calc
  ######### matrix Z3'?
  Z3p=(epsilon*(Ap%*%Y3)-(zpIg)%*%X)

  ###[,1]                 [,2]
  ##[1,] -0.57187500000000013  0.57187499999999991
  ##[2,]  0.38124999999999998 -0.38124999999999987
  
  
  
  # step 17 Matrix Z4' calc
  ######### matrix Z4'?
  Z4p=epsilon*Ap
  ##[,1]                [,2]
  ##[1,] 0.75 0.90000000000000013
  ##[2,] 0.50 0.89999999999999991
  
  ######################### comparison of X' using psuedo inverse
  ######################### with normalized format
  # step 18 Matrix X' calc

  ######### matrix X'? using normalized def with Q 
  gR=R*gamma
  Q= matrix(rep(gR,each=J), ncol=J, byrow=FALSE)
  Xp=solve(IJ-(Z2+Z2P)-Q)
  
  
  ######### matrix X'? using psuedo inverse 
  ###############[,1]                 [,2]
  #####[1,] -0.075414012738853453 -0.72458598726114654
  #####[2,] -0.483057324840764357 -0.31694267515923563
  Xpp= ginv(IJ-Z2-Z2P)
  pinv(IJ-Z2-Z2P) #same result as ginv

  ##              [,1]                 [,2]
  ##[1,]  0.23517883390494862 -0.15678588926996570
  ##[2,] -0.23517883390494870  0.15678588926996576

  XZ1=Xp%*%(Z1+Z1p)
  XZ3=Xp%*%(Z3+Z3p)



  ######### result C_hat_1 using X' psuedo inverse checked correct
  C_hat_1=((XZ1-(Y1+(Y2%*%XZ1)))%*%(T_hat))+((Y2%*%Xp-Xp)%*%diag(B%*%t(tau_hat_I)))+((Y2%*%Xp-Xp)%*%diag(Bp%*%t(tau_hat_F)))+
    ((XZ3-(Y3+(Y2%*%XZ3)))%*%(diag(crossprod(pi_I,tau_hat_I))))+(((Xp%*%Z4p)-(IJ+(Y2%*%Xp%*%Z4p)))%*%diag(crossprod(pi_F,tau_hat_F)))

  ######Results of C_hat_2 with normal matrix Xp=solve(IJ-(Z2+Z2P)-Q)

  #################[,1]
  ##[1,] -10.078556263269640
  ##[2,] -14.169851380042463







  XZ1pp=Xpp%*%(Z1+Z1p)
  XZ3pp=Xpp%*%(Z3+Z3p)
  C_hat_2=((XZ1pp-(Y1+(Y2%*%XZ1pp)))%*%(T_hat))+((Y2%*%Xpp-Xpp)%*%diag(B%*%t(tau_hat_I)))+((Y2%*%Xpp-Xpp)%*%diag(Bp%*%t(tau_hat_F)))+
    ((XZ3pp-(Y3+(Y2%*%XZ3pp)))%*%(diag(crossprod(pi_I,tau_hat_I))))+(((Xpp%*%Z4p)-(IJ+(Y2%*%Xpp%*%Z4p)))%*%diag(crossprod(pi_F,tau_hat_F)))



  return(data.frame(C_hat_1,C_hat_2))
  ####################results with pseudo inverse Xpp= ginv(IJ-Z2-Z2P)

  #                 [,1]
  #[1,] -10.078556263269636
  #[2,] -14.169851380042463

}




pi_I=matrix(c(0.6,0.4,0.2,0.8),ncol=2, nrow=2)
pi_F=matrix(c(0.5,0.5,0.4,0.6),nrow=2,ncol=2)






J=2;
gamma=c(0.5,0.5)
epsilon=3;
R=matrix(c(1,1.5),nrow=1)
T_hat=matrix(c(1,2),nrow=2)

tau_hat_I=matrix(c(3,4,5,6),byrow = TRUE, nrow=2)
tau_hat_F=matrix(c(7,8,9,10),byrow = TRUE, nrow=2)


calc_cf_log_m_1(J,R,pi_I,pi_F,gamma,T_hat,tau_hat_I,tau_hat_F,epsilon)

