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
#' @importFrom MASS ginv
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

  # step 2 Matrix A calc

  Rrecip=c(sapply(R, function(x) 1/x))  #getting the reciprocal of R vector (1/Ri)


  RR=crossprod(t(Rrecip),R)

  A=(crossprod(t(Rrecip),R)*(pi_I%*%Igammadiag))   #### matrix A def

  # step 3 Matrix diagZ calc
  diagz=diag(c(epsilon*Rrecip*(R%*%(diag(c(col))%*%t(pi_I%*%Igammadiag))))) ###matrix(diag(z))

  # step 4 Matrix A' calc
  Ap=(crossprod(t(Rrecip),R)*(pi_F%*%diaggamma))  ##### matrix A'

  # step 5 Matrix Y2 calc

  ###### matrix Y2 checked. correct
  Y2=crossprod(pi_F,(diaggamma+crossprod(Igammadiag,crossprod(t(X),crossprod(pi_I,diaggamma)))))

  # step 6 Matrix Z' calc
  #####matrix diag(z')
  diagzp=diag(c(epsilon*Rrecip*(R%*%(diag(c(colF))%*%t(pi_F%*%diaggamma)))))

  # step 7 Matrix Z2 calc
  ddzx=(diagz*Igammadiag)*(X%*%t(pi_I))
  ddzxp=(diagzp*Igammadiag)*(X%*%t(pi_I))

  zg=diagz*diaggamma
  zIg=diagz*Igammadiag
  zpg=diagzp%*%diaggamma
  zpIg=diagzp%*%Igammadiag

  #######matrix Z2 correct
  Z2=(A-zg)+((epsilon*A-zIg)%*%X%*%t(pi_I)%*%diaggamma) ##checked Matrix Z2

  # step 8 Matrix Z2' calc
  ########matrix Z2' checked correct
  Z2P=(Ap+(epsilon*(Ap%*%Y2))-zpg-(zpIg%*%X%*%t(pi_I)%*%diaggamma))

  # step 9 Matrix Z1 calc
  ########matrix Z1 checked
  Z1=(1/epsilon)*(diagz-((epsilon*A-zIg)%*%X%*%t(pi_I)))

  # step 10 Matrix Y1 calc

  ######## matrix Y1 checked correct
  Y1=((-1)/epsilon)*(t(pi_F)+crossprod(pi_F,crossprod(Igammadiag,crossprod(t(X),t(pi_I)))))

  # step 11 Matrix Z1' calc
  ######## matrix Z1' Checked
  Z1p=((1/epsilon)*(diagzp+(zpIg%*%X%*%t(pi_I))))+(epsilon*Ap%*%Y1)

  # step 12 Matrix B calc
  ######### matrix B

  #Check

  B=epsilon*(crossprod(t(Rrecip),R))*((pi_I)%*%(Igammadiag)%*%diag(c(colSums(pi_I))))

  # step 13 Matrix B' calc
  ###check Bp

  Bp=epsilon*(crossprod(t(Rrecip),R))*(pi_F%*%diaggamma%*%diag(c(colSums(pi_F))))

  # step 14 Matrix Z3 calc
  ######### matrix Z3 checked correct
  Z3=(epsilon*A-zIg)%*%X

  # step 15 Matrix Y3 calc
  ######### matrix Y3 checked correct
  Y3=crossprod(t(pi_Fgamma),X)

  # step 16 Matrix Z3' calc
  ######### matrix Z3'?
  Z3p=(epsilon*(Ap%*%Y3)-(zpIg)%*%X)

  # step 17 Matrix Z4' calc
  ######### matrix Z4'?
  Z4p=epsilon*Ap

  ######################### comparison of X' using psuedo inverse
  ######################### with normalized format
  # step 18 Matrix X' calc

  ######### matrix X'? using normalized def with Q
  gR=R*gamma
  Q= matrix(rep(gR,each=J), ncol=J, byrow=FALSE)
  Xp=solve(IJ-(Z2+Z2P)-Q)


  ######### matrix X'? using psuedo inverse

  Xpp= ginv(IJ-Z2-Z2P)
# IS THIS NECESSARY:  pinv(IJ-Z2-Z2P) #same result as ginv

  XZ1=Xp%*%(Z1+Z1p)
  XZ3=Xp%*%(Z3+Z3p)

  ######### result C_hat_1 using X' psuedo inverse checked correct
  C_hat_1=((XZ1-(Y1+(Y2%*%XZ1)))%*%(T_hat))+((Y2%*%Xp-Xp)%*%diag(B%*%t(tau_hat_I)))+((Y2%*%Xp-Xp)%*%diag(Bp%*%t(tau_hat_F)))+
    ((XZ3-(Y3+(Y2%*%XZ3)))%*%(diag(crossprod(pi_I,tau_hat_I))))+(((Xp%*%Z4p)-(IJ+(Y2%*%Xp%*%Z4p)))%*%diag(crossprod(pi_F,tau_hat_F)))

  XZ1pp=Xpp%*%(Z1+Z1p)
  XZ3pp=Xpp%*%(Z3+Z3p)
  C_hat_2=((XZ1pp-(Y1+(Y2%*%XZ1pp)))%*%(T_hat))+((Y2%*%Xpp-Xpp)%*%diag(B%*%t(tau_hat_I)))+((Y2%*%Xpp-Xpp)%*%diag(Bp%*%t(tau_hat_F)))+
    ((XZ3pp-(Y3+(Y2%*%XZ3pp)))%*%(diag(crossprod(pi_I,tau_hat_I))))+(((Xpp%*%Z4p)-(IJ+(Y2%*%Xpp%*%Z4p)))%*%diag(crossprod(pi_F,tau_hat_F)))

  data.frame(C_hat_1,C_hat_2)
}
