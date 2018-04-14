# Öknonometrie
function(){
  # Eigenschaften der Gleichung
  #   Invertible if: |b1/a1|<1                             MA nur invertierbar, wenn Gewichtungen gegen 0 konvergiert
  #                  abs(eigen(Bhat_mat)$values) < 1       AR nur Invertierbar, wenn station??r
  abs(eigen(A_mat)$values)
  
  #   stationary if: MA (MA is always stationary)
  #                  AR (abs(eigen(Ahat_mat)$values) if all Eigenvalues (or only a1) are smaller |1|)
  #                  AR ist NICHT station??r, wenn abh??ngig von der Zeit (t)
  
  # MA Prozess
  B_mat <- function(x){
    if(length(x)==1){
      B_matt<-rbind(c(-b_1,1),c(0,0))
    }else{
      a_vec<-t(x)
      B_matt<-matrix(ncol=length(a_vec),nrow=length(a_vec))
      B_matt[1,]<--a_vec
      B_matt[2:length(a_vec),1:(length(a_vec)-1)]<-diag(rep(1,length(a_vec)-1))
      B_matt[2:length(a_vec),length(a_vec)]<-0
    }
    return(B_matt)
  }
  # (B_matARMA <- t(B_mat(c(b_vec,0))))
  t((Bhat_mat <- B_mat(a_vec)))
  abs(eigen(Bhat_mat)$values)
  b_vec<-c(0.2,0.1,-0.3,0,0.1)
  
  # AR Prozess
  A_mat <- function(x){
    if(length(x)==1){
      A_matt<-rbind(c(a_1,1),c(0,0))
    }else{
      a_vec<-t(x)
      A_matt<-matrix(ncol=length(a_vec),nrow=length(a_vec))
      A_matt[1,]<-a_vec
      A_matt[2:length(a_vec),1:(length(a_vec)-1)]<-diag(rep(1,length(a_vec)-1))
      A_matt[2:length(a_vec),length(a_vec)]<-0
    }
    return(A_matt)
  }
  # (A_matARMA <- t(A_mat(c(a_vec,0))))
  (Ahat_mat <- A_mat(a_vec))
  abs(eigen(Ahat_mat)$values)
  a_vec<-c(0.2,0.1,-0.3,0,0.1)
  
  
  # Compute weights of AR/MA-representation
  # Mit welchem AR prozess l??sst sich dieser MA Prozess beschreiben?
  # dort wo die punkte konstant bleiben wird der Wert genommen => AR(30)
  
  ## Gegeben MA
  MAtoARweights <- function(x,XdimPlot){
    B_mat <- function(x){
      a_vec<-t(x)
      A_mat<-matrix(ncol=length(a_vec),nrow=length(a_vec))
      A_mat[1,]<--a_vec
      A_mat[2:length(a_vec),1:(length(a_vec)-1)]<-diag(rep(1,length(a_vec)-1))
      A_mat[2:length(a_vec),length(a_vec)]<-0
      return(A_mat)
    }
    (Bhat_mat <- B_mat(x))
    B_k<-diag(rep(1,dim(Bhat_mat)[1]))
    Bk_11<-0:XdimPlot
    for (i in 0:XdimPlot)
    {
      Bk_11[i+1]<-B_k[1,1]
      B_k<-B_k%*%Bhat_mat
    }
    plot(Bk_11)
    return(Bk_11)
  }
  MAtoARweights(b_vec,200) # falls nur ein Wert c(b_vec,0)
  
  ## Gegeben AR
  ARtoMAweights <- function(a_vec,XdimPlot){
    A_mat <- function(x){
      a_vec<-t(x)
      A_mat<-matrix(ncol=length(a_vec),nrow=length(a_vec))
      A_mat[1,]<-+a_vec
      A_mat[2:length(a_vec),1:(length(a_vec)-1)]<-diag(rep(1,length(a_vec)-1))
      A_mat[2:length(a_vec),length(a_vec)]<-0
      return(A_mat)
    }
    (Ahat_mat <- A_mat(x))
    B_k<-diag(rep(1,dim(Ahat_mat)[1]))
    Bk_11<-0:XdimPlot
    for (i in 0:XdimPlot)
    {
      Bk_11[i+1]<-B_k[1,1]
      B_k<-B_k%*%Ahat_mat
    }
    plot(Bk_11)
    return(Bk_11)
  }
  ARtoMAweights(a_vec,200) # falls nur ein Wert c(a_vec,0)
  
  # ARMA-Inversion mit nur 1 Koeff
  function(){
    len<-100
    a_1<-0.9
    b_1<-0.9
    (A_matt <- A_mat(a_1))
    (B_matt <- B_mat(b_1))
    
    len<-100
    A_k<-diag(rep(1,2))
    B_k<-A_k
    a_vec<-c(1,-a_1)
    b_vec<-c(1,b_1)
    AR_weights<-1:len
    MA_weights<-AR_weights
    for (i in 1:len)
    {
      AR_weights[i]<-(B_k%*%a_vec)[1]
      MA_weights[i]<-(A_k%*%b_vec)[1]
      A_k<-A_k%*%A_matt
      B_k<-B_k%*%B_matt
    }
    ts.plot(AR_weights,col="red",xlab="",ylab="",main="AR- and MA-inversions of ARMA(1,1)",ylim=c(-1,1))
    lines(MA_weights,col="blue")
    mtext("MA-inversion", side = 3, line = -1,at=len/2,col="blue")
    mtext("AR_inversion", side = 3, line = -2,at=len/2,col="red")
  }
  
  # ab welchem h0 wird der h-Schritt gr??sser als  (bsp. 90%) akkurat.
  MA_weights # wird ben??tigt
  accuracy <- 0.9 # 90% Genauigkeit
  which(sapply(1:20, FUN = function(i) sum(MA_weights[1:i]^2)/sum(MA_weights^2))>accuracy,arr.ind = T)[1]
  
  # weitere m??glichkeit
  ARMAtoMA(ar = a_vec, ma = b_vec, lag.max=100) # ARMAtoAR, wenn a_vec und b_vec vertauscht sind
  ARMAacf(ar = a_vec, ma = b_vec, lag.max = 10, pacf = FALSE)
  
  polyroot(c(a_vec)) # ergibt immer (1-erg)*(1-erg)
  
  
  # Erwartungswert eines Prozesses berechnen (NUR AR TEIL, MA weglassen!)
  mu<-Konstantte/(Standartabweichung- sum(c(a1,a2,a3))) # kein MA
  
  
  # acf zum sch??tzen von MA() => die wievielte ist deutlich ausserhalb des Bereichs (1. Strich ist der 0. nicht 1.)
  # x ist ein Datensatz
  acf(x)                    # beginnt bei 0
  # partieller acf
  acf(x,type="partial")
  ?pacf()                   # beginnt bei 1!
  ARMAacf(ar = a_vec, ma = b_vec, lag.max = 10, pacf = FALSE)
  
  # ein Modell sch??tzen
  #   stationarit??t (p,1,q) => 1. Differenzen (maximal 2)
  x_obj <- arima(x,order = c(AR-Anteil = 0, station??r = 0, MA-Anteil = 2), method="CSS")
  arima(x,order=c(0,0,0),method="CSS")
  # Saisonelle Modelle berechnen SARMA
  arima(x, include.mean=FALSE, order=c(0,1,1), method="ML",seasonal=list(order=c(0,1,1),period=12))
  
  # Sch??tzen der Koeffizienten
  x_obj <- arima(x[1:100], order=c(3,0,3))
  
  x_obj$coef # Koeffizienten der Gleichung anzeigen
  summary(x_obj)
  
  # Qualit??t betrachten
  tsdiag(x_obj, gof.lag=10)
  
  # Werte in Zukunft sch??tzen ($se = Standarterror)
  pval <- predict(x_obj, n.ahead=10)
  pval$pred-2*pval$se
  
  # Simulieren von Daten
  # Xt= epsilon[t] + 0.5 * epsilon[t-5]
  # [t-5] => Vektor l??nge = 4  &&  0.5 = letze der 4 Zahlen ist 0.5
  #sigma <- 1; length <- 100
  arima.sim(n = length, list(ar=c(0,0),ma = c(0,0,0,0,0.5)), sd = sigma)
  arima.sim
  
  #Zeitreihen Plotten
  ts.plot(x , lty=1:2)
  
  # Bei Station??ren AR()
  # R(0) = var(x) berechnen
  R0 <- function(a_vec,sigma){
    # sigma<-1
    # a_vec<-c(1.6,-1.4,0.5)
    A_mat<-matrix(ncol=length(a_vec),nrow=length(a_vec))
    A_mat[1,]<-a_vec
    A_mat[2:length(a_vec),1:(length(a_vec)-1)]<-diag(rep(1,length(a_vec)-1))
    A_mat[2:length(a_vec),length(a_vec)]<-0
    A_mat
    
    AkA<-kronecker(A_mat,A_mat)
    AkA
    vec_sigma<-c(sigma^2,rep(0,dim(AkA)[2]-1))
    vec_sigma
    # Ricatti equation
    R_mat_0_vec<-solve(diag(rep(1,dim(AkA)[2]))-AkA)%*%vec_sigma
    R_mat_0<-matrix(R_mat_0_vec,ncol=dim(A_mat)[2])
    return(R_mat_0[1,1])
  }
  
  # R(k) = cov(x_n , x_n-1) berechnen (mit Hilfe von R(0)!)
  # Auto-Covarianzen zu Verz??gerung k
  Rk <- function(a_vec,sigma,k){
    # R0 berechnen??
    A_mat <- function(x){
      a_vec<-t(x)
      A_mat<-matrix(ncol=length(a_vec),nrow=length(a_vec))
      A_mat[1,]<-a_vec
      A_mat[2:length(a_vec),1:(length(a_vec)-1)]<-diag(rep(1,length(a_vec)-1))
      A_mat[2:length(a_vec),length(a_vec)]<-0
      return(A_mat)
    }
    (Ahat_mat <- A_mat(a_vec))
    R1 <- function(a_vec,sigma,Ahat_mat){
      # Exercise 1.a
      # sigma<-1
      # a_vec<-c(1.6,-1.4,0.5)
      
      AkA<-kronecker(Ahat_mat,Ahat_mat)
      AkA
      vec_sigma<-c(sigma^2,rep(0,dim(AkA)[2]-1))
      vec_sigma
      # Ricatti equation
      R_mat_0_vec<-solve(diag(rep(1,dim(AkA)[2]))-AkA)%*%vec_sigma
      R_mat_0<-matrix(R_mat_0_vec,ncol=dim(Ahat_mat)[2])
      return(R_mat_0)
    }
    len_R<-k
    R_mat_k<-R1(a_vec,sigma,Ahat_mat)
    R_k<-rep(0,len_R)
    for (i in 0:(len_R-1))
    {
      R_k[i+1]<-R_mat_k[1,1]
      R_mat_k<-Ahat_mat%*%R_mat_k
    }
    return(R_k[k])
  }
  
  R0(c(1.6,-1.4,0.5), 1)
  Rk(c(1.6,-1.4,0.5), 1, 2) # K-1 = Rk !!!!!
  
  # Plotten der R(k) funktion  
  function(){
    par(mfrow=c(1,1))
    # Exercise 1.b
    # a_vec<-c(0.2,0.1,-0.3,0,0.1)
    len<-95
    x<-arima.sim(n=len,list(ar=a_vec))
    
    # Exercise 1.c
    
    acf_true<-R_k/R_mat_0[1,1]
    
    ts.plot(acf(x,plot=F,lag.max=100)$acf,xlab="",ylab="Time")
    lines(acf_true,col="blue")
    lines(rep(2/sqrt(len),len),lty=2)
    lines(rep(-2/sqrt(len),len),lty=2)
    lines(rep(0,len))
    mtext("Sample acf", side = 3, line = -1,at=length(acf(x,plot=F)$acf)/2,col="black")
    mtext("True acf", side = 3, line = -2,at=length(acf(x,plot=F)$acf)/2,col="blue")
  }
  
  
  # Absch??tzen der Anzahl Parameter (sobald gen??gend nah bei 0 => gen??gend Parameter)
  # Bsp.
  #  Plot2order <- 2
  #  maxorder <- 9
  #  y <- arima.sim(list(ma=0.1,ar=c(0,0,0.5)),n=100)
  AIC1 <- function(y,maxorder,Plot2order){
    
    sigma_k<-rep(0,maxorder+1)
    sigma_k[1]<-var(y)
    
    for (k in 1:maxorder) #k<-10
    {
      AR_obj<-arima(y,order=c(k,0,0))
      sigma_k[k+1]<-mean(AR_obj$res^2)
      # This is the same as
      sigma_k[k+1]<-AR_obj$sigma
      print(paste("AR-Order=",k," AIC=",AR_obj$aic,sep=""))
    }
    
    log_sigma_k <- log(sigma_k)
    par(mfrow=c(2,1))
    anf<-1
    plot((log_sigma_k)[anf:(maxorder+1)],col="blue",xlab="AR-order",ylab="",
         main="Sample mean-square error of AR(2)",type="l",axes="F")
    axis(1,at=anf:(maxorder+1),labels=-1+anf:(maxorder+1))
    axis(2)
    box()
    anf<-Plot2order+1
    plot((log_sigma_k)[anf:(maxorder+1)],col="blue",xlab="AR-order",ylab="",
         main="Same as above for orders=2,3,...,10 (larger or equal than true order 2)",type="l",axes="F")
    axis(1,at=1:(maxorder-anf+2),labels=-1+anf:(maxorder+1))
    axis(2)
    box()
  }
  AIC1(y,9,3)
  
  # ARMAorder
  function(){
    # Ben??tigte Werte (Eingabe in Funktion)
    y                   # Datenreihe
    maxarorder <- 5     #
    maxmaorder <- 5     #
    d <- 0              # Stationarit??t I(0)
    
    # Funktion, Schritt f??r Schritt ausf??hren
    sigma_jk<-matrix(rep(0,(maxmaorder+1)*(maxarorder+1)),
                     ncol=maxmaorder+1,nrow=maxarorder+1)
    sigma_jk[1,1]<-var(y)
    for (j in 0:maxarorder)
    {
      
      for (k in 0:maxmaorder)         #k<-4  j<-3
      {
        # Avoid interrupt if optimization fails    
        try_obj<-try(arima(y,order=c(j,d,k)),silent=T)
        if (is.character(try_obj[1]))
        {
          # If interrupted: sigma is set to a large value (`bad' model)      
          sigma_jk[j+1,k+1]<-1.e+99      
        } else
        {
          ARMA_obj<-arima(y,order=c(j,d,k))
          print(c(j,k))
          sigma_jk[j+1,k+1]<-ARMA_obj$sigma
          print(paste("  AR-order=",j,"MA-Order=",k,"  AIC=",ARMA_obj$aic,sep=""))      
        }
      }
    }
    log_sigma_jk<-log(sigma_jk)
    aic<-sigma_jk
    dimnames(aic)[[2]]<-paste("MA-order",0:maxmaorder,sep=" ")
    dimnames(aic)[[1]]<-paste("AR-order",0:maxarorder,sep=" ")
    bic<-aic
    for (j in 0:maxarorder)
    {
      for (k in 0:maxmaorder)
      {
        aic[j,k]<-log_sigma_jk[j,k]+2*(j+k)/len
        bic[j,k]<-log_sigma_jk[j,k]+log(len)*(j+k)/len
      }
    }
    aic
    bic
    which(aic == min(aic), arr.ind = TRUE) 
    which(bic == min(bic), arr.ind = TRUE) 
    AICarorder<-which(bic == min(aic), arr.ind = TRUE)[1]-1 
    AICmaorder<-which(bic == min(aic), arr.ind = TRUE)[2]-1 
    BICarorder<-which(bic == min(bic), arr.ind = TRUE)[1]-1
    BICmaorder<-which(bic == min(bic), arr.ind = TRUE)[2]-1
    results <- matrix(c(AICarorder,AICmaorder,BICarorder,BICmaorder), nrow=2)
    colnames(results) <- c("AIC", "BIC")
    rownames(results) <- c("AR-Order", "MA-Order")
    results
    
    
    
    
    
    
    
    
    # Plotten des Predictes
    AICBIC <- "BIC"  # AIC oder BIC eingeben
    schritte <- 20 # Anzahl Schritte in die Zukunft
    
    arorder <- results[1,AICBIC]
    maorder <- results[2,AICBIC]
    y_obj<-arima(y,order=c(arorder,d,maorder))
    
    # Ist der tsdiag() in Ordnung?!?
    tsdiag(y_obj)
    
    par(mfrow=c(2,1))
    ts.plot(y)
    y_pred<-predict(y_obj,n.ahead=schritte)
    
    ymin<-min(c(as.vector(y),y_pred$pred-2*y_pred$se))
    ymax<-max(c(as.vector(y),y_pred$pred+2*y_pred$se))
    
    ts.plot(c(as.vector(y),y_pred$pred),ylim=c(ymin,ymax))
    abline(v=length(y))
    pre <- y_pred$pred
    lines(length(y):(length(y)+length(pre)),c(y[length(y)],pre+2*y_pred$se),col="2")
    lines(length(y):(length(y)+length(pre)),c(y[length(y)],pre-2*y_pred$se),col="2")
  }
  
  ### MSFE
  # Out of sample -> modell erstellen, predicten und dann vergleichen
  # in sample -> modell erstellen, und mit wissen von bekannter "Zukunft" anpassen.
  function(){
    ## out-of-sample
    function(){
      len<-300              # l??nge der Datens??tze
      simanz<-100           # wie viele simulazionen
      h<-10                 # anz forecast werte
      set.seed(10)
      ARorderBlau<-8        # gesch??tze anz parameter 1
      ARorderRot<-6         # gesch??tze anz parameter 2
      a_vec<-c(0.1,0.1,0.1) # unser Orginalprozess
      #mu<-6/(1-sum(a_vec))
      
      
      # Function computes out-of-sample squared forecast errors
      for_sim_out<-function(h,len,j,ARorderBlau,ARorderRot,a_vec)
      {
        mu<-6/(1-sum(a_vec))
        set.seed(j)
        x<-mu+arima.sim(n=len+h,list(ar=a_vec))
        y_ar3<-arima(x[1:len],order=c(ARorderBlau,0,0))
        y_ar9<-arima(x[1:len],order=c(ARorderRot,0,0))
        forecast_error_true_out<-((mean(x[1:len])-x[(len+1):(len+h)])^2) 
        forecast_error_ar3_out<-(predict(y_ar3,n.ahead=h)$pred-x[(len+1):(len+h)])^2 
        forecast_error_ar9_out<-(predict(y_ar9,n.ahead=h)$pred-x[(len+1):(len+h)])^2
        forecast_error<-c(forecast_error_true_out,forecast_error_ar3_out,forecast_error_ar9_out)
        return(forecast_error)
      }
      
      # Generate 1000 realizations and compute squared forecast errors
      sample_id_out <- foreach(j = 1:simanz,.combine=rbind) %dopar% for_sim_out(h,len,j,ARorderBlau,ARorderRot,a_vec)
      # Compute mean-square forecast errors for h=1,...,10
      out_of_sample<-cbind(apply(sample_id_out[,1:h],2,mean),apply(sample_id_out[,h+1:h],2,mean),apply(sample_id_out[,2*h+1:h],2,mean))
      
      ymin<-min(apply(out_of_sample,1,min))
      ymax<-max(apply(out_of_sample,1,max))
      ts.plot(out_of_sample[,1],xlab="Forecast horizon",ylab="MSFE",main="Out-of-sample",ylim=c(ymin,ymax))
      lines(out_of_sample[,2],col="blue")
      lines(out_of_sample[,3],col="red")
      mtext("True (noise model)", side = 3, line = -1,at=h/2,col="black")
      mtext(c("AR order                     ",ARorderBlau), side = 3, line = -2,at=h/2,col="blue")
      mtext(c("AR order                     ",ARorderRot), side = 3, line = -3,at=h/2,col="red")
    }
    
    ## in-sample
    function(){
      set.seed(10)
      len<-100
      a_vec<-c(0.2,0.2)
      mu<-6/(1-sum(a_vec))
      simanz<-100
      h<-10
      or<-6
      ob<-2
      # Function computes in-sample squared forecast errors
      for_sim_in_ar3<-function(h,len,a_vec,mu,j,or,ob)
      {
        set.seed(j)
        x<-mu+arima.sim(n=len+h,list(ar=a_vec))
        # AR(2)
        p<-ob
        y_ar<-arima(x,order=c(p,0,0))
        x_for<-x[(len-p+1):(len)]
        a_for<-y_ar$coef[1:p]
        mu_hat<-y_ar$coef[p+1]
        for (i in 1:h)
        {
          x_for<-c(x_for,a_for%*%(x_for[length(x_for):(length(x_for)-p+1)]-mu_hat)+mu_hat)
        }
        x_for_h<-x_for[p+1:h]
        forecast_error_true_in<-(x_for_h-x[(len+1):(len+h)])^2
        # AR(3)
        p<-3
        y_ar<-arima(x,order=c(p,0,0))
        x_for<-x[(len-p+1):(len)]
        a_for<-y_ar$coef[1:p]
        mu_hat<-y_ar$coef[p+1]
        for (i in 1:h)
        {
          x_for<-c(x_for,a_for%*%(x_for[length(x_for):(length(x_for)-p+1)]-mu_hat)+mu_hat)
        }
        x_for_h<-x_for[p+1:h]
        forecast_error_ar3_in<-(x_for_h-x[(len+1):(len+h)])^2
        # AR(9)
        p<-or
        y_ar<-arima(x,order=c(p,0,0))
        x_for<-x[(len-p+1):(len)]
        a_for<-y_ar$coef[1:p]
        mu_hat<-y_ar$coef[p+1]
        for (i in 1:h)
        {
          x_for<-c(x_for,a_for%*%(x_for[length(x_for):(length(x_for)-p+1)]-mu_hat)+mu_hat)
        }
        x_for_h<-x_for[p+1:h]
        forecast_error_ar9_in<-(x_for_h-x[(len+1):(len+h)])^2
        forecast_error<-c(forecast_error_true_in,forecast_error_ar3_in,forecast_error_ar9_in)
      }
      # Generate 1000 realizations and compute squared forecast errors
      sample_id <- foreach(j = 1:simanz,.combine=rbind) %dopar% for_sim_in_ar3(h,len,a_vec,mu,j,or,ob)
      # Compute mean-square forecast errors for h=1,...,10
      in_sample<-cbind(apply(sample_id[,1:h],2,mean),apply(sample_id[,h+1:h],2,mean),
                       apply(sample_id[,2*h+1:h],2,mean))
      ymin<-min(apply(in_sample,1,min))
      ymax<-max(apply(in_sample,1,max))
      ts.plot(in_sample[,1],xlab="Forecast horizon",ylab="MSFE",
              main="True model AR(3): in-sample MSFE",ylim=c(ymin,ymax))
      lines(in_sample[,2],col="blue")
      lines(in_sample[,3],col="red")
      mtext(c("Richtiger AR                                ",length(a_vec)), side = 3, line = -1,at=h/2,col="black")
      mtext(c("AR ordnung                            ",ob), side = 3, line = -2,at=h/2,col="blue")
      mtext(c("AR ordnung                            ",or), side = 3, line = -3,at=h/2,col="red")
    }
    
  }
  
  
  
  # Saisonalit??t kontrollieren
  
  #Beispiel Daten
  function(){
    xdat<-read.table(file = "Users/davidkuchelmeister/Studium/2. Jahr/4. Semester/OESSY - Server HSNG WMAR/WILDI/Unterlagen_Daten_OESSY2//first_long.txt"
                     ,header = T)
    x<-as.ts(na.exclude(xdat))
  }
  
  per<-function(x,plot_T){
    len<-length(x)
    per<-0:(len/2)
    DFT<-per
    
    for (k in 0:(len/2))
    {
      cexp <- complex(arg=-(1:len)*2*pi*k/len)
      DFT[k+1]<-sum(cexp*x*sqrt(1/(2*pi*len)))
    }
    per<-abs(DFT)^2
    if (plot_T){
      par(mfrow=c(2,1))
      plot(per,type="l",axes=FALSE,xlab="Frequency",ylab="Periodogram",
           main="Periodogram")
      axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
                                      "4pi/6","5pi/6","pi"))
      axis(2)
      box()
      
      plot(c(NA,log(per)[2:length(per)]),type="l",axes=FALSE,xlab="Frequency",ylab="Log-periodogram",
           main="Log-periodogram")
      axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
                                      "4pi/6","5pi/6","pi"))
      axis(2)
      box()
    }
    return(list(DFT=DFT,per=per))
  }
  # falls es fehler gibt axes="F" probieren
  
  # falls es ein Peek bei pi/6 oder einem vielfachen davon gibt ist es eine Saisonale Reihe
  # pi/6 = Widerholung 1mal pro Jahr   => pi = Widerholung 6mal pro Jahr
  per(x-mean(x),T)
  
}
# STOP
function(){
  # Matrix quadrieren
  install.packages("expm")
  library(expm)
  
  # P = Matrix ; T = Anzahl Simulationen bis Zeitpunkt (t) ; x0 = Startpunkt
  markovN <- function(P, T, x0){
    vektor <- rep(NA,T)
    vektor[1] <- x0
    for(i in 2:T){
      vektor[i] <- sample(1:dim(P)[1],size = 1,prob=c({P[vektor[i-1],1:dim(P)[1]]}))
    }
    vektor
  }
  
  # Zeichnen der Markov Kette
  library(markovchain)
  kette <- new("markovchain", transitionMatrix=P)  # states= anzahl Dimensionen in Matrix P
  summary(kette)
  plot(kette)
  
  initialState <- c(0, 1, 0)
  after7Days <- initialState * (mcWeather ^ 7)
  
  transitionProbability(mcWeather, "cloudy", "rain")
  
  # mittlere Dauer ausrechnen
  {
  (P4 = P)
  Dauerbis <- 3 # Mittlere Dauer bis man im Zustand 3 ist.
  P4[, Dauerbis] <- 0  #NICHT ??NDERN
  solve(diag(dim(P4)[1]) - P4) %*% rep(1, dim(P4)[1])
  # Startwert auslesen als Position im Vektor!
  }
  
  
  ??t <- function(P, ??0 ,t){
    ??t <- ??0%*%P^t
    ??t
  } # Wahrscheinlichkeiten Zeitpunkt t
  Kt <- function(??t, c){
    Kt <- ??t*t(c)
    Kt
  } # Kosten Zeitpunkt t
  sum(Kt_erg) # total erwartete Kosten
  K <- function(??0, Pt, c){
    K <- ??0 * sum(Pt) * t(c)
    K
  }
  
  
  # Konvergenz der Wahrscheinlichkeiten mit EIGENVEKTOR
  eigen(P) # an der Stelle wo eigenWERT = 1 => Vektor auslesen
  eigV <- eigen(t(P))$vectors[,1]
  eigV/sum(eigV) # Normieren des Vektors
  
  
  # Weibul-Verteilung
  plot(dweibull(x = 1:1000,shape = 1,scale = 100),type="l") # Beispiel
  plot(dweibull(x = , shape = Beta, scale = 1/Lambda))
  pweibull()
  qweibull()
  rweibull()
  
  # Berechnung des Mittelwertes
  1/lambda * gamma(1 + 1/beta)
  
  
  
  
  
  ### Stetiger Markov-Prozess
  ### Zeitkontinuierlich
  
  # Simulieren und Plotten von Zust??nden (Daten k??nnen nicht ausgelesen werden)
  {
  P <- matrix(c(0, 0.2, 0.8, 0.1, 0.0, 0.9, 0.5, 0.5, 0.0), ncol=3, byrow=TRUE)
  lambda <- c(0.5, 0.05, 0.1);
  # Anfangszustand nicht vergessen!
  
  par(mfrow=c(4,4), mar=c(3,2,0,0), oma=c(0,0,0,0))
  for (run in 1:16) {
    plot(NULL, xlim=c(0, 20), ylim=c(0, 3), ylab="Zustand", xlab="Zeit", main=NULL)
    z <- 0;
    state <- 1; #Anfangszustand
    repeat{
      oldState <- state
      state <- 1
      sum
      #Siehe auch Aufgabe zur Simulation einer Markov-Kette
      cum <- 0
      num <- runif(1)
      for (i in 1:dim(P)[1]) {
        cum <- cum + P[oldState,i]
        if (num <= cum) {
          state <- i;
          break
        }
      }
      delta <- rexp(1, rate=lambda[state]);
      lines(matrix(c(z, oldState, z + delta, oldState, z + delta, state), ncol=2, byrow=T))
      z <- z + delta;
      if (z > 20) {
        break
      }
    }
  }
  }
  
  
  P <- matrix(c(0, 0.2, 0.8, 0.1, 0.0, 0.9, 0.5, 0.5, 0.0), ncol=3, byrow=TRUE) # Sprungwahrscheinlichkeiten
  lambda <- c(0.5, 0.05, 0.1) # Wartezeitparameter
  runs <- 10000
  Zeit <- 10.2    # Zeitpunkt wann ausgewertet wird
  Anfangszustand <- 1
  
  # Falls nicht funktioniert, Zeilen einzeln ausf??hren
  SimZust <- function(P,lmabda,Zeit,Anfangszustand,runs){
    states_at_102 <- c(0, 0, 0)
    for (run in 1:runs) {
      z <- 0;
      state <- Anfangszustand; #Anfangszustand
      repeat{
        oldState <- state
        delta <- rexp(1, rate=lambda[state]);
        cum <- 0;
        num <- runif(1)
        for (i in 1:dim(P)[1]) {
          cum <- cum + P[state,i]
          if (num <= cum) {
            state <- i;
            break
          }
        }
        z <- z + delta;
        if (z > Zeit) {
          states_at_102[oldState] <- states_at_102[oldState] + 1;
          break
        }
      }
    }
    states_at_102 / sum(states_at_102)
  }
  SimZust(P,lmabda,Zeit,Anfangszustand,runs)
  
  
  # Generator-Matrix berechnen Q
  {
  # Ben??tigte Ausgangsinfomrationen
  library(expm)
  P <- matrix(c(0, 0.2, 0.8, 0.1, 0.0, 0.9, 0.5, 0.5, 0.0), ncol=3, byrow=TRUE)
  lambda <- c(0.5, 0.05, 0.1);
  
  # Berechnen der Generatoren Matrix
  P <- matrix(c(0,0.3,0.7, 0.2,0, 0.8, 0.1,0.9,0), byrow=TRUE, ncol=3)
  R <- P / lambda
  Q <- R - diag(rowSums(R))
  
  # Zeitpunkt ausrechnen
  expm(Q * 10.2)
  c(1,0,0) %*% expm(Q * 10.2) # mit Startvektor
  
  }
  
  # Ausfallraten bestimmen:
  lambda_M = 1 / (Mittlere Laufzeit bis zum Ausfall)
  # lambda_M nicht das selbe Lambda wie der Wartezeitparameter
  R <- matrix(c(0,0.006,0.1,0), nrow=2, byrow=T)
  
  # mittlere Aufenthaltszeiten in den verschiedenen Zust??nden
  Q # Q-Matrix => Betrag der Diagonalen    # falls R gegeben => Q berechnen
  
  # Wahrsch von Zustand i zu Zustand j
  R_ij = lambda_i * P_ij
  P_ij = R_ij / lambda_i
  
  # asymptotische L??sung / Konvergenz der Zeitreihen
  # Asymptotische Entwicklung
  es = eigen(t(Q))
  # Zeile mit dem Eigenwert nahe 0 wird gew??hlt f??r k
  vv = es$vectors[,k]
  (ev = vv/sum(vv)) # Normieren!
  
  # Anzahl ??berg??nge pro Zeiteinheit
  R
  # Anzahl ??berg??nge pro Zeiteinheit auf lange Sicht
  R[von,nach] * ev[]
  
}
# STMO
function(){
  # Formelsammlung Regression in Befehle3 anschauen!!!
  library(leaps); library(car)
  
  # Modell exakt beschreiben
  #   Betas (??) und sigma mit HUT!
  #   wenn sie aus dem summary(lm.erg) gelesen werden
  
  #   y = ?? + ??x + E (mit E unabh??nhig N(0,sigma) verteilt
  
  # Hypothesentest
  #   H0 ist ??_n = 0 (n = 1-n)    Alle Betas ausser Beta0!
  #   ist das P-Value kleiner als das Signifaktsniveau => H1 wird angenommen / H0 verworfen
  
  #   Fehler 1. Art => H1 angenommen und H0 w??hr korrekt
  #   Eintrittswahrscheinlichkeit = z.B. 2% = Hypothesentest Niveau
  #   Fehler 2. Art => H0 behalten obwohl H1 korrekt w??hr
  #   Eintrittswahrscheinlichkeit nicht berechenbar
  
  # Tukeys first Aid
  function(){
    1. Logarythmus
    falls Konzentrationen oder Betr??ge [ log(y) ]
    2. Wurzelziehen
    falls Z??hlvariablen [ sqrt(y) ]
    3. Arcus-Sinus-Tranformation
    falls y ein Anteil   [ log((y+0.005)/(1.01-y)) ]
  }
  
  # Ausschliessverfahren der unn??tigen Parameter
  function(){ 
    head(AKW <- read.table("Studium/2. Jahr/4. Semester/STMO - Server rkst/StMo/Daten/AKW.dat",header = T))
    AKW$lgK <- log(AKW$K); AKW$sqrtN <- sqrt(AKW$N); AKW$lgG <- AKW$G; head(AKW)
    CEDHEC <- read.table(paste("Studium/2. Jahr/4. Semester/STMO - Server rkst/StMo/Daten/CastleEDHEC.dat",sep=""), header=T)
    
    # F-Wert verfahren
    function(){ 
      # R??ckwertsselektioniert
      AKW.lmV <- lm(lgK ~ lgG + D + WZ + BZ + Z + NE + KT + BW + sqrtN + KG, data=AKW)
      drop1(AKW.lmV, test="F")
      AKW.lmVnxt <- update(AKW.lmV, . ~ . - WZ )
      AKW.lmVnxt <- update(AKW.lmV, . ~ . - WZ - BW - Z - BZ - sqrtN - KT)
      drop1(AKW.lmVnxt, test="F")
      
      step(AKW.lmV, direction="backward")
      
      # Vorw??rtsselektioniert
      AKW.lmA0 <- lm(lgK ~ 1, data=AKW)
      AKW.lmA0 <- update(AKW.lmA0, . ~ . + WZ)
      add1(AKW.lmA0, test="F")
      
      step(AKW.lmA0, direction="forward", scope=list(upper = ~ lgG + D + WZ + BZ + Z + NE + KT + BW + sqrtN + KG, lower = ~ 1))
      
      # Beidseitig
      lm.erg <- lm(lgK ~ KG, data=AKW)
      step(AKW.lmA0, scope=list(upper = ~ lgG+ D + WZ + BZ + Z + NE + KT + BW + sqrtN + KG, lower = ~ 1), direction="both")
      update(AKW.lmV, . ~ . - BW) # - BW = abzeihen     + BW = zuz??hlen
      update(AKW.lmV, . ~ . - WZ - BW - Z - BZ - sqrtN - KT)
      
      # Analysieren der P-Werte (entweder mit summary oder...)
      add1(AKW.lmV, test="F")
      drop1(AKW.lmV, test="F")
    }
    
    # AIC-Verfahren (BIC)
    function(){
      # Modelle der 3 Varianten sind teils unterschiedlich, m??ssen aber alle ber??cksichigt werden
      
      # Ausschliessverfahren wird Schritt f??r Schritt automatisch gel??st
      x1 <- step(AKW.lmV, direction="backward")
      x2 <- step(AKW.lmA0, direction="forward",
                 scope=list(upper = ~ lgG + D + WZ + BZ + Z + NE + KT + BW + sqrtN + KG, lower = ~ 1))
      x3 <- step(AKW.lmA0, direction="both" ,
                 scope=list(upper = ~ lgG+ D + WZ + BZ + Z + NE + KT + BW + sqrtN + KG, lower = ~ 1))
      # Falls der AIC-Wert von <none> nur wenig gr??sser ist als der AIC der tiefsten erkl. Var gibt es nur ein "kleines" Minimum
      
      # BIC eine andere Berechnung f??r den step
      mfm.lm1 <- lm(FoHF ~ ., data=CEDHEC); require(leaps)
      x4 <- step(mfm.lm1, k=log(nrow(CEDHEC))) # k = Anzahl Beobachtungen der Datenmatrix
      
      # Zusamenfassung, was gemacht wurde
      x2$anova
    }
    
    # Cp-Verfahren
    function(){
      library(leaps)
      library(car)
      
      # nbest = bis 6 erkl. Var. runter berechnen,   nvmax = bis zur 10 erkl. Var. berechenn (falls so viele vorhanden)
      AKW.Cp <- regsubsets(lgK ~ lgG + D + WZ + BZ + Z + NE + KT + BW + sqrtN + KG, nbest=6, nvmax=10, data=AKW)
      summary(AKW.Cp)
      plot(AKW.Cp, scale="Cp")
      
      h <- subsets(AKW.Cp, statistic="cp", legend=TRUE, min.size=4, main="Mallow Cp", cex.subsets=0.7, las=1)
      abline(a=1, b=1, lty=2)
      
      # Modelle    => Das kleinste Modell soll beachtet werden
      #               Das Modell mit den wenigsten erkl. Var. (am weitesten links) ist jedoch das wichtigste, wenn es funktioniert
      #               Modelle in der N??he der Linie k??nnen auch gew??hlt werden, falls die restlichen Modelle zu schlecht sind
      
      # Die Abline => ist der Erwartungswert des Modells
      #               Wir erwarten dort den Wert, wenn das Modell alle n??tigen Variablen verwendet.
      
      
    }
    
    # Zusamenfassung und Vergleich
    function(){   
      Beruhend auf F-Wert, schrittweise r??ckw??rts:
        Resultat: lm(lgK ~ lgG + D + NE + KG)
      
      Beruhend auf F-Wert, schrittweise vorw??rts:
        Resultat: lm(lgK ~ lgG + D + NE + KG )
      
      Beruhend auf AIC, schrittweise r??ckw??rts:
        Resultat: lm(lgK ~ lgG + D + BZ + Z + NE + KT + sqrtN + KG)
      
      Beruhend auf AIC, schrittweise vorw??rts:
        Resultat: lm(lgK ~ KG + lgG + D + NE + KT + sqrtN)
      
      Beruhend auf AIC, "both":
        Resultat: lm(lgK ~ KG + lgG + D + NE + KT + sqrtN)
      
      Beruhend auf Cp:
        lm(lgK ~ lgG + D + NE + KT + KG)
      oder
      lm(lgK ~ lgG + D + NE + KT + sqrtN + KG)  ## hat minimales Cp
      
      Wie die Zusammenfassung zeigt, kommt man mit den verschiedenen Verfahren zu
      zwei unterschiedlichen Modellen.
      
      Potentiell "gute" Modelle sind:
        i)  lm(lgK ~ lgG + D + NE + KT + sqrtN + KG)
und
ii) lm(lgK ~ lgG + D + NE + KG)

Das Modell i) ist das "beste" Modell aus dem Verfahren mit dem Cp-Kriterium,
wie auch mit dem AIC (direction="both").
Das Modell ii) ist das kleinste Modell und das "beste" Modell, wenn der
F-Wert als Kriterium mit Stoppwert 5% genommen wird. Die Koeffizienten sind
dann alle "signifikant verschieden von Null" (Diese Aussage ist aber
                                              unzuverl?ssig!!). Es ist auch gem?ss dem Cp ein valides Modell, da es
immer noch einen Cp-Wert um die Winkelhalbierende und am wenigsten
Parameter hat.
Wir k?nnen Modell ii) leichte Priorit?t geben?ber Modell i) geben. Welches
soll nun gew?hlt werden? - Auf jeden Fall sollten die Resultate einer
Residuenanlyse in die Entscheidung einbezogen werden. Vielleicht gibt es
auch noch fachliche Argumente f?r order gegen eine Variante.

# Mit dem Gefundenen Modell muss auch eine Residuen Analyse durchgef??hrt werden
# und beschrieben werden, ob alle Messpunkte in Ordnung sind.
# plot(lm.erg), plot.lmSim(lm.erg) und plot(r.lm.erg)
    }
    
    # VIF (variance inflation factor, Qualit??t nachkontrollieren)
    function(){
      # Wie korreliert die Varianz der einzelnen Variablen?
      # korrelieren sie zu stark miteinander, schaukelt sich Ihre Varianz auf, was zu st??rkeren Schwankungen f??hrt.
      
      library(car)
      vif(lm.erg)
      #      Temp.Tank        Temp.Gas      Vapor.Tank Vapor.Dispensed
      #      12.997379        4.720998       71.301491       61.932647
      
      # Alle bis auf Temp.Gas haben einen VIF gr?sser als 10. Gem?ss Faustregel
      # sind also Probleme mit Kollinearit?ten vorhanden. Die Variable 'Vapor.Tank'
      # ist gem?ss VIF am meisten davon betroffen.
      
      # kann verbessert werden indem man die Daten berabeitet:
      # 1.
      # Residuen Analyse der Daten betrachten, m??gliche Ausreisser ausschliessen
      #
      # 2.
      # par(daten), ??ber die am st??rksten korrellierten Daten
      # werden der mean und die Differenz genommen und als neue Datenzeilen hinzugef??gt
      #   mean <- mean(dat[1:100,1],dat[1:100,2])
      #   diff <- dat[,1]-dat[,2]
      #
      # 3.
      # eine neues lm()-Modell wird mit den signifikanten daten von vorhin und
      # den neu erstellten Daten gemacht werden.
      #
      # 4.
      # danach wird noch step() ??ber das neue Modell angewendet um die unn??tigen
      # erkl??renden Variablen auzuschliessen
    }
  }
  
  # Splines
  function(){
    library(MASS) # Daten f??r Beispiel
    setwd("~")
    volDrop <- read.table("Studium/2. Jahr/4. Semester/STMO - Server rkst/StMo/Daten/VoltageDrop.dat", sep="", header=TRUE)
    vierC <- read.table("Studium/2. Jahr/4. Semester/STMO - Server rkst/StMo/Daten/vierC.dat", sep="", header=TRUE)
    
    # Polynomiale Regression
    function(){
      plot(mcycle)
      for(i in c(3,6,12)){
        lm.erg <- lm(accel~poly(times,i),mcycle)
        x <- data.frame(times = seq(min(mcycle$times),max(mcycle$times)))
        a <- predict(lm.erg,x)
        lines(a,col=2)
      }
    }
    
    # Nat??rliche Splines
    function(){
      plot(mcycle)
      x <- times
      for(i in c(5,10,20)){
        library(splines)
        lm.erg2 <- lm(accel~ns(times, df=i), data=mcycle)
        b <- predict(lm.erg2,x)
        lines(b,col=2)
      }
    }
    
    
    ## ACHTUNG: Default Einstellungen bei loess unterscheiden sich von
    ##          jenen bei scatter.smooth.
    ## Defaults bei loess: span=0.75, degree=2, family="gaussian"
    ## Defaults bei scatter.smooth: span=2/3, degree=1, family="symmetric"
    
    # Smoothing Splines
    function(){ 
      #   Smoothing Parameter  spar= 0.5261326  lambda= 1.197378e-05 (12 iterations)    => spar = 0-1
      #   Equivalent Degrees of Freedom (Df): 20.00239                                  => Optimale anzahl Freiheitsgrade
      #   Penalized Criterion: 34687.98                                                 => 
      #   GCV: 604.8662                                                                 => 
      plot(mcycle)
      for(i in 8:12) lines(smooth.spline(x=mcycle$times,y=mcycle$accel,df = 12.20876,spar = 0.6598558),col=2)
    }
    
    # LOWESS-Verfahren Lokale Regression
    function(){
      plot(VolDrop ~ time, data=volDrop)
      
      vD.lr1s <- loess(VolDrop ~ time, data=volDrop, span=0.1, degree=1, family="symmetric")
      hx <- seq(min(volDrop$time), max(volDrop$time), length=100)
      lines(hx, predict(vD.lr1s, newdata=data.frame(time=hx)), col="darkgreen")
      
      scatter.smooth(volDrop$time, volDrop$VolDrop, span=0.3, degree=2, family="symmetric")
      
      
      # Falls Ausreisser heruntergewichtet werden sollen: (Unter der Annahme das Ausreisser vorhanden sind)
      scatter.smooth(vierC$lCarat, vierC$lPrice, span=1, degree=1, family="symmetric")
    }
    
    # GAM (Regression mit Splines) additive Splines
    function(){
      library(gam)
      gam.erg <- gam(Price ~ lo(Carat, span=0.25, degree = 1), data=vierC) # Optimales lCarat heraus finden mit LOWESS
      par(mfrow=c(1,1)) # anzahl Plots = anzahl Variablen
      plot(gam.erg, resid=TRUE, se=TRUE)
      
      # das gam Modell mit dem normalen lm Modell vergelichen
      source("Studium/2. Jahr/4. Semester/STMO - Server rkst/StMo/RSourceCode/RFn_Plot-lmSim.R")
      vierC.lm <- lm(Price ~ Carat + Colour + Clarity + CBody , data=vierC)
      vierC.gam <- gam(Price ~ lo(Carat, span=0.25) + Colour + Clarity + CBody, data=vierC)
      par(mfrow=c(2,4))
      plot(vierC.lm)
      stats:::plot.lm(vierC.gam)
      
      par(mfrow=c(2,2))
      plot(vierC.gam, resid=TRUE, se=TRUE)
      summary(vierC.gam)
      
      # 2. m??glichkeit VORSICHT!!!
      function(){
        library(mgcv)
        vierC.gam <- gam(lPrice ~ s(lCarat) + Colour + Clarity + CBody , data=vierC)
        par(mfrow=c(1,1))
        plot(vierC.gam, residuals=TRUE, se=TRUE)
        summary(vierC.gam)
      }
    }
    
    
  }
  
  # Regressionsalayse
  function(){
    # Leverigen => Geld aufnehmen um damit zu investieren
    
    lm(Y-Achse ~ . , data=Datensatz) # der Punkt bedeuted, dass alle Daten zur erkl??rung der Zielvariable benuzt werden
    lm(Y-Achse ~ 1-X-Achse + 2-X-Achse + Stelle, data=Datensatz)
    
    lm(Y-Achse ~ X-Achse, data=Datensatz, weights = Datensatz[,i])
    # gewichtete Regression (plot Analyse nicht mehr m??glich)
    # die gewichtete Regression darf nicht logarythmiert werden
    # ergibt die lineare Gleichung der Linie=> y = b0 + b1*x + b2*x_2 + ...
    
    summary(lm.erg)
    coef(lm.erg) # genaue Werte zu ??0,??1,...,??N  (VORSICHT wegen sigma)
    
    # p-value < 0.001     =>  suggest a highly significant slope
    # P(>|t|)             =>  H0 : ??1=0      HA : ??1 ??? 0
    # F-Statistic         =>  Unterscheiden sich die zus??tzlichen Messungen signifikan vom intercept? F<0.05 JA
    #                         gibt es mind eine Erkl??rende Variable die, die Zielvariable erkl??ren k??nnte?
    # R-Squared           =>  1 ist optimal 0 ist extrem schlecht
    #                         z.B 0.8 =>  80% von der variabilit??t (die Streuung) der Zielvariable wird erklr??rt, 
    #                                     20% rest bleibt ??brung, variabilit??t steckt in den residuen.
    
    # Vergleichen von zwei Regressionen
    anova(lm.erg1, lm.erg2)
    
    # Werte voraussagen (Predict)
    function(){
      # im Data-Frame m??ssen die Koordinaten (oder Bereich) der Erkl??renden Variablen
      #   drin sein um an diesem Wert eine Voraussage machen zu k??nnen
      
      newdat = data.frame(Ort.der.1.erkl.Var = 100.5, Ort.der.2.erkl.Var = 4.5) # Werte sollten den ganzen X-Achsenbereich abdecken. (z.B. mit seq())
      predict(lm.erg, newdata = newdat, interval="prediction", level=0.95)
      predict(lm.erg, newdata = data.frame(1:100), interval="prediction", level=0.95)
      
      predict(lm.erg, newdata=data.frame(x-Achsen-benennung=50))
      predict(lm.erg, newdata=data, interval="confidence", level=0.95)
    }
    
    # Beschreiben von Plots
    #   Welcher Plot?
    #   Was wird gepr??ft?
    #   Was f??llt auf?
    #   Wie ist das Urteil?
    
    # Residuen- und Sensitivit??ts-Analyse mit Plots
    function(){
      par(mfrow=c(2,4))
      plot(lm.erg)
      # plot(lm.erg ,which=1:6)
      # Bedeutung der Plots:
      #   Plot 1: Tukey-Anscombe-Plott, Muss möglichst waagrecht sein (Bananenform, Badewannenform, usw. = schlecht)
      #           Dieser Plot beschreibt den ob der Erwartungswert konstant ist.
      #           => Ist der Erwartungswert konstant?
      #               plot(lm.erg$residuals)
      #   Plot 2: QQPlot = Sind die Residuen Normal verteilt? 
      #           Langschänzig Ausreisser:(links-unten,rechts-oben); kurzschw??nzigkeit (Gegenteil)
      #           Linksschiefschief Ausreisser:(links-unten,rechts-unten); Rechtsschief (Form von x^2)
      #           Liegen die Punkte innerhalb der stochastischen Fluktuation?
      #           => Gibt es Evidenz gegen eine Annahme der Normalverteilung
      #               qqnorm(lm.erg$residuals)
      #   Plot 3: Scale-Location Plot ist die Gleichverteiltheit des Tukey-Anscombe-Plotts.
      #           Ist die Fehlervariabilit?t konstant?
      #           Zeigt der Gl??tternen einen Anstieg oder Ausschlag der Fehlervariabilit??t?
      #           liegt dieser Ausschlag noch innerhlab der stochastischen Fluktuationen?
      #           => Also keine Evidenz gegen konstante Fehlervariabilität.
      #   Plot 4: Gibt es zu einflussreiche Beobachtungen? (bösartig, Hebelpunkt, gutmütig)
      #           Cooks distance, calculates how big is the influence on the regression if this measurement would be removed.
      #           High value influneces the lm a lot
      #   Plot 5: Gibt es zu einflussreiche Beobachtungen? (bösartig, Hebelpunkt, gutmütig)
      #           horizontal = Hebelwirkung , Vertikal = Einflussstärke
      #           x-Achse ist die Hebelwirkung und y-Achse ist der Einfluss
      #           Faustregel 1 = ab Cooksdistance 1 ist der Einfluss zu gross
      #           Faustregel 2 = x-Achse über 0.2 = Hebelpunkt! aber nicht unbedingt zu einflussreich
      #   Fazit:  gibt es verstösse gegen die Annahmen? Sind wir zufireden mit den Plots?
      
      source("/Users/davidkuchelmeister/Documents/Studium/2. Jahr/4. Semester/STMO - Server rkst/StMo/RSourceCode/RFn_Plot-lmSim.R")
      plot.lmSim(lm.erg,which=1:3)
      # Die Roten linien m??ssen die Grauen ??berdecken und m??ssen in die gleiche Richutung zeigen.
      #   (Sreuung ist konstant oder nicht)
      
      # Bedeutung der Plots:
      #   Plot 5: Tukey-Anscombe-Diagramm Simuliert
      #           => f??llt der erwartete Erwartungswert aus den Simulationen
      #   Plot 6: 
      #   Plot 7: Scale-lokation-Plot
      #           => steigt oder f??llt die Varianz aus den Simulierten Werten?
      #   Fazit:  Ist die Anpassung befriedigend? wegen was ja oder nein?
      
      # korrekte Beschreibung
      # Wichtige W??rter: (b??sartig, Hebelpunkt, langschw??nzigkeit, Rechtsschief, keine Auff??lligkeiten)
      {
        # Kurzfassung:
        # TA-Diagramm: Erwartungswert i.O.; allerdings zwei Ausreisser sichtbar: 1, 20
        # Normal-Plot: Zwei Ausreisser 1 und 20 sichtbar (wegen dem Verlauf der
        #              zentralen Punkte nach OBEN konnte man auch als rechts langschw??nzig
        #              bezeichnen), sonst i.O.
        # Scale-Location: i.O. (obwohl die rote Kurve in der Simulation z.T. am Rande liegt
        # Res. vs Leverage: Beob. 1 ist ein b?sartiger Hebelpunkt (Cook's Distanz > 1)
        #                   Zwei weitere Beobachtungen sind Helbelpunkte (Hii>> 0.2),
        #                   jedoch gutartig
      }
      
      # Anlayse plots:
      plot(Ox$erklärender Variablen, resid(lm.erg.opt)); abline(h=0, lty=3)
      
      
      ############# Check autocorrelation
      acf(residuals(lm.erg))
      # lines are over the doted line?
      # if yes, use the 
    }
    
    # Generalized Least Squares Regression
    function(){
      library(nlme)
      mdl.ac <- gls(birds ~year, data=dat, 
                    correlation = corAR1(form=~year),
                    na.action=na.omit)
      summary(mdl.ac)
      
      # AR(1) errors within each Mare
      fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
                 correlation = corAR1(form = ~ 1 | Mare))
      # variance increases as a power of the absolute fitted values
      fm2 <- update(fm1, weights = varPower())
    }
    
    # Robuste Regression
    function(){
      # Messungen werden gewichtet, stark ausschlagende werden weniger gewichtet als weniger ausschlagende
      library(robustbase)
      r.lm.erg <- lmrob(y ~ x1 + x2, data=D.synt, setting="KS2011")
      summary(r.lm.erg)
      
      plot(r.lm.erg) # ergibt 5 Grafiken
      # Nr. 1 = Daten m??ssen in Quadrat drin liegen (Es gibt keine Hebelwirkungen)
      # Nr. 2 = muss um Kurve streuen
      # Nr. 3 = Daten sollen um Winkelhalbierende streuen
      # Nr. 4 = 
      # Nr. 5 = Scale-Location Plot (Wurzel aus Daten)
    }
    
    # confidenz Intervall
    confint(lm.erg, level=0.95)
    
    # gefittete Werte von y
    lm.erg$fitted.values
    lm.erg$coefficients[2]
    
    #linie in Scatterplot einzeichnen
    abline(a=lm.erg[1],b=lm.erg[2])
    
    # Faktor Variablen erstellen
    #   Faktorvariable (x <- as.factor(x))
    #     Regressionsgerade = Hyperebene
    #     Die Hyperbene wird verschoben
    #     Nicht signifikante Faktorstufe => keine ??nderung des Achsenabschnittes
    farm1$industry <- as.factor(farm1$industry)
    
    # Regression mit 2 Parametern analysieren!
    scatter3d(y ~ x1 + x2, data=,  axis.scales=FALSE)
    
  }
  
}
# Mathematik der Finanzm??rkte
function(){
  r <- c(0.0288, 0.0362, 0.0430, 0.0493, 0.0549, 0.0599) #Bsp
  s <- c(0.04,0.03) # in 100tel schreiben
  
  diskontfaktor <- function(prozent,dauer){
    return((1/(1+prozent))^dauer)
  }
  diskontfaktor(s,1:length(s)) #normaler diskontfaktor
  diskontfaktor(r,1) #normaler diskontfaktor
  
  # Kassazins
  sfc <- function(rates){
    sapply(1:length(rates), function(i) (prod(rates[1:i]+1)^(1/i))-1)
  }
  (s <- sfc(r))
  
  fij <- function(si,i,j){
    (((1+si[j])^j / (1+si[i])^i)^(1/(j-i)))-1
  }
  (f <- fij(s, von, bis))
  
  (sstrich <- fij(s,1,1:length(s)))
  
  # Short-Rates
  rfc <- function(s){
    x <- sapply(1:length(s), function(i) (fij(s,i,i+1)))
    return(c(s[1],x[1:length(x)-1]))
  }
  rfc(s)
  
  # lambda R??ckverzinsung
  bond.price = function(n,lambda,c=0,m=1,N=100) {
    lambda[lambda > 1] = 0.01*lambda 
    q = 1/(1+lambda/m)
    cPer = c/m
    nPer = floor(n*m)
    #print(paste(q,cPer,nPer))
    price = cPer*q*(1-q^nPer)/(1-q) + N*q^nPer
    round(price,2)
  }
  bond.coupon = function(price, lambda, n=1, m=1, N=100, annualized=T) {
    if (lambda > 1) { lambda = 0.01*lambda }
    nPer = floor(n*m)
    #print(nPer)
    c = price*lambda / 
      ( 1 - 1 / (1 + lambda/m)^nPer )
    if (!annualized) {
      c = c/m
    }
    c
  }
  bond.lambda = function(price,n,c=0,m=1,N=100, percentage=F) {
    library(polynom)
    
    nPer = n*m
    cPer = c/m
    p = polynomial( c(-price, rep(cPer,nPer-1), cPer+N) )
    d = solve(p)
    d = max(Re(d[Im(d)==0]))
    lambda = m*(1/d-1)
    if (percentage) { lambda = 100*lambda }
    lambda
  }
  ir.futures = function (i,j,ir.spot) {
    sr.len = length(ir.spot)
    i.len = length(i)
    j.len = length(j)
    
    if ( j[length(j)]>sr.len ) { 
      warning("Future period truncated according to spot rates") 
      j = j[j<=sr.len]
      j.len = length(j)
    }
    if ( i[length(i)]>sr.len ) { 
      warning("Future period truncated according to spot rates") 
      i = i[i<=sr.len]
      i.len = length(i)
    }
    if (i[i.len] > j[j.len] ) {
      warning("Index vector i truncated because of incompatibilities between 
    i and j")
      i = i[i<=j(j.len)]
    }
    if (j[1] == 0) {
      warning("Index j should not be less than 1")
      j = j[j>0]
    }
    
    x = (1+ir.spot[j])^j
    #print(x)
    y = 1/(1+ir.spot[i[i>0]])^i[i>0]
    if (i[1]==0) { y = c(1,y) }
    
    #print(y)
    res = t(as.matrix(x) %*% t(as.matrix(y)))
    dimnames(res) = list(i,j)
    #print(res)
    idx = t(matrix(rep(j,length(i)),length(j),length(i)))-i
    dimnames(idx) = list(i,j)
    idx[idx<0] = NA
    idx[idx==0] = 1
    #print(idx)
    
    res^(1/idx)-1
  }
  bond.duration = function(n,lambda,c=0,m=1,N=100, type="normal") {
    lambda[lambda > 1] = 0.01*lambda 
    bw = bond.price(n,lambda,c,m,N)
    q = 1/(1+lambda/m)
    cPer = c/m
    nPer = 1:floor(n*m)
    d = q^nPer
    #print(paste(d))
    dur = ( cPer*sum(nPer*d)/m + n*N*q^(m*n) )/bw
    if (type=="modif") {
      dur = dur*q
    }
    round(dur,2)
  }
  # type = ("normal"|"modif")
  
  # Zahlungsstr??me eingeben
  cfStream.npv = function(cf, irSpot, per=1, digits=2) {
    nPer = 1:length(cf)
    ir = irSpot[1:length(cf)]
    ir[ir > 1] = 0.01*ir 
    discount = 1/(1+ir/per)^nPer
    round(sum(cf*discount),digits)
  }
  cfStream.duration = function(cf, irSpot, per=1, type="quasimodif") {
    nPer = 1:length(cf)
    ir = irSpot[nPer]
    ir[ir > 1] = 0.01*ir 
    irp1inv = 1/(1+ir)
    
    # Quasi-modifizierte oder Fisher-Weil Duration
    if (type=="quasimodif") {
      ee = nPer + 1
      dur = sum( nPer/per*cf*irp1inv^(ee) ) / sum(cf*irp1inv^nPer)
    } else if (type =="Fisher-Weil") {
      ee = nPer
      dur = sum( nPer/per*cf*irp1inv^(ee) ) / sum(cf*irp1inv^nPer)
    } else if (type =="continuous") {
      ee = nPer
      dur = sum( nPer/per*cf*exp(-ir*(nPer/per)) ) / sum(cf*exp(-ir*(nPer/per)))
    } else {
      stop("Unkown type of duration")
    }
    dur
  }
  # type = "quasimodif|Fisher-Weil|continuous"
  
  (cf.b1 <- c(rep(6,11),106)) # die letzte Zahlung mit der R??ckgabe zusammenz??hlen
  (s <- c(7.67, 8.27, 8.81, 9.31, 9.75, 10.16, 10.52, 10.85, 11.15, 11.42,11.67, 11.89))
  cfStream.npv(cf.b1,s)
  cfStream.duration(cf.b1,s,type = "quasimodif")
  
  
  
  # Konvexit??t berechnen (Immunisierung mit sovle Rechner)
  function(){
    #Funktion f??r die Konvexit??t
    Konvex = function(Lambda, P, C, N, n){
      
      Disc1 = 1/(P*(1+Lambda)^2)
      temp = 0
      
      for(i in 1:n){
        
        diskont1 = (1+Lambda)^(-i)
        if(i == n){
          c=C+N
        }else{
          c = C
        }
        Z??hler = (i * (i+1)) * c
        
        res = Z??hler*diskont1
        temp = temp + res
        
      }
      
      return(Disc1 * temp)
      
    }
    Konvex(Lambda = Umlaufrendite, P = Preis, C = Coupon, N = 100, n = Laufzeit)
    #Anleihe 1 bis 3
    Konvex(Lambda = 0.045, P = 103.3, C = 5, N = 100, n = 8)
    # 52.98409
    Konvex(Lambda = 0.045, P = 64.39, C = 0, N = 100, n = 10)
    # 100.7346
    Konvex(Lambda = 0.045, P = 108.97, C = 7, N = 100, n = 4)
    # 16.13463
    
    #Verpflichtung
    Konvex(Lambda = 0.045, P = 80245.1, C = 0, N = 100000, n = 5)
    # 27.4719
    
    # eintragen der Daten!
    A = matrix(c(103.3,64.39,108.97,
                 103.3*6.81,64.39*10,108.97*3.64,
                 52.98,100.74,16.13), ncol = 3, nrow = 3, byrow = T)
    A
    
    b = c(80245.1, 80245.1 * 5, 27.47)
    
    #Ax=b
    
    A_inv = solve(A)
    
    A_inv %*% b
    solve(a = A, b = b)
  }
  
  
  # Kovarianz-Matrix
  A <- cbind(x1,x2,x3)
  cov_mat <- var(A[,1:ncol(A)]); cov_mat
  
  # Portfoliovarianz NOCH NICHT GETESTED SKRIPT S 105
  Por_Var <- function(w,cov_mat){
    w <- t(w)
    Por_Var <- t(w) %*% cov_mat %*% w
    return(Por_Var)
  }
  
  kronecker(A_mat,A_mat)
  
  # Erwartungswert und Varianz eines optimierten Portfolios
  #   cov_mat <-E
  #   r_ptf <- 1
  EV <- function(cov_mat, r, r_ptf){
    A <- cov_mat
    A <- cbind(cov_mat,t(t(-r)),t(t(rep(-1,nrow(cov_mat)))))
    A <- rbind(A,c(r,0,0),c(rep(1,c(nrow(cov_mat))),0,0))
    
    b <- t(t(c(rep(0,nrow(cov_mat)),r_ptf,1))); b
    solve1 <- t(solve(A,b))
    
    solve1 <- as.numeric(solve1)
    names(solve1) <- c(rep("weight",nrow(cov_mat)),"lambda","mu"); solve1
    
    return(solve1)
  }
  # Bei minimaler Varianz
  EVm <- function(cov_mat, r, r_ptf){
    A <- cov_mat; A
    A <- cbind(A,rep(-1,nrow(cov_mat))); A
    A <- rbind(A,c(rep(1,c(nrow(cov_mat))),0)); A
    
    b <- c(rep(0,nrow(cov_mat)),1); b
    solve1 <- solve(A,b)
    
    solve1 <- as.numeric(solve1)
    names(solve1) <- c(rep("weight",nrow(cov_mat)),"lambda","mu"); solve1
    
    return(solve1)
  }
  
  
  # Begriffe
  # Preisvariabilitaet (Volatilitaet)
  
  # pft Berechnen und aufzeichnen
  # Mit Leerverk??ufen gehen die Geraden ins Unendlichen nach Oben und Unten heraus
  PTF(c(0.1,0.18),c(0.15,0.3),0)
  #   renditen <- c(0.1,0.18,0.2)
  #   sigma <- c(0.15,0.3,0.2)
  #   korrelation <- 0
  PTF <- function(renditen,sigma,korrelation){
    renditeA <- renditen[1]
    renditeB <- renditen[2]
    sigmaA <- sigma[1]
    sigmaB <- sigma[2]
    
    ptfRendite <- function(alpha, korrelation) {
      sigmaAB <- korrelation*sigmaA*sigmaB
      
      renditeP <- alpha*renditeA + (1-alpha)*renditeB
      varP <- sigmaA^2*alpha^2 + 2*sigmaAB*alpha*(1-alpha) + sigmaB^2*(1-alpha)^2
      sigmaP <- sqrt(varP)
      
      return(c(renditeP, sigmaP))
    }
    
    alpha <- seq(0,1,0.02)
    yachsenabschnitt <- (renditeA*sigmaB+renditeB*sigmaA)/(sigmaA+sigmaB)
    
    mat <- matrix(rep(NA, 51*2), nrow=51)
    #View(mat)
    
    for(i in 1:51) {
      mat[i,1] <- ptfRendite(alpha[i], korrelation)[1]
      mat[i,2] <- ptfRendite(alpha[i], korrelation)[2]
    }
    
    plot(mat[,2], mat[,1], xlab="Sigma", ylab="Rendite", type="l", xlim=c(0,0.3))
    lines(c(ptfRendite(0,1)[2], 0), c(ptfRendite(0,1)[1], yachsenabschnitt)) 
    lines(c(ptfRendite(1,1)[2], 0), c(ptfRendite(1,1)[1], yachsenabschnitt))
    
    x1 <- min(mat[,2]); x2 <- mat[order(mat[,2])[1],1]; x2.5 <- yachsenabschnitt;
    w1 <- (x2-renditeB)/(renditeA-renditeB); w2 <- 1-w1; x3 <- c(x1,x2,x2.5,w1,w2)
    names(x3) <- c("min VARIANZ","RENDITE min Var","RENDITE Risikofrei","GEWICHTUNG 1","GEWICHTUNG 2")
    points(x1,x2,col=2,lwd = 3,pch=4)
    text(x1,x2,labels = "min. Var",pos = 4,col=2)
    return(x3)
  }
  
  #Mehrere Anleihen
  PTFX <- function(renditen,sigma,korrelation){
    PTF2 <- function(renditen,sigma,korrelation){
      renditeA <- renditen[1]
      renditeB <- renditen[2]
      sigmaA <- sigma[1]
      sigmaB <- sigma[2]
      
      ptfRendite <- function(alpha, korrelation) {
        sigmaAB <- korrelation*sigmaA*sigmaB
        
        renditeP <- alpha*renditeA + (1-alpha)*renditeB
        varP <- sigmaA^2*alpha^2 + 2*sigmaAB*alpha*(1-alpha) + sigmaB^2*(1-alpha)^2
        sigmaP <- sqrt(varP)
        
        return(c(renditeP, sigmaP))
      }
      
      alpha <- seq(0,1,0.02)
      yachsenabschnitt <- (renditeA*sigmaB+renditeB*sigmaA)/(sigmaA+sigmaB)
      
      mat <- matrix(rep(NA, 51*2), nrow=51)
      #View(mat)
      
      for(i in 1:51) {
        mat[i,1] <- ptfRendite(alpha[i], korrelation)[1]
        mat[i,2] <- ptfRendite(alpha[i], korrelation)[2]
      }
      
      lines(mat[,2], mat[,1], xlab="Sigma", ylab="Rendite", type="l", xlim=c(0,0.3))
      
      #     x1 <- min(mat[,2]); x2 <- mat[order(mat[,2])[1],1]; x2.5 <- yachsenabschnitt; x3 <- c(x1,x2,x2.5)
      #     names(x3) <- c("min Varianz","Rendite bei min Var","Rendite min Var")
      #     points(x1,x2,col=2,lwd = 3,pch=4)
      #     text(x1,x2,labels = "min. Var",pos = 4,col=2)
      #     return(x3)
    }
    plot(x=0,y=0,xlim = c(0,max(sigma)),ylim = c(min(renditen),max(renditen)),col="white")
    for(j in 1:length(renditen)){
      for(i in 1:length(renditen)){
        if(i>j){
          PTF2(c(renditen[j],renditen[i]),c(sigma[j],sigma[i]),korrelation)
        }
      }
    }
  }
  PTFX(c(0.1,0.18,0.15,0.05),c(0.15,0.3,0.4,0.35),0)
  PTFX(c(0.1,0.18,0.15),rep(0.3,3),0)
  
}
# Datum in R
function(){
  #Der Vektorstrom wird in einen Zeitstrom umgewandelt
  data <- ts(data, start=c(startJahr,1), frequency=Anzahl pro Jahr)
  
  Aust.10Y
  
  rendite <- function(Daten,ZeilenReihe){
    sapply(1:(dim(Daten)[1]-1), function(i) Daten[i,ZeilenReihe]-Daten[i+1,ZeilenReihe])
  }
  rendite(Aust.10Y,4)
  
  head(Aust.10Y)
  Aust.10Y[,"Settle"]
  
  head(Aust.10Y)
  ?ts()
  (data <- ts(Aust.10Y ,start = 1990,frequency = 360))
  ts.plot()
  
  plot(diff(diff(log(y), lag=1)))
  ?diff()
  diff(1:100, lag = 15)
  ?stl() # aufspalten der einzelnen Effekte von ts()
  
  date()
  Sys.Date()
  as.character()
  
  # Jahre aufz??hlen
  nYears(valueDate, -1)
  
  # erh??hung des Tages und benennung des Orts
  nBusinessDays(Sys.Date(),2,"Zurich")
  
  'ACT/360' = in Tagen und Monaten
  '30E/360' = in Jahren
  
  # Ausschnitt aus Daten anzeigen
  window(Daten,von,bis)
}