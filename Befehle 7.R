function("OEKO 3"){
  
  function("ARCH / GARCH"){
    
    
    ######## ######## ######## ########.
    # QUESTIONS
    # 
    # Error Analysis: Was bedeutet shape?
    # 
    # Log Likelihood: google es nochmals und schreib es in die Formelzusammenfassung
    #
    # Was sagt: LM Arch Test
    #
    ######## ######## ######## ########.
    
    
    library(tseries)
    library(fGarch)
    
    #### dat sind log-returns ####.
    dat <- rnorm(1000,mean = 2,sd = 3)
    
    ### einfacher GARCH Befehl
    x.garch <- garch(dat,order=c(1,1))
    summary(x.garch)
    ### erweitereter GARCH Befehl
    y.garch_1 <- garchFit(~garch(2,0), data = dat, delta=2, include.delta = F, include.mean=T, trace=F)
    # Delta schätzen = include.delta=T     if(include.delta == T) include.delta = T #Delta wird weg gelassen
    # delta=2, include.delta = F
    y.garch_1 <- garchFit(~arma(1,1) + garch(1,1), data = dat, cond.dist = "std", include.mean = T,trace=F)
    y.garch_1 <- garchFit(~garch(1,1), data = dat, cond.dist = "std", include.mean = T,trace=F)
    summary(y.garch_1)
    
    ### APARCH
    y.aparch_delta <- garchFit(~arma(1,1)+aparch(1,1),data=y,include.delta=T,include.mean=T)
    summary(y.garch_delta)
    
    ## Modellgleichung: sigma_t^2 = omega + alpha1 * epsilon_t-1^2 + alpha2 * epsilon_t-2^2
    ## mu + Epsilon_t = log(X_t) + log(X_t-1)
    ## d entspricht omega aus den Formeln
    
    
    
    ### Error Analysis ###
    ## Bezeichnung: Ergebnis des Fittes
    # Summe der Aplha grösser 1 => Parameter restrikion wird nicht erfüllt
    # Summe der Betas über 1 => FALSCH: GARCH Prozess nicht stationär
    # P-Value grösser als 0.05 => Zahl unnötig, Zahl ist sehr klein => sie hat grossen Einfluss
    
    
    ### Standardised Residuals Tests ###
    # Jarque-Bera Test   =>  Hypothesentest für die Normalverteilung (p-Value klein => keine Normalverteilung)
    # Shapiro-Wilk Test  =>  Hypothesentest für die Normalverteilung
    #   Nullhypothese    =>  U ist Normalverteilt (Die Schiefe und Tails werden kontrolliert)
    
    ## Ljung-Box
    # Falls einer der Ljung-Box Test unter 0.05 kommt: eine Abhängigkeit existiert immernoch
    # Wenn die Ljung-Box Test R (P-Value) unter 0.05 sind => ARMA Modell erhöhen
    #     Kontrolle ob ARMA genauer erklärt werden muss
    # Wenn die Ljung-Box Test R^2 (P-Value) unter 0.05 sind => ARCH/GARCH Ordnung erhöhen
    #     Kontrolle ob immernoch Vola Cluster existieren
    
    ## Diagnosestatisktiken: überprüfen ob u_t symetrie und ob es zu viele Ausreisser hat.
    ## keie Autokorrelation wird mit den ersten Ljung-Box Tests (R) geprüft
    
    
    ### Information Criterion Statistics ###
    ## Wie gut ist das Modell?
    # BIC sollte möglichst klein sein (unter Null)
    #     Wenn das Modell zu komplex oder zu simpel ist, wird BIC steigen.
    # AIC sollte möglichst klein sein
    #     Wählt im Gegensatz zu BIC jedoch lieber zu komplizierte Modelle. BIC ist wichtiger.
    
    
    # stationär => Varianz ist endlich und konstant
    
    # Die bedingte Vola konvergiert irgendwann gegen die unbedingte Volatität wenn die die Summe der Beta kleiner
    # gleich 1 sind erfüllt sind sonst divergiert es Konvergenzgeschwindigkeit der prognostizierten Volatilität
    # hängt von alpha ab.
    # alpha bestimmt das Memorie...
    # Das Vertrauensintervall steigt, wenn man in einer tief Vola phase ist.
    # Es steigt, wenn man in einer Hohen ist.
    
    function("zeitreihen Analyse"){
      dowjones<-read.table(paste(Pfad.dat,"dowjones.txt",sep=""),sep="\t",header=F, na.strings = "NA")
      x<-as.ts(dowjones)
      layout(matrix(c(1,3,1,4,2,5),ncol=3))
      ts.plot(x); y<-diff(log(x)); ts.plot(y)
      ### ACF von log-returns und quadrierten log-returns
      acf(y,main = "data")             # Strukturen sind nur schlecht sichtabr   (ARMA bezieht sich darauf)
      acf(y^2,main="(data)^2")         # Strukturen werden deutlicher            (GARCH bezieht sich darauf)
      acf(abs(y),main="abs(data)")     #                                         (GARCH bezieht sich darauf)
      par(mfrow=c(1,1))
    }
    
    
    function("Modelle Schätzen"){
      # GARCH(1,1)-Modell schätzen
      y.garch_11 <- garchFit(~garch(1,1),data=dat,delta=2,include.delta=F,include.mean=F,trace=F)
      summary(y.garch_11)   # shiehe Beschreibung von oben
      y.garch_11@sigma.t    # Bedingte standardabweichung # ts.plot(y.garch_11@sigma.t)
      y.garch_11@residuals  # Residuen epsilon_t (das sind nicht die u_t's)
      y.garch_11@fit$coef   # Modell-Parameter
      
      # If sum of coefficients is smaller 1 => estimated.var is better than var(data)
      sum(y.garch_11@fit$coef)
      unbedingte.var <- omega / (sum(alpha) - sum(beta))  # Parameter aus GARCH-Modell
      estimated.var <- (y.garch_11@fit$coef[1]/(1-y.garch_11@fit$coef[2]-y.garch_11@fit$coef[3]))
      var(dat)
      
      # standardized residuals
      ts.plot(residuals(y.garch_11, standardize=T))
      ts.plot(y.garch_11@residuals)
      
      eps<-y.garch_11@residuals
      u<-eps/y.garch_11@sigma.t
      ts.plot(eps/u)
      
      # wird die bedingt Varianz momentan untershcätz oder überschätzt
      y.garch_11@sigma.t[length(y.garch_11@sigma.t)]
      var(dat)
    }
    
    
    function("Prognose"){
      y <- dat # die Zeitreihe
      y.garch_11 # das Modell
      fors<-1000 # Anzahl Schritte
      std_ml<-sqrt(y.garch_11@fit$coef[2]/(1-y.garch_11@fit$coef[3]-y.garch_11@fit$coef[4]))
      y.garch_11.pred<-predict(y.garch_11,n.ahead=fors)
      # Plot Garch Prognoseintervalle vs Standard Prognoseintervalle
      ts.plot(cbind(y.garch_11.pred[,1],
                    y.garch_11.pred[,1] + 2*y.garch_11.pred[,3],
                    y.garch_11.pred[,1] - 2*y.garch_11.pred[,3],
                    rep(mean(y),fors) + 2*sqrt(var(y)),
                    rep(mean(y),fors) - 2*sqrt(var(y)),
                    y.garch_11.pred[,1]+2*std_ml,y.garch_11.pred[,1]-2*std_ml),lty=c(1,2,2,3,3,4,4)
      )
    }
    
    
    #### Nicht symetrische GARCH
    # wenn xi = 0 ist, dann ist die Verteilung nicht schief
    
    # APARCH
    # ...kann den TGARCh EGARCH und normalen GARCH darstellen
    # siehe Zusamenfassung
    
  }
  
  function("State Space"){
    
    function("KF_level"){
      
      ############### input
      # 
      ###### parma
      # R = 1, Q = 0 => standart
      R<-1
      Q<-0
      xi10<-y[1]
      parma<-c(sqrt(R),sqrt(Q),xi10)
      
      ###### y Timeseries
      # input Timeseries
      y
      
      ###### opti
      # T/F Logical
      # optimieren (T), nicht optimieren (F)
      # if True then KF_level returns the criterion value only;
      # otherwise it returns also the state-vector and the variances
      opti<-F
      
      ###### specify the optimization criterion
      ## outofsample
      # if T then the out-of-sample forecast errors are used; otherwise the in-sample forecast errors
      # T/F Logical
      outofsample<-T
      ## maxlink
      # maximumlikelihood (T), chiquadrat (F)
      # T/F Logical
      maxlik<-T
      
      ###### P10
      # Variance of y, mutiplied to a big number
      # muss einfach eine grosse Zahl sein
      P10<-100000*var(y)
      
      
      
      
      ## R und Q können beliebig sein, nur das Verhältnis zählt
      # Q/R  gross    -   hohe adaptivität      (>0.1 ist adaptiv)
      # Q/R  klein    -   kleine adaptivität    (<0.1 ist starr, 10E-6 => starr)
      # R is a scalar, Q is a diagonal matrix
      
      ts.plot(cumsum(y))
      (y_obj <- KF_level(parma,y,opti,outofsample,maxlik,P10))
      
      ############### output
      # Mean over time
      y_obj$xitt
      ts.plot(y_obj$xitt)
      # Variance
      y_obj$Ptt
      plot(y_obj$Ptt)
      
    }
    KF_level<-function(parma,y,opti,outofsample,maxlik,P10) {
      len<-length(y)
      # note that we parametrize the variances such that they are always positive!
      R<-parma[1]^2
      Q<-parma[2]^2
      xi10<-parma[3]
      Pttm1<-0:len
      Ptt<-1:len
      Pttm1[1]<-P10
      xtt<-xi10
      # initialization of loglikelihood
      logl<-0.
      # we collect the state vectors xi_{t|t-1} and xi_{t|t}
      xittm1<-1:len
      xitt<-xittm1
      # Start of KF recursions
      for (i in 1:len)    #i<-2
      {
        M<-Q
        K<-Pttm1[i]
        # The Kalman Gain
        Kg<-1/(K+R)
        # epsilon
        epshatoutofsample<-y[i]-xtt
        # xi_{t|t-1}
        xittm1[i]<-xtt
        # up-date xi_{t|t}
        xtt<-xtt+Pttm1[i]*Kg*epshatoutofsample
        # in-sample forecast error (after y_t has been observed)
        epshatinsample<-y[i]-xtt
        # compute P_{t|t}
        Ptt[i]<-Pttm1[i]-Pttm1[i]*Kg*Pttm1[i]
        # compute P_{t+1|t}
        Pttm1[i+1]<-Ptt[i]+M
        # trace xi_{t|t}
        xitt[i]<-xtt
        # The optimization criterion
        if (outofsample)
        {
          if (maxlik)
          {
            logl<-logl+log(K+R)+epshatoutofsample^2*Kg
          } else
          {
            logl<-logl+epshatoutofsample^2
          }
        } else
        {
          if (maxlik)
          {
            logl<-logl+log(K+R)+epshatinsample^2*Kg
          } else
          {
            logl<-logl+epshatinsample^2
          }
        }
      }
      if (opti)
      {
        return(logl/len)
      } else
      {
        return(list(logl=logl/len,
                    xitt=xitt,xittm1=xittm1,Pttm1=Pttm1,Ptt=Ptt))
      }
    }
    
    function("Kalman_gen"){
      ###### parma
      # Nur Q und R benutzen!
      parma <- c(sqrt(Q),sqrt(R))
      
      ###### y
      # y is the data
      y
      
      ###### opti
      # T/F Logical
      # optimieren (T), nicht optimieren (F)
      # if True then KF_level returns the criterion value only;
      # otherwise it returns also the state-vector and the variances
      opti<-F
      
      ###### specify the optimization criterion (nicht ändern)
      ## outofsample
      # if T then the out-of-sample forecast errors are used; otherwise the in-sample forecast errors
      # T/F Logical
      outofsample<-T
      ## maxlink
      # maximumlikelihood (T), chiquadrat,least-squares (F)
      # T/F Logical
      maxlik<-T
      
      ###### xi10 and P10
      # variance start value
      P10<-100000*var(y)
      # state-vector start value
      xi10<-y[1]
      
      # Fm and H correspond to the system matrices
      Fm<-1
      H<-1
      
      # time_varying<-T then the system matrix H is filled with the data
      time_varying<-F
      
      # In case of a REGRESSION enter datamatrix of the parameters
      # it corresponds to the series of explanatory variables.
      x_data<-NULL
      
      
      (obj <- KF_gen_modified(parma,y,opti,outofsample,maxlik,xi10,P10,Fm,H,time_varying,x_data))
      
      ############### output
      # Mean over time
      obj$xitt # col1 = Level, col2 = slope
      ts.plot(obj$xitt)
      # Variance
      obj$Ptt
    }
    KF_gen<-function(parma,y,opti,outofsample,maxlik,xi10,P10,Fm,H,time_varying,x_data) {
      len<-length(y)
      if (length(parma)>2)
      {
        Q<-diag(parma[1:(length(parma)-1)]^2)
      } else
      {
        Q<-as.matrix(parma[1]^2)
      }
      R<-parma[length(parma)]^2
      Pttm1<-array(dim=c(dim(Q),len+1))
      Pttm1[,,1:(len+1)]<-P10
      Ptt<-array(dim=c(dim(Q),len))
      xttm1<-xi10
      logl<-0.
      xittm1<-matrix(nrow=len,ncol=(length(parma)-1))
      xittm1[1,]<-xi10
      xitt<-xittm1
      # If time_varying==T then we fill data into H (SSM with time-varying coefficients)
      if (time_varying)
      {
        if (is.null(x_data))
        {
          # For an autoregressive model we fill in past y's
          H<-c(y[1:dim(Q)[2]])
          # We need the first y[1:p] in H: therefore the first equation will be for t=p+1
          anf<-dim(Q)[2]+1
        } else
        {
          # For a regression model we fill in the explanatory data
          H<-x_data[1,]
          anf<-1
        }
      } else
      {
        anf<-1
      }
      # Kalman-Recursion: starts in i=dim(Q)[2] for a time series model
      # and in i=1 for a regression model
      
      for (i in anf:len)        #i<-1     H<-c(1,0)       xitt[,2]
      {
        # Kalman-Gain
        He<-(H%*%(Pttm1[,,i]%*%H))[1,1]+R
        epshatoutofsample<-y[i]-(H%*%xttm1)[1,1]
        xittm1[i,]<-xttm1
        xtt<-xttm1+Pttm1[,,i]%*%H*epshatoutofsample/He
        epshatinsample<-y[i]-(H%*%xtt)[1,1]
        xitt[i,]<-xtt
        xttm1<-Fm%*%xtt
        Ptt[,,i]<-Pttm1[,,i]-((Pttm1[,,i]%*%H)%*%(H%*%Pttm1[,,i]))/He
        Pttm1[,,i+1]<-Fm%*%Ptt[,,i]%*%t(Fm)+Q
        if (time_varying)
        {
          if (is.null(x_data))
          {
            # For an autoregressive model we fill past y's in H
            H<-c(y[i-dim(Q)[2]+1:(dim(Q)[2])])
          } else
          {
            # For a regression model we fill the explanatory data in H
            H<-x_data[min(len,i+1),]
          }
        }
        # Here we specify the optimization criterion: least-squares
        # or maximum likelihood, in-sample or out-of-sample
        if (outofsample)
        {
          if (maxlik)
          {
            logl<-logl+log(He)+epshatoutofsample^2/He
          } else
          {
            logl<-logl+epshatoutofsample^2
          }
        } else
        {
          if (maxlik)
          {
            logl<-logl+log(He)+epshatinsample^2/He
          } else
          {
            logl<-logl+epshatinsample^2
          }
        }
      }
      if (opti)
      {
        return(logl/len)
      } else
      {
        return(list(logl=logl/len,xitt=xitt,xittm1=xittm1,Ptt=Ptt,Pttm1=Pttm1))
      }
    }
    
    KF_gen_modified<-function(parma,y,opti,outofsample,maxlik,xi10,P10,Fm,H,time_varying,x_data) {
      len<-length(y)
      if (length(parma)>2)
      {
        # MODIFICATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Q<-diag(c(parma[1],rep(0,length(parma)-2)))
      } else
      {
        Q<-as.matrix(parma[1]^2)
      }
      R<-parma[length(parma)]^2
      Pttm1<-array(dim=c(dim(Q),len+1))
      Pttm1[,,1:(len+1)]<-P10
      Ptt<-array(dim=c(dim(Q),len))
      xttm1<-xi10
      logl<-0.
      xittm1<-matrix(nrow=len,ncol=(length(parma)-1))
      xittm1[1,]<-xi10
      xitt<-xittm1
      # If time_varying==T then we fill data into H (SSM with time-varying coefficients)
      if (time_varying)
      {
        if (is.null(x_data))
        {
          # For an autoregressive model we fill in past y's
          H<-c(y[1:dim(Q)[2]])
          # We need the first y[1:p] in H: therefore the first equation will be for t=p+1
          anf<-dim(Q)[2]+1
        } else
        {
          # For a regression model we fill in the explanatory data
          H<-x_data[1,]
          anf<-1
        }
      } else
      {
        anf<-1
      }
      # Kalman-Recursion: starts in i=dim(Q)[2] for a time series model
      # and in i=1 for a regression model
      
      for (i in anf:len)        #i<-1     H<-c(1,0)       xitt[,2]
      {
        # Kalman-Gain
        He<-(H%*%(Pttm1[,,i]%*%H))[1,1]+R
        epshatoutofsample<-y[i]-(H%*%xttm1)[1,1]
        xittm1[i,]<-xttm1
        xtt<-xttm1+Pttm1[,,i]%*%H*epshatoutofsample/He
        epshatinsample<-y[i]-(H%*%xtt)[1,1]
        xitt[i,]<-xtt
        xttm1<-Fm%*%xtt
        Ptt[,,i]<-Pttm1[,,i]-((Pttm1[,,i]%*%H)%*%(H%*%Pttm1[,,i]))/He
        Pttm1[,,i+1]<-Fm%*%Ptt[,,i]%*%t(Fm)+Q
        if (time_varying)
        {
          if (is.null(x_data))
          {
            # For an autoregressive model we fill past y's in H
            H<-c(y[i-dim(Q)[2]+1:(dim(Q)[2])])
          } else
          {
            # For a regression model we fill the explanatory data in H
            H<-x_data[min(len,i+1),]
          }
        }
        # Here we specify the optimization criterion: least-squares
        # or maximum likelihood, in-sample or out-of-sample
        if (outofsample)
        {
          if (maxlik)
          {
            logl<-logl+log(He)+epshatoutofsample^2/He
          } else
          {
            logl<-logl+epshatoutofsample^2
          }
        } else
        {
          if (maxlik)
          {
            logl<-logl+log(He)+epshatinsample^2/He
          } else
          {
            logl<-logl+epshatinsample^2
          }
        }
      }
      if (opti)
      {
        return(logl/len)
      } else
      {
        return(list(logl=logl/len,xitt=xitt,xittm1=xittm1,Ptt=Ptt,Pttm1=Pttm1))
      }
    }
    
    # Einfluss auf Xitt
    #   - out of sample is sinnvoll, insample ist sinnlos
    #   - maximum likelihood ist interesannter als kleinste quadrate
    #       kleinster quadrate kann nure das Verhältniss von Q und R berechnen, aber Q und R selber geht nur mit maximum likelihood
    
    
    # optimising for parma
    function(){
      ## must be given
      # 
      # function: KF_level
      # y
      # opti
      # outofsample
      # maxlik
      # P10
      
      ## KF_gen
      # parma<-c(rep(0,dim(diag(Q))),sqrt(R))
      #
      ## KF_level
      # parma <- sqrt( c(var(y)/len,var(y)) )
      
      objopt<-nlminb(
        start=parma
        ,objective=KF_level
        ,y=y
        ,opti=opti
        ,outofsample=outofsample
        ,maxlik=maxlik
        ,P10=P10
      )
      parma<-c(abs(objopt$par[1:2]),objopt$par[3]) # WICHTIG abs Umwandeln
      
      # The first two parameters are (positive) standarddeviations
      parma<-c(abs(parma[1:2]),parma[3])
      
      # 3. Q/R ratio (model is fairly adaptive)
      parma[2]/parma[1]
      
      # 4.
      parma
      
    }
    
    # Plot the comparison
    function(){
      # input
      y
      obj$xitt
      obj$Ptt
      mu
      len
      
      
      ts.plot(obj$xitt,ylim=c(min(y),max(y)),col="blue",
              main=paste("maxlik=",maxlik,", outofsample=",outofsample,", R=",
                         round(parma[2],3),", Q=",round(parma[1],3),sep=""),
              xlab="",ylab="")
      lines(y,col="black")
      lines(mu,col="red")
      lines(cumsum(y)/(1:len),col="green")
      lines(obj$xitt-2*sqrt(obj$Ptt[1,1,]),col="blue",lty=2)
      lines(obj$xitt+2*sqrt(obj$Ptt[1,1,]),col="blue",lty=2)
      lines(rep(mu,len),col="red")
      mtext("Xi_{t|t}: 95% confidence intervall", side = 3, line = -2,at=len/2,col="blue")
      mtext("data", side = 3, line = -1,at=len/2,col="black")
      mtext("smooth level (target)", side = 3, line = -3,at=len/2,col="red")
      mtext("mean", side = 3, line = -4,at=len/2,col="green")
    }
    
    # Check mean-square performances
    #   of the arithmetic mean,
    #   of the state-space estimate and of the identity
    function(){
      mean((y-mu)^2)
      mean((obj$xitt[,1]-mu)^2)
      mean((cumsum(y)/(1:len)-mu)^2)
    }
    
  }
  
}

function("Risk Engineering"){
  library(QRM)
  library(R2HTML)
  library(FinTS)
  library(fBasics)
  library(copula)
  library("PerformanceAnalytics")
  library(ghyp)
  library(timeSeries)
  
  ## 9 & 10 extrem wichtig f??r Pr??fung
  
  DJ30daily <- getReturns(DJ, type = "continuous")
  par(mfrow = c(2, 3)); for (i in 1:5){ histPlot(DJ30daily[,i])}
  par(mfrow = c(2, 3)); for (i in 1:5){ qqnormPlot(DJ30daily[,i])}
  par(mfrow = c(1,1))
  colIds(DJ)  # colnames(DJ)
  
  function("Zeitdaten anpassen"){
    timeSequence(from = start(DJ30daily), to = end(DJ30daily), by = "month") # "week","month","quarter"
    
    # Zeitschrittgrösse ändern
    DJ30weekly <- aggregate(DJ30daily, by = by$quarter, FUN = sum)
    DJ30weekly2 <- daily2weekly(DJ30daily) # daily2monthly()
    
    # Amerikanisch Schreibweise des Datum!!!
    window(DJ, timeDate("01/01/1993"), timeDate("12/31/2000"))
  }
  
  function("Eigenschaften über ZeitreihenVerteilung auslesen (explore funktion)"){
    VolAnnual <- sqrt(colVars(series(DJ30daily))) # * sqrt(250) * 100
    Skewness <- colSkewness(DJ30daily)
    Kurtosis <- colKurtosis(DJ30daily) # kurtosis(DJ30daily)
    
    AutocorTest(Lt, type = "Ljung-Box")$p.value
    normalTest(Lt, method = "jb", na.rm = TRUE)@test$p.value
    
    # Alles zusammen von einer Zeitreihe
    explore <- function(data){
      explo <- list()
      explo[["Standartabweichung"]] <- sqrt(colVars(series(data))) # * sqrt(250) * 100
      explo[["Skewness"]] <- colSkewness(data)
      explo[["Kurtosis"]] <- colKurtosis(data) # kurtosis(DJ30daily)
      explo[["AutocorTest"]] <- AutocorTest(data, type = "Ljung-Box")$p.value
      explo[["normalTest"]] <- normalTest(data, method = "jb", na.rm = TRUE)@test$p.value
      return(explo)
    }
  }
  
  function("VaR (Value at Risk) / ES (Expected Shortfall)"){
    # Mittelwert, Standartabweichung, Quantile und ...
    
    # Die Datenreihe wird der Grösse nach sortiert, danach wird kontrolliert wie höufig, dass 
    
    # Value at Risk => mit einer Wahrscheinlichkeit von 95% (alpha) wird dieser Wert nicht unterschritten
    # Expected Shortfall => Der Schnitt wo der Wert im Schnitt liegen wird
    
    X <- sort(series(DJ[,1]))
    len <- length(X)
    alpha <- c(0.95, 0.99)
    
    VaR1 <- X[ceiling(alpha * len)]
    
    ES1 <- rep(NA,2)
    ES1[1] <- mean(X[ceiling(alpha[1] * len):len])
    ES1[2] <- mean(X[ceiling(alpha[2] * len):len])
    
    # par(mfrow=c(2,1))
    plot(X, (1:len) / len, type = "l")
    abline(v = c(VaR1, ES1), col = c(2, 2, 3, 3), lty = rep(1:2, 2))
    list("Value at Risk" = VaR1, "Expected Shortfall" = ES1)
    par(mfrow=c(1,1))
    
    ghyp()
    ghyp::qghyp(p = ,)
    
    VaR(getReturns(DJ[,3]),p = 0.95, method = "gaussian")
    ES(getReturns(DJ[,3]),p = 0.95)
  }
  
  function("Test for normalise and autocorrelation"){
    normalTest(returns[,1], method="jb", na.rm = TRUE)@test$p.value
    AutocorTest(returns[,1], type="Ljung-Box")$p.value
    AutocorTest(abs(returns), type="Ljung-Box")$p.value
  }
  
  # Portfolioverlust L_t+1
  # Vt = investitionssumme
  # wi = gewichte der Eingegebenen Zeitreihen
  # rf.delta = Data returns (X_t)
  Verlust <- function(Vt,wi,rf.delta, type="full"){
    x <- rf.delta
    is.ts <- FALSE
    if (class(x) == "timeSeries") {
      x <- series(x)
      is.ts <- TRUE
    }
    if (type == "linear") {
      res <- - Vt * (x %*% wi)
    } else {
      res <- - Vt * (exp(x) - 1) %*% wi
    }
    if (is.ts) {
      res <- timeSeries(pos = positions(rf.delta), data = res)
    }
    return(res)
  }
  type = "linear"
  
  function("Verteilungen fitten"){
    # Verteilungen Fitten
    fit.gaussuv(data = DJ30daily[,1])
    fit.tuv(data = DJ30daily[,1])
    fit.norm(rnorm(100000))
    fit.st(rt(10000,3))
    
    # GH-Verteilung Hyperbolische Verteilung 
    mv.hyp.symm = ghyp(sigma=matrix(c(1,-0.7,-0.7,1),ncol=2), mu=c(0,0), gamma=c(0,0), psi=1, chi=1, lambda=1)
    
    ## Verteilung anpassen
    # conddist - 
    mod4 <- garchFit(formula = ~ garch(1,1), data = ptf.loss1, cond.dist = "std", include.mean = TRUE,trace=FALSE)
  }
  
  function("Skewness und Kurtosis erstellen"){
    TwoNormDist <- function(data, sig1=1,sig2=1,mu1=0, mu2=0,p=0.5){
      p*dnorm(data,mean=mu1,sd=sig1) + (1-p)*dnorm(data,mean=mu2, sd=sig2)
    }
    
    par(mfrow=c(2,1),mar=c(3,5,2,4))
    plot(seq(-4,4,0.1),TwoNormDist(seq(-4,4,0.1),sig1=2,sig2=0.5),type="l",
         ylim=c(0,0.6),ylab="Density")
    lines(seq(-4,4,0.1),dnorm(seq(-4,4,0.1),sd=1),lty=2)
    lines(seq(-4,4,0.1),dt(seq(-4,4,0.1)*sqrt(2),df=4)*sqrt(2),lty=3)
    
    plot(seq(-4,4,0.1),TwoNormDist(seq(-4,4,0.1),mu1=-0.5,mu2=0.5,
                                   sig1=2,sig2=0.5),type="l",ylab="Density")
    lines(seq(-4,4,0.1),dnorm(seq(-4,4,0.1),sd=1),lty=2)
  }
  
  function("Kopula"){
    # ein Kopula ist:
    # Ist eine Verteilungsfunktionö im d-dimensionalen Eiheitswürfel mit uniformen Randverteilungen
    
    # Bedeutung
    # C(0.1,0.2) = 5%
    # Die Warhsch. das die Werte der 1. Variable kleiner 0.1 und der 2. kleiner 0.2 sind steht bei 5%
    
    
    
    ### rGleichverteilt mit korrelation unter den erzeugten daten
    ### Funktion zur Simulation einer GaussCopula
    # n   - anzahl Messungen
    # dim - anzahl dimensionen
    # rho - korrelation
    rcopula.gauss =function(n, dim, rho){
      library(ghyp)
      Sigma<-matrix(rho,dim,dim)
      diag(Sigma)<-1
      N<-gauss(mu=rep(0,dim),sigma=Sigma)
      X<-rghyp(n,N)
      u<-apply(X,2,FUN=pnorm)
      return(u)
    }
    
    ### GEV fitten
    ### Berechnen von H_xi,mu,sigma
    # benutzen von gev()
    GEV.week <- gev(series(sp500.week),block = T)
    rlevel.year <- rlevel.gev(GEV.week, k.blocks = 10)
    
    
    ### Fit Pareto Modell
    ### verallgemeinertes Pareto modell
    out <- gpd(data,10) # Schwellwert 10M
    
    # Plot the results
    evirplot(out,1)
    # Riskmeasure
    riskmeasures(out, c(0.99,0.995,0.999))
    # Schätzer VaR
    gpd.q(out.plot, 0.995)
    # Schätzer ES
    gpd.sfall(out.plot, 0.995)
    
  }
  
  library(evir)
  setwd("~")
  source("Studium/6. Semester/RiskEngineering - Moodle/R/Funktionen.R")
  # Wichtig suchen!
  (RAW.Data,pick = 3)
  # x - Achse => die echten Werte auf einer Log-Scala
  # der erste Bogen ist der VaR, der vierte Bogen ist der ES
  
}
