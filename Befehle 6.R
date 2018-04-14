# OEKO 2
function(){
  # Create sample data
  set.seed(1)
  x <- arima.sim(n = 100, list(ar=c(0.9),ma = c(0,0,0,0,0.5)), sd = 1)
  x <- as.timeSeries(x)
  ### Diskrete Fourier Transformation
  function(){
    DFT <- function(x_k){
      len <- length(x_k)
      DFT_b<-0:(len/2)
      DFT_c<-0:(len/2)
      for (k in 0:(len/2))
      {
        cexp <- complex(arg=-(1:len)*2*pi*k/len); cexp
        # Complex conjugate: this will be used for computing the IDFT
        cexpt <- complex(arg=(1:len)*2*pi*k/len); cexpt
        four<-sum(cexp*x_k*sqrt(1/(2*pi*len)))
        # Complex conjugate: this will be used for computing the IDFT
        fourc<-sum(cexpt*x_k*sqrt(1/(2*pi*len)))
        DFT_b[k+1]<-four
        # Complex conjugate: this will be used for computing the IDFT
        DFT_c[k+1]<-fourc
      }
      par(mfrow=c(3,1))
      ts.plot(x_k,main="Data: SAR(1)")
      acf(x_k)
      plot(abs(DFT_b),type="l",axes=F,col="blue",main="DFT")
      axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
                                      "4pi/6","5pi/6","pi"))
      axis(2)
      box()
      return(abs(DFT_b))
    }
    # komisch aussehendes E(w_k) bedeuted DFT
    # x_k = Zahlenreihe
    # Verschiebung der Sinuskurve kann ausgelesen werden
    # Daten können von DFT auf Beobachtungen zurückgeführt werden
  }
  
  ### Periodogram-Funktion
  # Periodigramm = I_TX(x_t) = (|DFT|)^2  (Betrag und quadrat von DFT)
  # plot_T = T oder F
  # Analytisch = I_tx(w_k)
  # peak dort wo die eingegebene Reihe saisonal ist.
  # Periodogram
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
    if (plot_T)
    {
      par(mfrow=c(2,1))
      plot(per,type="l",axes=FALSE,xlab="Frequency",ylab="Periodogram",
           main="Periodogram")
      axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
                                      "4pi/6","5pi/6","pi"))
      axis(2)
      box()
      plot(log(per),type="l",axes=FALSE,xlab="Frequency",ylab="Log-periodogram",
           main="Log-periodogram")
      axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
                                      "4pi/6","5pi/6","pi"))
      axis(2)
      box()
    }
    return(list(DFT=DFT,per=per))
  }
  # length(x)/2+1 = length(per$b)
    periodogram <- per(x,T)
  
  # Anzahl Schritte pro Saison = K
  K<-ceiling(length(periodogram$per)-1) # bei ungerader Zahl 1 abgezogen und aufgerundet (das selbe wie einmal abrunden)
  K<-as.integer(len/2)
  
  ### Gamma1 (amplitude truncated ideal Filter)
  function(){
    K<-as.integer(len/2)
    ### A Lowpass Filter
    cutoff<-pi/6
    Gamma<-((0:K)*pi/K)<cutoff
    
    plot(Gamma,type="l",axes=F,xlab="Frequency",main="Target Gamma",ylab="")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
  }
  ### Gamma2 (amplitude punctual filter)
  function(){
    # K<-ceiling(length(periodogram$per)-1)
    K<-as.integer(len/2)
    # Building your own Filter
    Gamma<-rep(1,K+1)
    Gamma[1+(K/6)+(-1:1)]<-0 #Diese Zeile abändern!!!! diese Zeile kann mehrmals auf das Gamma angewendet werden
    plot(Gamma,type="l",axes=F,xlab="Frequency",main="Target Gamma",ylab="")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
  }
  ### classical seasonal adjustment Flter B^12 (bad, better make your own, looks like a sinus)
  function(){
    b12<--1; omega_k<-(0:K)*pi/K
    trffkt<-1+b12*complex(arg=12*omega_k)
    ts.plot(abs(trffkt))
  }
  ### Bandpass-Filter
  function(){
    # bei 0 muss die zahl extrem gross gewählt werden !
    cutoff_1 <- pi/12
    cutoff_2 <- pi/6
    len <- 1000
    K <- as.integer(len/2)
    
    Gamma<-rep(0,K+1)
    Gamma<-((0:K)*pi/K)<cutoff_2
    Gamma[((0:K)*pi/K)<cutoff_1]<-FALSE
    plot(Gamma,type="l",axes=F,xlab="Frequency",main="Target Gamma",ylab="")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
  }
  
  ### realer Filter1 symmetric NUR WENN SYMETRIC GEFRAGT (Finite truncated filter, symmetric lowpass Filter)
  function(){
    
    K<-as.integer(len/2)
    # grösse der Schwankungen
    ord<-800 # K = -800...0...800 => ord = 800  !AND! K muss nicht 800 sein (k wird gleich wie immer berechnet)
    gammak<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))
    # Compute finite sum
    # len1: Number of frequency ordinates (resolution of discrete grid in [-pi,pi])
    len1<-K
    Gamma_hat<-0:len1
    for (k in 0:len1)#k<-0
    {
      omegak<-k*pi/len1
      Gamma_hat[k+1]<-gammak%*%(cos(omegak*0:ord)*c(1,rep(2,ord)))
      # alternativ mit gammak gespiegelt
      Gamma_hat[k+1]<-c(gammak[(ord+1):2],gammak)%*%c(cos(omegak*ord:1),cos(omegak*0:ord))
    }
    
    plot(abs(Gamma_hat),type="l",axes=F,main="Gamma (blue) and finite approximation (black)",xlab="",ylab="")
    lines(Gamma,col="blue")
    axis(1, at=c(0,1+1:6*len1/6),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
    
    ### Exakte Amplitude
    Gamma_hat
    ### Funktion die angewendet wird
    gammak
    
  }
  ### SYMETRISCHER-Bandpassfilter Approximation nicht symetrisch => mit dfa()
  function(){
    cutoff_1 <- pi/12 # kleinerer wert für passband
    cutoff_2 <- pi/6 # grösserer wert für passband
    ord <- 900 # K = -900...0...900 => ord = 800  !AND! K muss nicht 800 sein (k wird gleich wie immer berechnet)
    
    gamma_1<-c(cutoff_1/pi,(1/pi)*sin(cutoff_1*1:ord)/(1:ord))
    gamma_2<-c(cutoff_2/pi,(1/pi)*sin(cutoff_2*1:ord)/(1:ord))
    # Compute finite sum (approximation)
    Gamma_hat_pb<-0:ord
    for (k in 0:ord)#k<-0
    {
      omegak<-k*pi/ord
      Gamma_hat_pb[k+1]<-(gamma_2-gamma_1)%*%(cos(omegak*0:ord)*c(1,rep(2,ord))) # skript seite 43 
    }
    plot(Gamma_hat_pb,type="l",axes=F,
         main="Passband Gamma (blue) and finite approximation (black)",xlab="",
         ylab="")
    axis(1, at=c(0,1+1:6*ord/6),labels=c("0","pi/6","2pi/6","3pi/6",
                                         "4pi/6","5pi/6","pi"))
    axis(2)
    box()
    
    Gamma_hat  # approximierter filter
    gammak <- gamma_2
    
  }
  ### realer Filter2 (punctual filter) machen, Berechnen der Parameter des Filters (dfa_ms)
  function(){
    ## Diese Methode kann NUR Meansquare berechnen
    # L = Länge des Filters bzw. Filtergewichte dfa$b (ehöhen verbessert die qualität)
    #     2π/cutoff = Dauer = guter Wert für L (weiteste Frequenz link wählen)
    # periodigramm = per(x)
    # Lag = Schritte bis in die Zukunft...?
    # Gamma = Angewendeter Filter
  } # Beschreibung
  dfa_ms<-function(L,periodogram,Lag,Gamma){
    
    K<-length(periodogram)-1
    X<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)*sqrt(periodogram)
    X_y<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)
    for (l in 2:L)          #l<-L<-21
    {
      X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                    1.i*sin((l-1-Lag)*pi*(0:(K))/(K)))*sqrt(periodogram))
      X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                        1.i*sin((l-1-Lag)*pi*(0:(K))/(K))))
    }
    xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
    # MA-Filtercoefficients
    b<-as.vector(solve(xtx)%*%(t(Re(X_y))%*%(Gamma*periodogram)))
    # Transferfunction
    trffkt<-1:(K+1)
    trffkt[1]<-sum(b)
    for (k in 1:(K))#k<-1
    {
      trffkt[k+1]<-(b%*%exp(1.i*k*(0:(length(b)-1))*pi/(K)))
    }
    return(list(b=b,trffkt=trffkt))
  }
    # die Länge vom Periodigramm und Gamma müssen gleich sein
    dfa_ms(20,periodogram$per,0,Gamma[1:length(periodogram$per)])
  
  ### Amplitude berechnen
  amp<-abs(dfa$trffkt)
  ### Time-Shift berechnen
  K<-ceiling(length(periodogram$per)-1)
  shift<-Arg(dfa$trffkt)[2:length(dfa$trffkt)]/(pi*(1:(K))/(K))
  shift[1]<-sum((0:(L-1))*dfa$b)/sum(dfa$b) # Verzögerung in Frequenz Null
  ### Plot der Amplitude und des Time-Shiftes
  function(){
    par(mfrow=c(1,2))
    plot(amp,type="l",axes=F,xlab="Frequency",main="Amplitude",ylab="")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
    plot(shift,type="l",axes=F,xlab="Frequency",main="Shift",ylab="")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
    par(mfrow=c(1,1))
  }
  
  ### Anwenden1 des Filters auf die Daten und Plotten der Auswirkung auf das Periodigramm
  function(){
    xx <- x//x_i
    (xx <- read.csv("Studium/5.Semester/OEKO 2 - Mail/Daten/bp_intraday.csv"))
    
    # ord => wenn k -900...0...900 => ord= 900
    ord<-10 # muss kleiner als die länge von xx gewählt werden
    
    ### Anwenden des Filters auf die Daten
    par(mfrow=c(2,1))
    gammak<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))
    yy<-xx
    # (ord+1):(len-ord) = 901:2000 abzudeckender Bereich
    len <- length(yy)
    for (i in (ord+1):(len-ord)) { c[i]<-gammak%*%xx[i:(i-ord)]+gammak[2:ord]%*%xx[(i+1):(i+ord-1)] }
    
    mplot<-cbind(xx,yy)[(ord+1):(len-ord),]
    ts.plot(mplot,lty=1:2)
    len1<-dim(mplot)[1]
    
    ### Plotten der Auswirkung auf das Periodigramm
    plot((per(mplot[,2],F)$per),type="l",axes=F,col="red",main="Output (red) vs. input (blue) vs. convolution (green)")
    lines((per(mplot[,1],F)$per),col="blue")
    axis(1,at=1+0:6*(len1+1)/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
    par(mfrow=c(1,1))
    
    xf_i//y <- yy
  }
  ### Anwenden2 des Filters mit dem  real-time DFA Filter
  function(){
    xx <- x_i//x
    
    par(mfrow=c(2,1))
    len <- length(xx)
    # die ersten L (anzahl Filtergewichte) Daten werden weggelassen, die Zahlenreihe wird verkürzt.
    yy<-xx
    for (i in L:length(xx)){ yy[i]<-dfa$b%*%xx[i:(i-L+1)] }
    ts.plot(cbind(yy,xx),col=c("black","grey")) #Plotten des Resultats
    lines(y_i,col="red") # Vorsicht!!!!!!!!!!!!
    
    # Periodigramm Vergleich
    # Ist die Saisonalität am gewünschten Ort zurück gegangen?
    # Falls nicht kann der Filter auch ein zweites mal darüber angewendet werden
    # oder das Periodigramm kan an der gewünschten Stelle proportional VERGRÖSSERT werden. (NICHT 0 SETZTEN!)
    # L <- length(dfa$b)
    plot((per(yy[L:len],F)$per),type="l",axes=F,col="red",main="Output (red) vs. input (blue)")
    lines((per(xx[L:len],F)$per),col="blue")
    axis(1,at=1+0:6*(len-L+1)/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
    par(mfrow=c(1,1))
    
    y//xf_i <- yy
  }
  ### Anwenden3 auf OUT-OF-SAMPLE
  function(){
    # Trick: wir kennen ja die vergangenen Beobachtungen von x. Wir setzen sie in einem verlängerten vektor z ein
    # x_i <- insample
    # x_o <- outsample
    # L <- übergeben von oben
    z_o<-c(x_i[(length(x_i)-(L-1)):length(x_i)],x_o)
    xf_o<-x_o
    for (i in 1:length(x_o))
    {
      xf_o[i]<-dfa$b%*%z_o[L+i:(i-L+1)]
    }
    
    ts.plot(cbind(xf_o,x_o),lty=1:2)
    
    per_x_o<-per(x_o,F)$per
    per_xf_o<-per(xf_o,F)$per
    # Jetzt kann das Periodogram ?ber den ganzen Zeitbereich angewendet werden
    plot(per_xf_o,type="l",axes=F,col="red",main="Output (red) vs. input (blue)",ylim=c(0,max(per_x_o)))
    lines(per_x_o,col="blue")
    axis(1,at=1+0:6*(length(x_o))/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
    
    
    # Out-of-sample MSE im Zeitbereich berechnen
    mean((y_o-xf_o)^2)
    # alles um die 25% sollte ok sein
    
    # Zum Vergleich: hier noch einmal die in-sample Werte
    # Time-domain MSE OUT OF SAMPLE
    mean(((y_i-xf_i)[24:120])^2)
    # Frequency-domain (the real-time filter with transfer
    # function trffkt minimizes this expression for all imaginable filters of length L)
    (2*pi/length(Gamma))*abs(Gamma-trffkt)^2%*%weight_func
  }
  
  ### MBA filter (model-based)
  function(){
    # Hypothetical sample length
    len<-120
    # Symmetric lowpass target
    cutoff<-pi/12
    # Frequency resolution
    K<-len*10
    # Order of approximation : 10, 100, 1000, 10000
    ord<-len
    # Compute coefficients gamma
    
    
    gamma_k<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))
    # sum(gamma_k)+sum(gamma_k[2:ord])
    omega_k<-(0:K)*pi/K
    
    # AR(1)-coefficient
    a1 <- 0.3
    a2 <- 0.6
    b1 <- 0.9
    # Model-based filter
    gamma_0<-gamma_k%*%(a1^(0:ord)) # AR(1)
    gamma_0<-gamma_k%*%(a1^(0:ord)+a2^(0:ord)) # AR(2)
    gamma_0<-gamma_k%*%(a1^(0:ord)/b1^(0:ord)) # ARMA(1,1)
    gamma_0<-gamma_k%*%((a1^(0:ord)+a2^(0:ord))/b1^(0:ord)) # ARMA(2,1)
    
    gamma_mba_rt<-c(gamma_0,gamma_k[2:(ord)])
    ts.plot(gamma_mba_rt)
    
    trffkt_mba<-rep(NA,K+1)
    for (i in 0:K)
    {
      trffkt_mba[i+1]<-gamma_mba_rt%*%exp(1.i*omega_k[i+1]*(0:(ord-1)))
    }
    
    amp<-abs(trffkt_mba)
    shift<-Arg(trffkt_mba)[2:length(trffkt_mba)]/(pi*(1:(K))/(K))
    par(mfrow=c(2,1))
    ts.plot(amp,main="Amplitude")
    ts.plot(shift,main="Shift")
    
    # gamma_mba_rt  = dfa$b
    # trffkt_mba    = dfa$trffkt
    
    
    #### dfa_analytics wird an MBA angepasst
    # a1, a2, b1 usw. die selben wie oben
    omega_k<-(0:K)*pi/K
    weight_func<-1/(abs(1-a1*exp(1.i*omega_k))^2*2*pi)
    weight_func<-abs(1+b1*exp(1.i*omega_k))^2/(abs(1-a1*exp(1.i*omega_k))^2*2*pi)
    # Estimate filter coefficients
    dfa_ar1<-dfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2)
    b_dfa<-dfa_ar1$b
    
    ### Generate white noise, use arima and In- Out of sample
    function(){
      len<-120
      set.seed(pi)
      #x_all<-rnorm(2*len)
      ar1<-0.9
      x_all<-arima.sim(n = 2*len, list(ar = ar1))
      x_in<-x_all[1:len]
      x_out<-x_all[len+1:len]
    }
    ### Compare Filters
    function(){
      # dfa muss vorher berechnet sein
      benchmark<-cbind(dfa$b,gamma_mba_rt)
      ts.plot(benchmark,lty=1:2)
    }
  }
  ### Convolution 
  function(){
    ## Convolution
    # Das Periodigramm der x Reihe kann mit dem Filter^2 multipliziert werden
    # = periodigramm von y (y = Reihe mit angewendetem Filter)
    
    # Theoretisch perfelt angewendeter Filter
    lines((per(mplot[,1],F)$per)*Gamma_hat[seq(from=1,to=length(Gamma_hat),length.out=length(per(mplot[,1],F)$per))]^2,
          col="green")
    lines((per(mplot[,2],F)$per)*Gamma_hat_2^2,col="green")
    # convolution Theorem (Approximations Fehler)
    conv <- (per(mplot[,1],F)$per)*
      Gamma_hat[seq(from=1,to=length(Gamma_hat),length.out=length(per(mplot[,1],F)$per))]^2-per(mplot[,2],F)$per
    lines(conv)
    per(conv,T)
  }
  ### Annulaized Sharpe Performance Mass
  function(){
    
    ### Anwenden des Filters (schon gemacht???)
    ## wird benötigt:
    # dfa$b
    # L <- 2π/cutoff = Dauer = L
    
    # Filter anwenden
    xf<-0*x
    for (i in L:length(x)) xf[i]<-dfa$b%*%x[i:(i-L+1)]
    
    xf <- xf0
    ###
    shift<-1
    perf_weight<-cumsum(xf[1:(length(xf)-shift)]*x[(1+shift):length(x)])/sqrt(var(xf[1:(length(xf)-shift)]))
    ts.plot(perf_weight)
    
    shift<-1
    perf_sign<-cumsum(sign(xf[1:(length(xf)-shift)])*x[(1+shift):length(x)])
    ts.plot(perf_sign)
    
    
    ### Annualized sharpe
    # in sample (länge angeben)
    span<-1:499
    # out of sample (länge angeben)
    span<-1:1000
    sharpe_long_only<-sqrt(96*250)*mean(x[span])/sqrt(var(x[span]))
    sharpe_weight<-sqrt(96*250)*mean(diff(perf_weight[span]))/sqrt(var(diff(perf_weight[span])))
    sharpe_sign<-sqrt(96*250)*mean(diff(perf_sign[span]))/sqrt(var(diff(perf_sign[span])))
    
    
    sharpe_long_only
    sharpe_weight
    sharpe_sign
  }
  
  ### IN and OUT of Sample unterteilen
  function(){
    len_i<-120 # IN-SAMPLE Span
    len_o<-dim(mplot)[1]-len_i # OUT-SAMPLE Span
    
    # define in- and out-of-sample data
    # Daten werden in in-sample und out-sample geteilt
    x_i<-mplot[1:len_i,1]
    y_i<-mplot[1:len_i,2]
    x_o<-mplot[(len_i+1):dim(mplot)[1],1]
    y_o<-mplot[(len_i+1):dim(mplot)[1],2]
    
    in_sample_data<-cbind(x_i,y_i)
    colnames(in_sample_data)<-c("data x","   (generally unobserved) target y")
    
    out_of_sample_data<-cbind(x_o,y_o)
    colnames(out_of_sample_data)<-c("data x","   (generally unobserved) target y")
    periodogram <- per(x_o,T)
    
    # xf_i = nach anwenden des Filters auf x_i
    # y_i = nach anwenden des Filters auf x ZUGESCHNITTEN auf x_i
    # y = nach anwenden des Filters auf x
  }
  
  ### MSE (=mean square errors),(Time-domain,Frequency-domain)
  function(){
    # Time-domain
    # y_i  = GANZER Filter wurde bereits angewendet
    # xf_i = IN-SAMPLE Filter wurde auf die Daten bereits angewendet
    mean(((y_i-xf_i)[L:length(xf_i)])^2)
    
    # Time-domain OUT OF SAMPLE
    mean((y_o-xf_o)^2)
    
    # Frequency-domain
    # (the real-time filter with transfer function trffkt minimizes this
    #  expression for all imaginable filters of length L)
    # Gamma  = Amplitude Theoretisch exact des IN-SAMPLE Filters
    # trffkt = Amplitude des IN-SAMPLE Filters
    (2*pi/length(Gamma))*abs(Gamma-trffkt)^2%*%per(x_i,F)$per
    
    
    # Out-of-sample MSE im Zeitbereich berechnen
    function(){
      # Time-domain
      mean(((y_i-xf_i)^2)[L:length(y_i)])
      # Frequency-domain (the real-time filter with transfer function
      #  trffkt minimizes this expression for all imaginable filters of length L)
      (2*pi/length(Gamma))*abs(Gamma-trffkt)^2%*%per(x_i,F)$per
    }
  }
  
  
  ### Erwiterte Parameter berechnugn (dfa_analytic)
  function(){
    # This function computes analytic DFA-solutions
    # L            - length of the MA-filter
    #                2π/cutoff = Dauer = guter Wert für L
    # weight_func  - periodogram
    # Gamma        - transferfunction of the symmetric filter (target)
    # cutoff       - Schrittgrösse pro Saison (π/...   ; π/12 = 12 Schritte)
    # Lag          - lag-parameter: Lag=0 implies real-time filtering, Lag=L/2 implies symmetric filter
    # Constraints
    # i1           - amplitude-restriktionen (TRUE/FALSE, Standart = F, TRUE meistens uninteressant)
    #                erhöht Time-shift und amplitude leicht von grossen Frequenzen
    #                  Amplitude[0] = 1 festgesetzt
    # i2           - timeshift-restriktionen (TRUE/FALSE, Standart = F, TRUE kann manchmal nützlich sein)
    #                reduziert Time-shift für grosse Frequenzen, erhöht die amplitude
    #                  Shift[0] = 0 festgesetzt
    # Customization
    # alles 0 => Mean-square (kleinste Quadrate)
    # lambda        - Time-shift verkleinern      smoothness      Erfahrungswert: 5-10    selten >30
    #                wird Lambda erhöht, verkleinert sich der timeshift und die Amplitude erhöht sich.
    #                besonders in grossen Frequenzen
    # eta          - Amplitude verkleinern        timeliness      Erfahrungswert: 1-1.5   nie >1.5
    #                wird eta erhöht, verkleinert sich die Amplitude und der Time-shift vergrössert sich
    #                Amplitude wird überall etwa gelich verändert, Time-Shift vorallem in grossen Frequenzen
    # The function returns the weights of the MA-Filter as well as its transferfunction
  } # Beschreibung
  dfa_analytic<-function(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2){
    K<-length(weight_func)-1
    # Define the amplitude weighting function weight_h (W(omega_k))
    omega_Gamma<-as.integer(cutoff*K/pi)
    if ((K-omega_Gamma+1)>0)
    {
      weight_h<-weight_func*(c(rep(1,omega_Gamma),(1:(K-omega_Gamma+1))^(eta)))
    } else
    {
      weight_h<-weight_func*rep(1,K+1)
    }
    # First order filter restriction: assigning a large weight to frequency zero
    if (i1)
      weight_h[1]<-max(1.e+10,weight_h[1])
    
    X<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)*sqrt(weight_h)
    X_y<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)
    if (i2)
    {
      # Second order restriction: time shift in frequency zero vanishes
      for (l in 2:(L-1))
      {
        X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
                      cos((L-1-Lag)*pi*(0:(K))/(K))+
                      sqrt(1+Gamma*lambda)*1.i*(sin((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
                                                  sin((L-1-Lag)*pi*(0:(K))/(K))))*sqrt(weight_h))
        X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
                          cos((L-1-Lag)*pi*(0:(K))/(K))+
                          1.i*(sin((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*sin((L-1-Lag)*pi*(0:(K))/(K)))))
      }
      xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
      # MA-Filterweights
      b<-as.vector(solve(xtx)%*%(t(Re(X_y))%*%(Gamma*weight_h)))
      # the last weight is a function of the previous ones through the second order restriction
      b<-c(b,-sum(b*(0:(length(b)-1)))/(length(b)))
    } else
    {
      for (l in 2:L)
      {
        X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                      sqrt(1+Gamma*lambda)*1.i*sin((l-1-Lag)*pi*(0:(K))/(K)))*sqrt(weight_h))
        X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                          1.i*sin((l-1-Lag)*pi*(0:(K))/(K))))
      }
      xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
      # MA-Filterweights
      b<-as.vector(solve(xtx)%*%(t(Re(X_y))%*%(Gamma*weight_h)))
    }
    # Transferfunction
    trffkt<-1:(K+1)
    trffkt[1]<-sum(b)
    for (k in 1:(K))#k<-1
    {
      trffkt[k+1]<-(b%*%exp(1.i*k*(0:(length(b)-1))*pi/(K)))
    }
    return(list(b=b,trffkt=trffkt))
  }
  
  
  
  
  ####### Ablauf zum Aufbau eines ARMA-Filters:
  
  # avec kann auf 0 gesetzt werden
  # bvec muss IMMER mit c(1,...) anfangen
  # len1 = anzahl Unterbrüche (Qualität)
  # omega = Ausgabe der transferfuntion an diesem Punkt. (ändert nichts an den Daten)
  arma.u <- function(avec,bvec,len1,omega=pi/12){
    
    par(mfrow=c(2,1))
    a<-c(1,rep(0,length(avec)))
    b<-rep(0,length(bvec))
    a[1:length(avec)+1]<-avec
    b[0:length(bvec)]<-bvec
    
    bv <- b[1]
    if(length(b)>1) for(i in 2:length(b)){ bv <- bv+b[i]*complex(arg=-(0:len1)*(i-1)*pi/len1) }
    av <- a[1]
    if(length(a)>1) for(i in 2:length(a)){ av <- av-a[i]*complex(arg=-(0:len1)*(i-1)*pi/len1) }
    trffkt <- bv/av
    
    plot(abs(trffkt),type="l",axes=FALSE)
    title("amplitude")
    axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
    
    omega_grid<-((0:len1)*pi/len1)
    plot(-Arg(trffkt)/omega_grid,type="l",axes=FALSE)
    title("time-shift")
    axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
    par(mfrow=c(1,1))
    
    # Amplitudeneffekt in Omega
    amp_eff_ma1<-abs(trffkt)[(len1)/(pi/omega)+1]
    amp_eff_ar1<-abs(trffkt)[(len1)/(pi/omega)+1]
    
    # Shifteffekt
    shift_eff_ma1<--(Arg(trffkt)/omega_grid)[(len1)/(pi/omega)+1]
    shift_eff_ar1<--(Arg(trffkt)/omega_grid)[(len1)/(pi/omega)+1]
    
    # Output
    shifted_scaled_series<-amp_eff_ma1*amp_eff_ar1*cos(((1:len1)-shift_eff_ma1-shift_eff_ar1)*omega)
    
    # Shift übergeben
    Shift <- -Arg(trffkt)/omega_grid
    
    back <- list(shifted_scaled_series,Shift,trffkt,amp_eff_ma1,shift_eff_ma1)
    names(back) <- c("KorrigierterOutput","shift","amplitude", "Amplitudeneffekt_in_Omega", "Shifteffekt_in_Omega")
    return(back)
  }
  arma.u(avec = 0,bvec = c(1,0.9),len1 = 20)
  
  ### Anwenden von ARMA-Filter auf Daten
  function(){
    # Der Effektive Filter ist die ARMA Gleichung
    # dh abbildern der Gleichung  y_t = a_1*y_t-1 + a_2*y_t-2 + b_0*x_t + b_1*x_t-1
    x<-rnorm(len)
    y<-x
    
    # for Schlaufe startwert anpassen!
    # MA(3)
    startwert <- 13 # eins grösser als der grösste lookback
    y <- x # WICHTIG!!! MUSS FEST GELEGT SEIN
    b12 <- 1
    for (i in startwert:length(x)) y[i]<-b12*x[i-12]+x[i]
    
    # AR(12)
    function(){
      for (i in 13:len) y[i]<-y[i-12]+x[i]
    }
    # ARMA(1,1)
    function(){
      a1 <- 0.9: b1 <- 0.5
      for (i in 2:len) y[i]<-a1*y[i-1]+x[i]+b1*x[i-1]
    }
    
    # Zum vergleich: die Reihen vor/nach Filterung
    per(x,T)
    per(y,T)
  }
  
  
  ### RainbowPlot
  function(){
    par(mfrow=c(2,2))
    amp<-abs(trffkt)
    shift<-Arg(trffkt)/omega_k
    for (i in 2:2)
    {
      ymin<-min(amp[,i,],na.rm=T)
      ymax<-max(amp[,i,],na.rm=T)
      plot(amp[,i,1],type="l",main=paste("Amplitude functions, a1 = ",a_vec[i],sep=""),
           axes=F,xlab="Frequency",ylab="Amplitude",col=colo[1],ylim=c(ymin,ymax))
      mtext("Lag=0", side = 3, line = -1,at=len/4,col=colo[1])
      for (j in 2:(L/2+2))
      {
        lines(amp[,i,j],col=colo[j])
        mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=len/4,col=colo[j])
      }
      axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                           "4pi/6","5pi/6","pi"))
      axis(2)
      box()
      ymin<-min(shift[,i,],na.rm=T)
      ymax<-max(shift[,i,],na.rm=T)
      plot(shift[,i,1],type="l",main=paste("Time-Shifts, a1 = ",a_vec[i],sep=""),
           axes=F,xlab="Frequency",ylab="Shift",col=colo[1],ylim=c(ymin,ymax))
      mtext("Lag=0", side = 3, line = -1,at=len/4,col=colo[1])
      for (j in 2:(L/2+2))
      {
        lines(shift[,i,j],col=colo[j])
        mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=len/4,col=colo[j])
      }
      axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                           "4pi/6","5pi/6","pi"))
      axis(2)
      box()
      ymin<-min(b[,i,],na.rm=T)
      ymax<-max(b[,i,],na.rm=T)
      plot(b[,i,1],col=colo[1],ylim=c(ymin,ymax),main=paste("Filter coefficients"),
           ylab="Output",xlab="lag",axes=F,typ="l")
      for (j in 2:(L/2+2))
      {
        lines(b[,i,j],col=colo[j],type="l")
        mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=L/2,col=colo[j])
      }
      axis(1,at=1:L,labels=-1+1:L)
      axis(2)
      box()
      
      ymin<-min(yhat_Lag[,i,],na.rm=T)
      ymax<-max(yhat_Lag[,i,],na.rm=T)
      ts.plot(yhat_Lag[,i,1],col=colo[1],ylim=c(ymin,ymax),
              main=paste("Output series"),ylab="Output")
      for (j in 2:(L/2+2))
      {
        lines(yhat_Lag[,i,j],col=colo[j])
        mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=len/4,col=colo[j])
      }
      
    }
    
  }
  ### Tentacle Plot
  function(){
    vintage<-array(dim=c(len,3,len))
    dim(vintage)
    # For each of the three AR(1)-processes We compute the vintage series
    for (i in 1:3)
    {
      for (j in L:len)#j<-L
      {
        vintage[(j-as.integer(L/2)):j,i,j]<-yhat_Lag[j,i,(as.integer(L/2)+1):1]
        vintage[1:(j-as.integer(L/2)-1),i,j]<-
          yhat_Lag[(as.integer(L/2)+1):(j-1),i,as.integer(L/2)+1]
      }
    }
    # We select the third DGP with a1=-0.9
    i<-3
    vintage_triangle<-vintage[,i,]
    dimnames(vintage_triangle)[[2]]<-paste("Publ. ",1:len,sep="")
    dimnames(vintage_triangle)[[1]]<-paste("Target ",1:len,sep="")
    
    
    
    
    colo<-rainbow(len)
    # file = paste("z_vintages.pdf", sep = "")
    # pdf(# file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
    par(mfrow=c(3,1))
    for (i in 1:3)
    {
      ymin<-min(vintage[,i,],na.rm=T)
      ymax<-max(vintage[,i,],na.rm=T)
      ts.plot(vintage[,i,L],col=colo[1],ylim=c(ymin,ymax),
              main=paste("Tentacle plot: vintages and full revision sequence,
                       a1 = ",a_vec[i],sep=""),ylab="Vintages")
      for (j in (L+1):len)
      {
        lines(vintage[,i,j],col=colo[j])
      }
      lines(vintage[,i,len],col="red",lwd=2)
    }
    
    # Weitere Beispiel
    function(){
      par(mfrow=c(2,1))
      i<-2
      ymin<-min(vintage[,i,],na.rm=T)
      ymax<-max(vintage[,i,],na.rm=T)
      ts.plot(vintage[,i,L],col=colo[1],ylim=c(ymin,ymax),
              main="Vintages: full revision sequence and final release (black)",ylab="Vintages")
      for (j in (L+1):len)
      {
        lines(vintage[,i,j],col=colo[j])
      }
      lines(vintage[,i,len],col="black",lwd=2)
      i<-2
      ymin<-min(vintage[,i,],na.rm=T)
      ymax<-max(vintage[,i,],na.rm=T)
      ts.plot(vintage[,i,L],col=colo[1],ylim=c(ymin,ymax),
              main="Vintages: full revision sequence and real-time initial release (black)",
              ylab="Vintages")
      for (j in (L+1):len)
      {
        lines(vintage[,i,j],col=colo[j])
      }
      lines(yhat_Lag[,i,1],col="black",lty=1)
    }
    
    
  }
  
}

# STDM
function(){
  # von Herr Ruckstuhl
  # Um wie viel ist die Berechnung besser, wenn man welches Quantil wählt
  # Der Lift chart zeigt an um wie viel mal es mit dieser Whal besser ist als eine Zufalls wahl
  plot.liftchart <- function(true, predScore, quant=(1:10)/10, add=FALSE, col="blue", cum=FALSE, ylim=NULL){
    ## true = 0 oder 1 (hohe Score --> 1)
    
    ## true sortieren nach den Scorewerten:
    trueSort <- true[order(predScore, decreasing=TRUE)]
    
    ## Indices der entsprechenden Quantile suchen
    QIndex <- round(quantile(1:length(trueSort), quant))
    
    ## Fehlerraten an den Quantilpunkten bestimmen
    rate <- cRate <- rep(NA, length(QIndex))
    for (i in 1:length(QIndex)){
      rate[i] <- mean(trueSort[1:QIndex[i]])
      cRate[i] <- sum(trueSort[1:QIndex[i]])
    }
    ##
    rRate <- rev(rate)[1]
    if(add){
      if(cum) lines(quant, cRate/rev(cRate)[1], col=col, type="b")
      else lines(quant, rate/rRate, col=col, type="b")
    }
    else{
      if(cum){
        plot(quant, cRate/rev(cRate)[1], lwd=2, col=col, type="b",
             xlim=c(0,1), ylim=c(0,1), las=1,
             xlab="Quantile", ylab="Kumulativer Lift")
        abline(a=0,b=1, lty=2)
      }
      else{
        plot(quant, rate/rRate, lwd=2, col=col, type="b", las=1,
             xlab="Quantile", ylab="Lift", ylim=ylim)
        abline(h=1, lty=2)
      }
    }
  }
  plot.liftchart(true = as.numeric(oil[,1]),predScore = as.numeric(oilT.knnA))
  
  ### logistischer Klassifizierer
  predScore <- predict(Ban.mn, type="probs")
  ### Klassifikations Bäume und Random Forest
  predScore <- predict(Ban.rpG, newdata=Ban, type="prob")[,"TRUE"]
  
  
  # Skript beispiel S22
  setwd("~")
  dat <- read.table("Studium/5.Semester/STDM - server rkst/Daten/Protein.dat", header=TRUE, row.names=1)
  Ban <- read.table("Studium/5.Semester/STDM - server rkst/Daten/Bankrott.dat", header=TRUE)
  
  # Hauptkomponentenanalyse (PCA)
  # 1. Arbeitsblatt PCA
  function(){
    
    library(cluster)
    library("MASS")
    library(vegan)
    
    # Allgemein
    # Hauptkomponenten n = PCn
    
    # Standartisieren
    # varianzen standartisierung
    ColVarDat <- apply(dat,2,var)
    ColMeanDat <- colMeans(dat)
    PrS <- scale(x = dat, center=ColMeanDat, scale=ColVarDat)
    PrS <- scale(dat, center=mean(dat), scale=sd(dat))
    # MAD standartisierung
    PrS <- scale(dat, center=mean(dat), scale=mad(dat))
    # DIFF standartisierung?
    scale(dat, center=median(dat), scale=diff(range(dat)))
    # Bsp. Robust Skaliert
    # Wenn es viele Ausreisser gibt
    x.l <- apply(dat, MARGIN=2, FUN=median)
    x.s <- apply(dat, MARGIN=2, FUN=mad) ## FUN = ausgewählte Standartisierung von oben
    PrS <- scale(dat, center=x.l, scale=x.s)
    # Wann müsssen Daten Standartisiert werden?
    #   Wenn die Grundeinheiten nicht überein stimmen Bsp. (Meter/Kilogramm)
    #   Wenn die Varianzen zu unterscheidlich sind diag(var(dat))  ## boxplot(dat)
    #   Wenn die Varianz keine Wichtige Infomrmation ist die beibehalten werden muss
    
    
    # Berechnen der Matrizen
    Sx <- prcomp(PrS) # WICHTIG "scale=F"! (standartmässig auf FALSE)
    # Berechenn der Hauptkomponenten von Sx (Covarianz-Matrix)
    # wird verwendet, wenn keine Skalierung durchgeführt werden muss oder schon geschehen
    Rx <- prcomp(PrS, scale=T)
    # Berechenn der Hauptkomponenten von Rx (Korrelations-Matrix)
    # benutzen wenn noch sklaiert werden muss
    
    # Hauptkomponentenanalyse
    PrS.pcR <- prcomp(PrS)
    P.pcR.p <- predict(PrS.pcR) # Berechnen der Hauptcomponentenscores
    library(MASS)
    # par(mar=c(5,4,1,1), mgp=c(2, 0.8, 0), las=1)
    ausschnitt <- 5
    eqscplot(c(-ausschnitt,ausschnitt), c(-ausschnitt, ausschnitt), type="n",
             xlab="Erste Hauptkomponente", ylab="Zweite Hauptkomponente")
    text(P.pcR.p[,1], P.pcR.p[,2], labels=rownames(P.pcR.p), cex=0.7)
    # Beschreibung:
    #   Cluster: Anhäufungen von daten-Punkten
    #   Maveriks: Ausreisser
    
    
    ### Anzahl benötigter Haupkomponenten auslesen das selbe wie VarianceIndebtednessAnalysis
    screeplot(PrS.pcR, type="l", cex=0.7, las=0)
    summary(PrS.pcR)
    ### Varianz verschuldungs Analyse
    # Es werden nur die Hauptkomponenten betrachted, welche genügend Einfluss haben  ## plot(Rx)
    # 80% der Varianz müssen erklärt sein mit den gewählten Hauptkomponenten
    # Welche Attribute tragen am stärksten zu der gesamten Varinaz bei?
    VarianceIndebtednessAnalysis <- function (dat.prcomp) {
      names(dat.prcomp$sdev) <- paste("PC", 1:length(dat.prcomp$sdev), sep="")
      # damit PC richtig angeschrieben sind
      par(mfrow=c(1,3))
      screeplot(dat.prcomp, main="", las=2)
      abline(h=1, lwd=2, col="gray")
      title("Varianzen")
      cat("Auslesen der Elbogens. Ab wo sind die Varianzen zu klein.
          Atribute mit der Grössten Streuung (Varianz) wird gesucht
          Faustregel 3 (Kaiser-kriterium): (nur bei Rx) verwende alle Komponenten mit  Eigenwert>1")
      cat("                                                                                    ")
      h.cvar <- cumsum(dat.prcomp$sdev^2); h.cvar <- h.cvar/rev(h.cvar)[1]
      plot(h.cvar, type="b", ylab="Cumulative Proportion of Variance",
           xlab="Number of Components")
      title("Cumsum Varianzen")
      abline(h=0.8, lwd=2, col="gray")
      plot(dat.prcomp$sdev^2,type = "b")
      title("Varianzen")
      par(mfrow=c(1,1))
      return(summary(dat.prcomp))
    }
    (dat.prcomp <- VarianceIndebtednessAnalysis(PrS.pcR))
    # 80-20-Regel: 70-80% der Varianz müssen erklärt werden von den gewählten Variablen
    
    # Interpretieren von Hauptkomponenten
    # Skript beispiel S27
    structure(dat.prcomp$rotation, class="loadings")
    # leere stellen = Werte sehr nahe bei Null => Variable hat keinen Einfluss
    # Analysieren:
    #   Nuller?
    #   Vorzeichen  =>  Negative Vorzeichen = sehr gute Qualität veredelt
    #                   Positive Vorzeichen = schlechtere Qualität Unveredelt
    #   Spalten     =>  Hauptkomponenten n = PCn
    #                   1 Spalte: Veredelung
    #                   2 Spalte: Teuer oder Billig
    
    biplot(PrS.pcR)
    # durch welche Faktoren werden die Hauptkomponenten hauptsächlich bestimmt?
    # Horizintale Pfeile => hoher Einfluss auf die 1 Hauptkomponente
    # Vertikale Pfeile => hoher Einfluss auf die 2 Hauptkomponente
    
    
    # im Fall der unsalierten Variablen
    dat.cov <- cov(data.matrix(dat))
    eigen(dat.cov)
  }
  
  # Multidimensionale Skalierung
  # Distanzen und MDS
  function(){
    # install.packages("cluster")
    library("MASS")
    library(cluster)
    library(vegan)
    
    ## euklidische Distanzfunktion
    dist() # Distance Matrix Computation
    daisy() # Dissimilarity Matrix Calculation
    ## verschiedene Metriken
    dist(dat, method= "manhattan") # method = c("euclidean","canberra","binary","minkowski","manhattan","maximum")
    daisy(dat, metric = "gower") # metric = c("euclidean", "manhattan", "gower")
    ## dist von Zeitreihen
    dat.c <- cor(dat[,-(1:2)]) # Korrelation der Zeitreihen
    dat.dist <- sqrt(2*(1 - dat.c)) # Transformation
    
    ### Ähnlichkeit zu Distanzen umrechnen sqrt(s_rr + s_tt - 2*s_rt)
    D <- sqrt((matrix(diag(S), ncol=ncol(S), nrow=nrow(S)) +
                 matrix(diag(S), ncol=ncol(S), nrow=nrow(S), byrow=TRUE) - 2*S))
    
    ### Distanzen schnell auslesen:
    round(as.matrix(daisy(dat))[nach, von],3)
    round(as.matrix(daisy(dat))[c(3, 5, 7), 1],3)
    
    ### Standartisieren
    # mit Standartabweichung , scale=sd() kann zu allem geändert werden
    scale.dat <- data.frame(x1=scale(dat[,1], center=mean(dat[,1]), scale=sd(dat[,1])),
                            x2=scale(dat[,2], center=mean(dat[,2]), scale=sd(dat[,2])))
    # mean absolute deviation standartisierung
    daisy(dat, stand=TRUE) # direkt bei daisy
    # wurzel + asin (arcsin) - transformation
    asin(sqrt(dat/100)) # /100!
    
    ### 3 verschiedene Arten die Distanzen dar zu stellen
    ## 1. Ploten der Daten (Dendrogramm):
    f.p <- function(obj, names=dimnames(obj)[[1]]){
      eqscplot(obj, type="n", xlab="Erste Hauptkoordinate",
               ylab="Zweite Hauptkoordinate", cex=0.7, las=1)
      text(obj, labels=names, cex=0.7)
    }
    f.p(cmdscale(dist(dat), eig=T, k=2)$points)
    ## 2. Falls f.p nicht gut funktioniert
    # vorsicht Resultat häufig willkürlich gespiegelt!!!!!!!
    HauptKOORDINATENAnalyse <- function (dat,ausschnitt=0,names=as.character(rownames(dat))) {
      jobT.mds <- cmdscale(dist(dat), k=2)
      library("MASS")
      par(mfrow=c(1,1))
      if(ausschnitt==0){
        eqscplot(jobT.mds, type="n", xlab="Erste Hauptkoordinate",
                 ylab="Zweite Hauptkoordinate", cex=0.7, las=1)
      }else{
        eqscplot(c(-ausschnitt,ausschnitt), c(-ausschnitt, ausschnitt), type="n", xlab="Erste Hauptkoordinate",
                 ylab="Zweite Hauptkoordinate", cex=0.7, las=1)
      }
      text(jobT.mds, labels=names, cex=0.7)
    }
    HauptKOORDINATENAnalyse(dat)
    ## 3. Distanzen in 2 Dimensionen darstellen
    ausschnitt <- 0.5
    eqscplot(c(-ausschnitt,ausschnitt), c(-ausschnitt, ausschnitt), type="n", xlab="Erste Hauptkoordinate",
             ylab="Zweite Hauptkoordinate", cex=0.7, las=1)
    text(cmddata$points, labels=attr(cmddata$points, "dimnames")[[1]], cex=0.7)
    
    
    ## Ordinale Variablen
    ordered(cut(dat[,1], c(0,25,50,75)))
    # right = TRUE oder FALSE liefern das selbe Ergebnis
    # Bedeutung = Rangordnung der Variablen wird eingehlaten, jedoch nicht der metrische Wert.
    #   Ordinalskalierte Merkmale sind kategorielle Merkmale, bei denen sich die verschiedenen
    #   Ausprägungen in eine sinnvolle Reihenfolge bringen lassen.
    dat.ord <- data.frame(ordered(cut(dat[,1], c(0,25,50,75))) , dat[,2])
    daisy(dat.ord)
    
    # Ordinale MDS (ordinale multidimensionale Skalierung (nach Kruskal and Shepard))
    # wenn die Distanz unklar oder ungenau ist wird mit Rängen gerechnet (Ordinal)
    # das erste Geröll wird gewählt! nicht das erste vor dem geröll, wie bei Varianz verschuldung
    Ordinale_MDS <- function(dat.dist){
      f.p <- function(obj, names=dimnames(obj)[[1]]){
        eqscplot(obj, type="n", xlab="Erste Hauptkoordinate",
                 ylab="Zweite Hauptkoordinate", cex=0.7, las=1)
        text(obj, labels=names, cex=0.7)
      }
      library("MASS")
      library(cluster)
      par(mfrow=c(1,2))
      Ord.mds <- isoMDS(dat.dist)
      f.p(Ord.mds$points)
      D.stress <- rep(NA,15)
      for(i in 1:length(D.stress)) D.stress[i] <- isoMDS(d=dat.dist, y=cmdscale(dat.dist, k=i),k=i, maxit=500 )$stress
      plot(1:length(D.stress), D.stress, type="b", las=1,ylab = "STRESS")
      title("Güte Plot")
      abline(h=0.05*100, col=2, lty=3)
      show("Stresswert der Daten")
      show(Ord.mds$stress)
      show("muss unter 20% liegen! 10% ist relativ gut")
      show("")
      show("Wie viele Komponenten benötigt es für eine gute Repräsentation im 2-Dimensionalen Raum?")
      show("Ist ein deutlicher Ellbogen sichtbar?")
      show("wie viele Komponenten benötigt es um unter oder in die Nähe eines STRESS von 5% zu kommen?")
      par(mfrow=c(1,1))
    }
    # um wie viel verbessert sich der Wert, wenn ich die x-te Variable hinzu nehme
    Ordinale_MDS(dist(dat))
    # Ordinale_MDS (stressplot mit funktion von oben)
    function(){
      library("MASS")
      
      # k <- Anzahl Hauptkomponenten (stress wird verkleinert, falls erhöht)
      CD.kas <- isoMDS(CD.dis, k=2)
      # Konvergiert bei einem final value von 8.470026
      # Der Wert 8.47 (=8.47%) zeigt, dass gem?ss Faustregel (?ber 20% schlecht,
      # unter 5% ausgezeichnet) die zweidimensionale  Darstellung geeignet ist.
      
      ## Plotten ordinale MDS
      par(mfrow=c(1,1))
      eqscplot(CD.kas$points, type="n", xlab="", ylab="", cex=0.7, las=1,
               sub="Ordinale MDS nach Kruskal and Shepard")
      text(CD.kas$points, labels=labels(CD.dis), cex=0.7, las=1)
      
      
      
      #Spanning Tree
      library(vegan)
      CD.L2mst <- spantree(CD.dis)
      plot(CD.L2mst, CD.kas, type="t", labels=labels(CD.dis), cex=0.7,
           xlab ="Erste Hauptkoordinate", ylab = "Zweite Hauptkoordinate")
      mtext(side=3, line=1, text="Ordinale MDS")
      ## Er zeigt auf, dass vor allem IND und EGY lokal verzerrt dargestellt werden.
      ## Gem?ss MST ist der n?chste Nachbar von IND EGY, und der andere n?chste Nachbar
      ## von EGY ist USA. Gem?ss ordinaler MDS m?sste der n?chste Nachbar von IND ZAI
      ## sein, und einer der n?chsten Nachbar von EGY m?sste ZAI oder BRA sein.
      
    }
    
    # "nonlinear mapping"-Methode von Sammon
    # Elbogenplot => Fehleranteil, mit k gewählten Dimensionen
    # (wie weit ist es noch von der Realität entfernt, wenn k dim gewählt werden)
    Sammmon_nonlinear_mapping <- function(dat.dist){
      f.p <- function(obj, names=dimnames(obj)[[1]]){
        eqscplot(obj, type="n", xlab="Erste Hauptkoordinate",
                 ylab="Zweite Hauptkoordinate", cex=0.7, las=1)
        text(obj, labels=names, cex=0.7)
      }
      library("MASS")
      library(cluster)
      D.sam <- sammon(dat.dist, k=2)
      f.p(D.sam$points)
      D.stress <- rep(NA,15)
      for(i in 1:length(D.stress)) D.stress[i] <-sammon(dat.dist, k=i, y=cmdscale(dat.dist, k=i))$stress
      plot(D.stress, type="b", las=1)
      abline(h=0.05, col="red", lty=2)
      title("Wie viele Komponenten braucht es?")
    }
    Sammmon_nonlinear_mapping(dist(dat))
    
    # metrische MDS
    (cmddata <- cmdscale(dist(dat), eig=T, k=2))
    # Gibt es eindeutig (nicht wegen numerischen Fehlern) mehrere negative Eigenwerte
    # => dann kann die Messung nicht geeignet im Euklidischen Raum dargestellt werden
    
    # Güte Analyse (Scree-Diagramm) von einer Dimension zur nächsten wird man um so viel genauer
    # Eigenwerte aufkummuliert
    BenoetigteDimensionen <- function (cmddata) {
      x <- cmddata$eig; x[x<0] <- 0
      plot(cumsum(x)/sum(x), type="b")
      abline(h=0.8,col="red",lty=2)
      title("80% müssen erreicht werden",cex=0.5)
      y <- cumsum(x)/sum(x)
      return(list(Eigenwerte=x,Cumsum_Norm_Eigenwerte=y))
    }
    BenoetigteDimensionen(cmddata)
    
    # Minimum spanning Tree
    function(){
      # D = Distanzen-Matrix
      library(vegan)
      tree <- spantree(dist(dat))
      
      # Either one of these methods
      dist.method <- cmdscale(dist(dat), k=2)   # Metrisches MDS und MST
      dist.method <- isoMDS(dist(dat))          # Ordinales MDS und MST
      dist.method <- sammon(dist(dat), k=2)     # Sammon's nonlinear mapping und MST
      plot(tree, dist.method, labels= rownames(dat))
      
      # Plotting all possibilities of the Spanningtrees
      SpannTree <- function (D,names=rownames(D),namegroesse=0.3, einzelplotten = 1:3) {
        if(is.null(names)){names <- attr(dist(D),"Labels")}
        par(mfrow=c(1,1), las=1)
        D.L2mst <- spantree(D)
        
        if(sum(einzelplotten == 1)>0||einzelplotten=="metrisch"){
          ## Metrisches MDS und MST:
          D.mds <- cmdscale(D, k=2)
          plot(D.L2mst, D.mds, type="t", labels=names, cex=namegroesse,
               xlab ="Erste Hauptkoordinate", ylab = "Zweite Hauptkoordinate")
          mtext(side=3, line=1, text="Metrische MDS")
        }
        
        if(sum(einzelplotten == 2)>0||einzelplotten=="ordinal"){
          ## Ordinales MDS und MST:
          D.kas <- isoMDS(D)
          plot(D.L2mst, D.kas, type="t", labels=names, cex=namegroesse,
               xlab ="Erste Komponente", ylab = "Zweite Komponente")
          mtext(side=3, line=1, text="Ordinale MDS")
        }
        
        if(sum(einzelplotten == 3)>0||einzelplotten=="sammon"){
          ## Sammon's nonlinear mapping und MST:
          D.sam <- sammon(D, k=2)
          plot(D.L2mst, D.sam, type="t", labels=names, cex=namegroesse,
               xlab ="Erste Komponente", ylab = "Zweite Komponente")
          mtext(side=3, line=1, text="Sammon's Nonlinear Mapping")
          ## Es gibt Verzerrungen: in allen drei Darstellungen ist z.B. die
          ## Position von "V" deutlich verzerrt.
        }
      }
      SpannTree(dist(dat),einzelplotten = "sammon")
    }
    
    # skalierte grössen plot
    eqscplot(dat[,1],dat[,2],ratio = 1,tol = 0.5)
    
  }
  
  # Cluster Analyse
  function(){
    ## Bsp Daten
    DPfad <- "Studium/5.Semester/STDM - server rkst/Daten/"
    load(paste(DPfad ,"CountriesDis.RDA",sep = "")) ## --> CD.dis
    CD.dis
    
    # Cluster Analyse
    CD.con <- hclust(CD.dis, method="single") # method = "single", "average", "complete", "ward.D2"
    plot(CD.con)
    
    ### Cluster in Gruppen unterteilen und Anzeigen Spanningtree (von Gianfranco)
    function(){
      # Ausgebend er Cluster zuweisung in einem Vector
      # k = anzahl ausgewählter Cluster
      groups <- cutree(CD.con,k = 3)
      # einzeichnen von Gruppierungen in die Clusterplotts
      rect.hclust(CD.con,k = 3,border = "red")
      
      ### Beschriften der Cluster in einen Plot (Farbe zuteilung in einem Spanntree nach cluster)
      library(vegan)
      sammon <- sammon(CD.dis, k=2)
      span <- spantree(CD.dis)
      plot(span,sammon$points, type="t", labels=span$labels,
           xlab ="Hauptkoordinatenachsen", col=groups)
    }
    
    # a quick look
    ClusterAnalyse <- function (dat.dis) {
      par(mfrow=c(1,4))
      CD.con <- hclust(dat.dis, method="single")
      plot(CD.con, labels=labels(dat.dis), main="Single-linkage")
      #
      CD.av <- hclust(dat.dis, method="average")
      plot(CD.av, labels=labels(dat.dis), main="Average-linkage")
      #
      CD.com <- hclust(dat.dis, method="complete")
      plot(CD.com, labels=labels(dat.dis), main="Complete-linkage")
      #
      CD.w <- hclust(dat.dis, method="ward.D2")
      plot(CD.w, labels=labels(dat.dis), main="Ward")
      par(mfrow=c(1,1))
      
      show("Single-linkage:")
      show("Eher ein Exotisches Verfahren, welchem nicht immer getraut werden kann")
      show("Average-linkage:")
      show("")
      show("Complete-linkage:")
      show("")
      show("Methode von Ward:")
      show("")
      show(paste("Wo werden die Abstände vertikal immer grösser? Sobald sie sehr gross sind kann vertikal geschnitten werden und die daten in anzhal Cluster geteilt werden"))
      par(mfrow=c(1,1))
    }
    ClusterAnalyse(CD.dis)
    
    
    CD.mMDS <- cmdscale(CD.dis, k=nrow(as.matrix(CD.dis))-1, eig=T)
    CD.mMDS$eig
    # NUMERISCHE FEHLER, WENN EIGENWERTE EXTREM KLEIN
    plot(CD.mMDS$eig, type="b"); abline(h=0, lwd=2)
    ## Drei Eigenwerte sind kleiner 0, und zwar mit einem Betrag, der nicht nur auf
    ## numerische Rundungsfehler zurückzuführen ist. NUMERISCHE FEHLER, WENN EIGENWERTE EXTREM KLEIN
    ## Folglich lassen sich diese Unähnlichkeiten nicht korrekt im euklidischen Raum
    ## darstellen.
    
    ### Cluster unterscheiden
    dist.method <- cmdscale(dist(CD.dis), k=2)   # Metrisches MDS und MST
    dist.method <- isoMDS(dist(CD.dis))          # Ordinales MDS und MST
    dist.method <- sammon(dist(CD.dis), k=2)     # Sammon's nonlinear mapping und MST
    
    eqscplot(dist.method$points, xlab="Erste Hauptkoordinate", ylab="Zweite Hauptkoordinate",
             main="Single Linkage Method", sub="Metrische MDS", cex=0.7, las=1)
    
    # lim = aussschnittgrösse
    ClusterPlot <- function (data.dist,anz.cluster=3, lim =2){
      #       if(method == "metrisch") dist.method <- cmdscale(dist(data), k=2)   # Metrisches MDS und MST
      #       if(method == "ordinal") dist.method <- isoMDS(dist(data))          # Ordinales MDS und MST
      #       if(method == "sammon") dist.method <- sammon(dist(data), k=2)     # Sammon's nonlinear mapping und MST
      data.method <- cmdscale(data.dist, k=2)
      
      syBsp.clS <- hclust(data.dist, method="single")
      syBsp.av <- hclust(data.dist, method="average")
      syBsp.clC <- hclust(data.dist, method="complete")
      syBsp.W <- hclust(data.dist, method="ward.D2")
      
      par(mfrow=c(2,2))
      eqscplot(data.method, xlab="Erste Hauptkoordinate", ylab="Zweite Hauptkoordinate",
               main="Single Linkage Method", sub="Metrische MDS", cex=0.7,
               las=1,xlim = c(-lim,lim), ylim = c(-lim,lim))
      #         for(i in 1:anz.cluster){
      #           ok <- cutree(syBsp.clS, k=anz.cluster)==i
      #           points(data.method[ok,], col=i, pch=16, cex=2)
      #         }
      #         text(data.method, labels=dimnames(syBsp)[[1]], cex=0.7,col=0 )
      ##
      
      eqscplot(data.method, xlab="Erste Hauptkoordinate", ylab="Zweite Hauptkoordinate",
               main="Average Linkage Method", sub="Metrische MDS",type="n", cex=0.7,
               las=1,xlim = c(-lim,lim), ylim = c(-lim,lim))
      for(i in 1:anz.cluster){
        ok <- cutree(syBsp.av, k=anz.cluster)==i
        points(data.method[ok,], col=i, pch=16, cex=2)
      }
      text(data.method, labels=dimnames(syBsp.av)[[1]], cex=0.7,col=0 )
      ##
      
      eqscplot(data.method, xlab="Erste Hauptkoordinate", ylab="Zweite Hauptkoordinate",
               main="Complete Linkage Method", sub="Metrische MDS", type="n", cex=0.7,
               las=1,xlim = c(-lim,lim), ylim = c(-lim,lim))
      for(i in 1:anz.cluster){
        ok <- cutree(syBsp.clC, k=anz.cluster)==i
        points(data.method[ok,], col=i, pch=16, cex=2)
      }
      text(data.method, labels=dimnames(syBsp.clC)[[1]], cex=0.7,col=0 )
      ##
      
      eqscplot(data.method, xlab="Erste Hauptkoordinate", ylab="Zweite Hauptkoordinate",
               main="Method von Ward", sub="Metrische MDS", type="n", cex=0.7,
               las=1,xlim = c(-lim,lim), ylim = c(-lim,lim))
      for(i in 1:anz.cluster){
        ok <- cutree(syBsp.W, k=anz.cluster)==i
        points(data.method[ok,], col=i, pch=16, cex=2)
      }
      text(data.method, labels=dimnames(syBsp.W)[[1]], cex=0.7,col=0 )
      par(mfrow=c(1,1))
    }
    ClusterPlot(CD.dis,anz.cluster = 3,lim=10)
    
    
    ### Heat Maps
    function(){
      dat <- read.table("Studium/5.Semester/STDM - server rkst/Daten/abst.dat", header=T, sep=",")
      dat.dist <- dist(dat)
      # dunkel = tiefe Werte
      # hell = hohe Werte
      # scale = skalieren der Zeilen oder Spalten
      heatmap(data.matrix(dat), scale="none") # scale = "none", "column", "row"
      
      # Sind Farbwürfel erkennbar? wenn ja welche Gruppen vereinigen sie?
      heatmap(as.matrix(dat.dist), symm=TRUE)
    }
    
    ### K-means (ist schlecht)
    function(){
      abst <- read.table("Studium/5.Semester/STDM - server rkst/Daten/abst.dat", header=T, sep=",")
      abst.dist <- dist(abst)
      abst.av <- hclust(abst.dist, method="average")
      
      ### K-means
      # funktionert auch auf sehr grosse Datenmengen, weil die
      # Daten nicht notwendigerweise im Hauptspeicher vorhanden sein müssen
      
      # Nur klare Cluster können unterschieden werden, weil die Mitte der Cluster bestimmt wird
      # und höufig nicht mit der wirklichen Mitte überein stimmt.
      
      ### Optimale Anzahl CLuster für K-Mean berechenn
      x.ss <- rep(NA,16)
      for(i in 2:16){
        abst.km <- kmeans(abst, centers=i, nstart=11)
        x.ss[i] <- sum(abst.km$withinss)
      }
      x.ss[1] <- sum(sapply(abst, FUN=var))*(nrow(abst)-1)
      plot(1:16, x.ss, type="b", ylim=range(c(0,x.ss)))
      abline(h=0,lty=2)
      # Knick bei 4 und 10
      # Grundsätzlich muss es eine MONOTON FALLENDE KURVE sein (Falls nicht, fand der
      # Algorithmus keine gute Lösung!)
      # 4 oder sogar 10 Clusters sind erforderlich
      
      ### K-means einzeln berechenn
      # centers = Anzahl Mittelpunkte
      (kmean <- kmeans(abst, centers=3, nstart=11))
      kmean$withinss
      ### Berechnen des K-mean
      abst.km0 <- kmeans(abst, data.matrix(abst)[c("TG","SG","ZH"),]) # K = 3 => drei Zeilen auswählen
      abst.km0$tot.withinss
      abst.km0$centers
      abst.km0$cluster
      
      (h.ct <- cutree(abst.av,3))   # Identifizieren der Cluster
      (h.ex <- rep(h.ct, ncol(abst))) # Clusterzugehörigkeit für jedes Element der
      (h.list <- list(h.ex, col(abst))) # Spaltenzugehörigkeit für jedes Element
      
      (x.initial <- tapply(data.matrix(abst), h.list, mean)) # siehe help
      dimnames(x.initial) <- list(NULL, dimnames(abst)[[2]])
      
      abst.km <- kmeans(abst, x.initial)
      # Zugehörigkeit
      abst.km$cluster
      # ist das selbe wie
      cutree(hclust(abst.dist, method="average"),3)
    }
  }
  
  # Klassen Analyse (Nächster Nachbar und Naiver Bayes Methode)
  function(){
    library(class)
    library(mda)
    library(ROCR)
    library(klaR)
    ### Daten für Bsp.
    # Ban$Aktiv = Binärer Vektor
    DPfad <- "Studium/5.Semester/STDM - server rkst/Daten/"
    Ban <- read.table(paste(DPfad,"Bankrott.dat",sep=""), header=T)
    
    ### Grafische Darstellung
    pairs(Ban[,1:4], pch=(Ban$Aktiv==0)+1, col=(Ban$Aktiv==0)+3)
    h.cl <- as.integer(oil$zone) # NUR FALLS BUSCHSTABEN UNTERSCHEIDUNG!!!
    pairs(oil[,-1], pch=h.cl+1, col=h.cl+3)
    
    
    ### knn - Verfahren
    ### Nächste Nachbar Methode
    # komplizierte Grenzen, komplexere berechungen als mit Logistischer Analyse
    # k = anzahl Nachbarn (k=1 ist error immer 0, erst gültig ab k=3; keine gerade Zahl, weil nicht Eindeutig)
    Ban.knn1 <- knn(train=Ban[,-5], test=Ban[,-5], cl=as.factor(Ban$Aktiv), k=3)
    
    ### Konfusions-Matrix
    # Spalten = Vorhergesagte Zugehörigkeit : Zeilen = Tatsächliche Zugehörigkeit
    confusion(object=Ban.knn1, true=as.factor(Ban$Aktiv))
    
    ### Plot Nächste Nachbar Methode Confusion
    ## Welches ist die optimale Anzahl weggelassener Nachbarn? knn(k=?)
    function(){
      # wählen der ersten unter den tiefsten Zahlen (ausser 0)
      # y-Achse Fehlerrate mit leave-one-out-Methode
      # x-Achse Anzahl gewählte Nachbaren
      
      x.nfr <- rep(0, 15)
      x.k <- seq(1,2*length(x.nfr), by=2)
      for(i in 1:length(x.k)){
        Ban.knn <- knn(train=Ban[,-5], test=Ban[,-5], cl=as.factor(Ban$Aktiv),k=x.k[i])
        Ban.cm <- confusion(object=Ban.knn, true=as.factor(Ban$Aktiv))
        x.nfr[i] <-  attr(Ban.cm, "error")
      }
      plot(x.k, x.nfr, type="b")
    }
    
    
    ### leave-one-out / Kreuzvalidierung - Method
    ## Direkt error berechnen
    Ban.knn1.cv <- knn.cv(train=Ban[,-5], cl=as.factor(Ban$Aktiv), k=5)
    confusion(object=Ban.knn1.cv, true=as.factor(Ban$Aktiv))
    
    ## einfache Kreuz Validierung
    # Plot leave-one-out Method Confusion
    function(){
      # Vorsicht Resultate schwanken
      knn.loo <- function(train, class, k=1){
        class <- as.factor(class) ## class sollte immer eine Faktovariable sein
        x.cnt <- 0
        for(i in 1:nrow(train)){
          t.knn <- knn(train=train[-i,], cl=class[-i], test=train[i,,drop=F], k=k)
          x.cnt <- x.cnt + (class[i] != t.knn)
        }
        x.cnt/nrow(train)
      }
      x.loo <- rep(0, 15)
      x.k <- seq(1,2*length(x.loo), by=2)
      for(i in 1:length(x.k)){
        x.loo[i] <- knn.loo(train=Ban[,-5], cl=as.factor(Ban$Aktiv), k=x.k[i])
      }
      print("dort wo es am tiefsten ist oder abfalcht wird gewählt")
      print("dieser wert wird erneut in knn() eingesetzt und die confusion berechnet")
      plot(x.k, x.loo, type="b")
    }
    
    ## Kreuz Validierung mit auslassung von 3 oder mehr Beobachtungen
    function(){
      # um einen Bias zu vermeiden. Gleich wie die leave-one-out methode werden Objekte weggelassen.
      # Jedoch wird nicht nur eines sondern geleich mehrere Weggelassen, gemessen wird das in %. x%-iger Weglassung
      # cv = anzahl weggelassene
      
      knn.CV <- function(train, cl, k=1, cv=3){## Version rkst
        ag <- nrow(train) %/% cv
        rest <- nrow(train) %% cv
        if(rest == 0){
          x.rand <- rep(1:ag, rep(cv,ag))[sample(nrow(train))]
        } else{
          x.rand <- c(rep(1:ag, rep(cv,ag)),1:rest)[sample(nrow(train))]
        }
        x.cnt <- 0
        for(i in sort(unique(x.rand))){
          ok <- x.rand != i
          t.knn <-  knn(train=train[ok,], cl=cl[ok],
                        test=train[!ok,, drop=F], k=k)
          x.cnt <- x.cnt + sum(as.vector(cl[!ok] != t.knn))
        }
        x.cnt/nrow(train)
      }
      knn.CV(train=Ban[,-5], cl=as.factor(Ban$Aktiv), cv=3, k=3)
      
      ### Plotten Kreuz Validierung
      # Welches ist die optimale Anzahl wegelassener Beobachtungen
      x.cv3 <- rep(0, 15)
      x.k <- seq(1,2*length(x.cv3), by=2)
      for(i in 1:length(x.k)){
        x.cv3[i] <- knn.CV(train=Ban[,-5], cl=as.factor(Ban$Aktiv), k=x.k[i],cv = 3)
      }
      plot(x.k, x.cv3, type="b")
    }
    
    
    ### Naiver bayes
    function(){
      library(klaR)
      # Grunddaten nach Transfomation falls nötig sind Ban
      Ban.nb <- NaiveBayes(as.factor(Aktiv) ~ ., data = Ban)
      Ban.nb
      par(mfrow=c(2,2)) # Anzahl Attribute die geplotted werden müssen
      # zeichnet für jede Variable die bedingte Naive-Dichte, gegeben die Klassenzugehörigkeit.
      plot(Ban.nb)  ## zeichnet für jede Variable die bedingte NV-Dichte,
      ##               gegeben die Klassenzugehörigkeit.
      
      ## apparent error rate:
      Ban.nbP <- predict(Ban.nb)
      require(mda)
      confusion(object=Ban.nbP$class, true=as.factor(Ban$Aktiv))
      
      ## Naive Bayes (leave-one-out):
      NB.loo <- function(formula, data){
        y.p <- rep(NA, nrow(data))
        for(i in 1:nrow(data)){
          t.nb <- NaiveBayes(formula, data=data[-i,])
          y.p[i] <- predict(t.nb, newdata=data[i,, drop=F])$class
        }
        y.p
      }
      Ban.nbLOO <- NB.loo(as.factor(Aktiv) ~ ., data = Ban) # Aktiv = Name der Klassen Zeile
      confusion(object=Ban.nbLOO, true=Ban$Aktiv+1)
      
    }
    
    ### Spezivität/Sensitivität (Performance-Mass)
    function(){
      ### leave-one-out method
      Ban.knn1.cv <- knn.cv(train=Ban[,-5], cl=as.factor(Ban$Aktiv), k=1)
      conf <- confusion(Ban.knn1.cv, as.factor(Ban$Aktiv))
      
      Spez.Sens <- function(conf){
        Spezivitat <- conf[1]/(conf[1]+conf[2])
        Sensitivitat <- conf[4]/(conf[4]+conf[3])
        S <- list()
        S[[1]] <- Spezivitat; names(S[[1]]) <- "Spezivität"
        S[[2]] <- Sensitivitat; names(S[[2]]) <- "Sensitivität"
        return(S)
      }
      Spez.Sens(conf)
      
      # Beschreibung wie es funktioniert
      function(){
        Spezivität => (oben-links) / (oben-links + oben-rechts)
        # Spezivität = TN / N
        17/(17+4)  ## = 0.8095 d.h. ca 81%
        Sensitivität => (unten-rechts) / (unten-rechts + unten-links)
        # Sensitivität (=Trefferquote) = TP/P
        21/(4+21)  ## = 0.84, d.h. ca 84%
      }
    }
    ## Berechne Matthew's correlation coeffizient aus der Konfusionsmatrix:
    f.mcc <- function(pred, true){
      require(mda)
      cm <- confusion(object=pred, true=true)
      h1 <- sum(cm[,2])*sum(cm[2,])
      h2 <- sum(cm[1,])*sum(cm[,1])
      if(h1==0 | h2==0) return(0)
      else return(((cm[2,2]*cm[1,1] - cm[2,1]*cm[1,2]))/ sqrt(h1)/sqrt(h2))
    }
    f.mcc(pred=WWVersScore$Cl1class, true=WWVersScore$Purchase) ##  = 0.191353
    ## oder über Chiquadrat-Test:
    cv.test <- function(pred,true) {
      CV <- sqrt(chisq.test(pred, true, correct=FALSE)$statistic /
                   (length(pred) * (min(length(unique(pred)),length(unique(true))) - 1)))
      print.noquote("Cram?r V / Phi:")
      return(as.numeric(CV))
    }
    cv.test(pred=WWVersScore$Cl1class, true=WWVersScore$Purchase) ## = 0.191353
    
    ### Performance-Vergleich mit Liftchart:
    # Vergelich: der Höhere ist besser
    plot.liftchart <- function(true, predScore, quant=(1:10)/10, add=FALSE, col="blue", cum=FALSE, ylim=NULL){
      ## true = 0 oder 1 (hohe Score --> 1)
      
      ## true sortieren nach den Scorewerten:
      trueSort <- true[order(predScore, decreasing=TRUE)]
      
      ## Indices der entsprechenden Quantile suchen
      QIndex <- round(quantile(1:length(trueSort), quant))
      
      ## Fehlerraten an den Quantilpunkten bestimmen
      rate <- cRate <- rep(NA, length(QIndex))
      for (i in 1:length(QIndex)){
        rate[i] <- mean(trueSort[1:QIndex[i]])
        cRate[i] <- sum(trueSort[1:QIndex[i]])
      }
      ##
      rRate <- rev(rate)[1]
      if(add){
        if(cum) lines(quant, cRate/rev(cRate)[1], col=col, type="b")
        else lines(quant, rate/rRate, col=col, type="b")
      }
      else{
        if(cum){
          plot(quant, cRate/rev(cRate)[1], lwd=2, col=col, type="b",
               xlim=c(0,1), ylim=c(0,1), las=1,
               xlab="Quantile", ylab="Kumulativer Lift")
          abline(a=0,b=1, lty=2)
        }
        else{
          plot(quant, rate/rRate, lwd=2, col=col, type="b", las=1,
               xlab="Quantile", ylab="Lift", ylim=ylim)
          abline(h=1, lty=2)
        }
      }
    }
    # Bsp.
    h.true <- as.integer(WWVersScore$Purchase)-1 # true muss numerisch im bereich von 0 bis 1 sein!!!!
    plot.liftchart(true=h.true, predScore=WWVersScore$Cl1, quant=c(0.02, 0.05, (1:10)/10), ylim=c(0,8))
    ## Bis zu einer Selektion der besten 20% ist der Klassifizierer 1 klar besser
    ## als der Klassifizierer 2. Danach spielt es keine Rolle mehr.
  }
  
  # Klassische Klassifikationsmethoden (logistischer Klassifizierer , Klassifikations Bäume)
  function(){
    library(MASS); library(nnet); library(rpart)
    
    # Beispiel Daten
    DPfad <- "Studium/5.Semester/STDM - server rkst/Daten/"
    Ban <- read.table(paste(DPfad,"Bankrott.dat",sep=""), header=T)
    
    
    
    ### logistischer Klassifizierer
    # Grenzen werden einfach gezogen die Methode ist relativ simpel
    # simpler als knn Verfahren!
    # Es werden kategorisierte Variablen bevorzugt
    
    Ban.mn <- multinom(Aktiv ~ . , data=Ban) # Aktiv = klassifizierung
    ## final  value 13.721620     => wird ausgegeben nicht in Variable gespeichert
    confusion(true=Ban$Aktiv, object=predict(Ban.mn))
    
    ## Leave-one-out (Out of Sample)
    t.pr <- rep(NA, nrow(Ban))
    for(i in 1:nrow(Ban)){  ##
      t.mn <- multinom(Aktiv ~ . , data=Ban[-i,], trace=FALSE)
      t.pr[i]  <- as.character(predict(t.mn, newdata=Ban[i,,drop=FALSE]))
      # print(paste(i/nrow(Ban)*100,"%",sep = ""))
    }
    confusion(true=Ban$Aktiv, object=t.pr)
    
    ## AIC Variablen selektion (Modellvereinfachung)
    # NUR Coeffizienten am Schluss sollen benutzt werden (ohne intercept)
    step(Ban.mn)
    Ban.mn2 <- multinom(Aktiv ~ CFpVs + WspVk, data=Ban)
    # Resultat normal in confusion und leave-one-out einsetzten
    
    
    
    ### Klassifikations Bäume
    # Es werden stetige Rohdaten bevorzugt, (kategoriesieren vernichtet Informationen)
    
    # default einstellungen mit nur einer frage
    # basierend auf der Entropie (standart: split="information")
    Ban.rpG <- rpart(Aktiv ~ ., data=Ban, method="class", parms=list(split='information'))
    # Gini-Index
    Ban.rpG <- rpart(Aktiv ~ ., data=Ban, method="class", parms=list(split='gini'))
    
    ## Ausgewachsener Klassifikationsbaum:
    # Gini-Index                  # parms=list(split='gini'),
    # basierend auf der Entropie  # parms=list(split='information'),
    Ban.rp <- rpart(Aktiv ~ ., data=Ban, method="class",
                    parms=list(split='gini'),
                    control=rpart.control(minsplit=1, minbucket=1, cp=0.00))
    
    ## Plotten
    plot(Ban.rp, branch=0.6, uniform=T)
    text(Ban.rp, use.n=TRUE)
    
    # Conusion
    Ban.rpGp <- predict(Ban.rp, newdata=Ban, type="class")
    confusion(true=Ban$Aktiv, object=Ban.rpGp)
    
    ## leave-one-out
    t.prV <- rep(NA, nrow(Ban)) # anzahl Zeilen des Datensatzes
    i <- 1
    for(i in 1:nrow(Ban)){
      t.rpV <- rpart(Aktiv ~ ., data=Ban[-i,], method="class", parms=list(split='gini'),
                     control=rpart.control(minsplit=1,minbucket=1,cp=0.00)) # Baum erstellen aus Modell ohne Zeile i
      t.rpP <- predict(t.rpV, newdata=Ban[i,,drop=F], type="class") # predicten der Zeilte i
      t.prV[i] <- t.rpP
      # print(paste(i/nrow(Ban)*100,"%",sep = ""))
    }
    confusion(true=as.numeric(as.factor(Ban$Aktiv)), object=t.prV)
    
    ## Pruning zurückstutzen des Baumes
    # Linie wird bei der kleinsten Messung am oberen ende des Intervalls gezogen
    # die erste Messung, welche unterhalb dieser Linie liegt wird gewählt
    # x-Achse = Anzahl Abstufungen des Ausgangsmodells
    plotcp(Ban.rp) ## cp = Delta I im Skript
    
    Ban.rpIvp <- prune(Ban.rp, cp=0.27) # Ausgelesene x-Achse aus Plot eingeben
    par(mfrow=c(1,1), mar=c(1,1,1,1), xpd=T) ## xpd=T: otherwise on some devices the text is clipped
    plot(Ban.rpIvp, branch=0.6, uniform=T); text(Ban.rpIvp, use.n=TRUE)
    
  }
  
  ### randomForest
  function(){
    library(randomForest)
    
    <### randomForest
    # Stetige Variablen (Messungen) werden bevorzugt, weil kategorisierte Variablen Infomrationen vernichten
    
    (Ban.rf <- randomForest(as.factor(Aktiv) ~ ., data=Ban))
    # as.factor() => muss mit der Kategoriellen Variable gemacht werden (Einteilung in Gruppen)
    # Wenn nicht auf alle randomForest gemacht wird, PUNKT weglassen und
    #  randomForest(as.factor(Aktiv) ~ X1 + X2 + X3,data=Datensatz) schreiben
    
    Ban.rf
    # Gibt u.a. die Out-of-Bag-Fehlerrate (OOB) an ???
    
    ## die Wichtigkeit der Varibale im random Forest
    importance(Ban.rf)
    # grafische Darstellung (nur bei vielen Variablen nützlich)
    varImpPlot(Ban.rf)
    # Die Variablen mit hohen Werten sind wichtig
  }
  
  
  ### Boosting (gebooster Klassifikationsbaum)
  function(){
    library(gbm)
    # interaction.depth = ???
    # distribution = "adaboost"  => dann wird AdaBoost-Verfahren verwendet.
    Ban.gbm <- gbm(as.integer(Aktiv) ~ ., data=Ban, distribution="bernoulli",
                   interaction.depth = 3,  train.fraction = 0.8, cv.folds=5,
                   n.minobsinnode=5, n.trees = 5000, verbose=FALSE)
    
    
    summary(Ban.gbm)
    # Gibt Importance-Mass an und erzeugt entsprechende Grafik
    # Aus diesem Output muss man schliessen, welche Variablen wichtig sind, um die Klassenlabels vorherzusagen.
    # Was heiss die Klassenlabes vorher zu sagen????
    
    ## Optimale Anzahl Bäume (jeweils inklusive plot)
    par(mfrow=c(1,3))
    gbm.perf(Ban.gbm, method="OOB")   # unterschätzt gemäss Warnung die optimale Anzahl tendenziell (Warnung)
    gbm.perf(Ban.gbm, method="test")
    gbm.perf(Ban.gbm, method="cv")
    # --> verwenden gerundete Anzahl Bäume
    
    ## Anwenden der anzahl Bäume
    Ban.gbm.P <- predict(Ban.gbm, newdata=Ban, n.trees=2000, type="response")
    
    require(mda)
    confusion(true=Ban$Aktiv, object=round(Ban.gbm.P,0))
    
    ### Confusion-Matrix Vergleich mit der Leave on out Methode
    t.pr <- rep(NA, nrow(Ban))
    x.cnt <- 0
    for(i in 1:nrow(Ban)){  ##
      t.gbm <- gbm(as.integer(Aktiv) ~ ., data=Ban[-i,], distribution="bernoulli",
                   interaction.depth = 3,  train.fraction = 0.8,
                   n.minobsinnode=5, n.trees = 2000, verbose=FALSE)
      t.pr[i] <- round(predict(t.gbm, newdata=Ban[i,,drop=FALSE], n.trees=2000,
                               type="response"),0)
      x.cnt <- x.cnt + (Ban$Aktiv[i] != t.pr[i])
    }
    confusion(true=Ban$Aktiv, object=t.pr)
    
  }
}

# MF 2
function(){
  
  # Zinsrechnen
  function(){ 
    ###Legende
    # c = Coupon
    # N = Nominal wert (Anfangswert)
    # n = anz. Schritte
    # s = Kassazinsen
    # f_1,2 = Terminzinssatz
    
    r <- c(5.0, 5.3, 5.6, 5.8, 6.0, 6.1) #Bsp
    s <- c(0.05, 0.053, 0.056, 0.058, 0.06, 0.061) # in 100tel schreiben
    
    diskontfaktor <- function(prozent,dauer){
      return((1/(1+prozent))^dauer)
    }
    diskont.s <- diskontfaktor(s,1:length(s)) #normaler diskontfaktor
    diskont.r <- diskontfaktor(r,1) #normaler diskontfaktor 
    diskont.s <- cumprod(diskont.r) # Umwandeln
    
    # Kassazins
    sfc <- function(rates){
      sapply(1:length(rates), function(i) (prod(rates[1:i]+1)^(1/i))-1)
    }
    (s <- sfc(r))
    
    fij <- function(si,i,j){
      (((1+si[j])^j / (1+si[i])^i)^(1/(j-i)))-1
    }
    (f <- fij(s, von, bis))
    
    # Kassazinskurve
    (sstrich <- fij(s,1,1:length(s)))
    
    # Short-Rates (von Jahr(i) bis Jahr(i+1))
    rfc <- function(s){
      fij <- function(si,i,j){
        (((1+si[j])^j / (1+si[i])^i)^(1/(j-i)))-1
      }
      x <- sapply(1:length(s), function(i) (fij(s,i,i+1)))
      return(c(s[1],x[1:length(x)-1]))
    }
    rfc(s)

    # lambda Rückverzinsung
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
    
    # Zahlungsströme eingeben
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
  }
  
  # F = Terminpreis
  # S = Kassapreis
  
  # M = Endzeitpunkt
  # r = Zins (0.05) pro Jahr
  # k = Startzeitpunkt
  # m = unterjährig verzinst
  
  S <- 412; M <- 3; c <- 0.5; r <- 0.0225; m <- 1 # Bsp Werte
  
  # Abdiskontieren
  d <- function(M,r,k,m) return( 1/(1+(r/m))^(M-k) )
  d(M,r,5,m)
  
  # Terminpreis
  F.from.S <- function(S,M,c,r,m){
    if(length(c)==1){c <- rep(c,M)}
    d <- function(M,r,k,m) return( 1/(1+(r/m))^(M-k) )
    return( S/d(M,r,0,m) + sum(c/d(M,r,0:(M-1),m)) )
  }
  F <- F.from.S(S,M,c,r,m)
  
  # Kassapreis
  S.from.F <- function(F,M,c,r,m) {
    if(length(c)==1){c <- rep(c,M)}
    d <- function(M,r,k,m) return( 1/(1+(r/m))^(M-k) )
    return(d(M,r,0,m)*F - sum(d(M,r,0:(M-1),m)*c))
  }
  S.from.F(F,M,c,r,m)
  
  ### Optimieren nach Zins
  function(){
    ## Formel für Kassapreis aus Terminpreis
    S.from.F <- function(F,M,c,r) {
      d = 1/(1+r/12)^(0:M)
      d[M+1]*F - sum(d[1:M])*c/12
    }
    ## Preis-Varianz bei unterschiedlichen Maturitäten
    varS <- function(r) {
      F = c(406.5, 416.64, 423.48, 433.84)
      M = c(1,4,6,9)
      c = 20
      S = numeric(4)
      for (i in 1:4) {
        S[i] = S.from.F(F[i],M[i],c,r)
      }
      var(S)
    }
    ## Minimierung der Varianz mit Optimierungsalgorithmus
    res = optimize(f = varS, lower = 0, upper = 1)
    res$minimum # r = 0.05001842
    res$objective # varS(r) = 8.127062e-06
  }
  
  multi.model <- function(dT,title) {
    S0 = 10
    nu = 0.15
    sigma = 0.4
    mu = nu + sigma^2/2
    # dT = 1/52
    nStep = 10/dT
    nStep
    dZ = sqrt(dT) * rnorm(nStep)
    time = seq(0,nStep)*dT
    # Methode 1:
    # Mit It^o-Standardform
    # Drift: mu
    Klammer.eckig = 1 + mu*dT + sigma*dZ
    S1 = cumprod( c(S0,Klammer.eckig) )
    # Methode 2:
    # Mit Gl. fur log. Preise
    # Drift: nu
    Klammer.geschweift = nu*dT + sigma*dZ
    S2 = cumprod(c(S0,exp(Klammer.geschweift)))
    plot.ts(ts(cbind(S1,S2)),main=title,plot.type="single",col=c(1,2),ylab="S")
    legend("topleft", c("normal","log-normal"), col=c(1,2),lty=c(1,1))
    mean(S1)-mean(S2)
  }
  par(mfrow=c(1,2))
  multi.model(dT = 1/52,title="multiplikatives Modell, dt=1/52")
  multi.model(dT = 1/250,title="multiplikatives Modell, dt=1/250")
  
  ### Binomainal Tree
  function(){
    # binomial trees
    
    # u = up
    # d = down
    # per = anzahl Schritte
    # S = S(0) = Anfangsbetrag
    # ir = r (zins in 0.05 form)
    # K = Normalfall = 0 (pmax(tree[[as.character(per)]]-K,0))
    # sig = Sigma
    # dT = delta t = ∆t
    # payout = pmax(Tree$'per'-K,0)
    
    createBTree = function(u,d=1/u,per=1,S=1){
      tree = list()
      tree[["0"]] = S
      for (i in 1:per) {
        x = tree[[as.character(i-1)]]
        x = c(u*x,d*x[i])
        tree[[as.character(i)]] = x
      }
      tree
    }
    createBTreeSym = function(sig,dT,per=1,S=1){
      u = exp(sig*sqrt(dT))
      d = 1/u
      tree = list()
      tree[["0"]] = S
      for (i in 1:per) {
        x = tree[[as.character(i-1)]]
        x = c(u*x,d*x[i])
        tree[[as.character(i)]] = x
      }
      tree
    }
    
    # European und American, aber komisch
    function(){    
      # european option type valuation for arbitraty payout 
      # given as a payout vector
      option.european = function(tree,ir,payout,per){
        if (missing(per)) { per = length(tree)-1 }
        if (length(payout)!=length(tree))
          stop("Length of payout vector don't fit size of tree.")
        # compute risk neutral probabilities
        u = tree[[2]][1]/tree[[1]]
        d = tree[[2]][2]/tree[[1]]
        q = (1+ir-d)/(u-d)
        omq = 1-q
        df = 1/(1+ir)
        # print(c(q,df))
        otree = list()
        y = payout
        otree[[as.character(per)]] = y
        for (i in per:1) {
          x = df*(q * y[-(i+1)] + omq * y[-1])
          otree[[as.character(i-1)]] = x
          y = x
        }
        otree
      }
      
      # american option type valuation for arbitraty payout
      # given as payout lattice
      # tree=Tree; ir = 0.01;tree.payout=pmax(Tree$'5'-50,0);per=5;K=50
      option.call.american = function(tree,ir,tree.payout,per,K=0){
        if (missing(per)) { per = length(tree)-1 }
        
        # compute risk neutral probabilities
        u = tree[[2]][1]/tree[[1]]
        d = tree[[2]][2]/tree[[1]]
        q = (1+ir-d)/(u-d)
        omq = 1-q
        df = 1/(1+ir)
        # print(c(q,df))
        otree <- tree.payout
        y <- otree
        #         list()
        #         y = pmax(tree[[as.character(per)]]-K,0)
        #         otree[[as.character(per)]] = y
        for (i in per:1) {
          x = df*(q * y[-(i+1)] + omq * y[-1])
          x = pmax(x,tree[[ as.character(i-1)]] -K)
          otree[[as.character(i-1)]] = x
          y = x
        }
        otree
      }
      option.call.american = function(tree.price,ir,tree.payout,per){
        if (missing(per)) { per = length(tree)-1 }
        
        # compute risk neutral probabilities
        u = tree[[2]][1]/tree[[1]]
        d = tree[[2]][2]/tree[[1]]
        q = (1+ir-d)/(u-d)
        omq = 1-q
        df = 1/(1+ir)
        # print(c(q,df))
        otree = tree.payout
        #   list()
        #   y = pmax(tree[[as.character(per)]]-K,0)
        #   otree[[as.character(per)]] = y
        for (i in per:1) {
          x = df*(q * y[-(i+1)] + omq * y[-1])
          x = pmax(x,tree[[ as.character(i-1)]] -K)
          otree[[as.character(i-1)]] = x
          y = x
        }
        otree
      }
      
      ### Beispiel
      (Tree <- createBTree(u=1.07484,S=53,per=5))
      createBTreeSym(sig = 0.25,dT = 1/12,per = 5,S = 53)
      option.european(tree = Tree,ir = 0.01,payout = pmax(Tree$'5'-50,0),per = 5)
      option.call.american(tree= Tree,ir = 0.01,tree.payout = pmax(Tree$'5'-50,0),per = 5)
      
    }
    
    # Wann lohnt es sich zu verkaufen?
    function(){
      # ir = r / dt
      baum <- createBTree(u = 1.07484,d = 0.93037,per = 5,S = 53)
      baum.mat <- tree2matrix(baum)
      
      eu <- option.call(tree = baum,ir = 0.01,K = 50,per = 5)
      am <- option.call.american(tree = baum,ir = 0.01,K = 50,per = 5)
      di <- option.call.digital(tree = baum,ir = 0.01,K = 50,per = 5)
      
      (eu.mat <- tree2matrix(eu,type="valuation"))
      (am.mat <- tree2matrix(am,type="valuation"))
      (di.mat <- tree2matrix(di,type="valuation"))
      
      pay <- K-baum.mat-eu.mat
      pay[which(pay<0)] <- NA
      pay
      
      pay <- K-baum.mat-am.mat
      pay[which(pay<0)] <- NA
      pay
    }
    
    # ir: Zinssatz PRO Periode r/dt
    option.call = function(tree,ir,K,per){
      if (missing(per)) { per = length(tree)-1 }
      
      # compute risk neutral probabilities
      u = tree[[2]][1]/tree[[1]]
      d = tree[[2]][2]/tree[[1]]
      q = (1+ir-d)/(u-d)
      omq = 1-q
      df = 1/(1+ir)
      # print(c(q,df))
      otree = list()
      y = pmax(tree[[as.character(per)]]-K,0)
      otree[[as.character(per)]] = y
      for (i in per:1) {
        x = df*(q * y[-(i+1)] + omq * y[-1])
        otree[[as.character(i-1)]] = x
        y = x
      }
      otree
    }
    option.call.digital = function(tree,ir,K,per){
      if (missing(per)) { per = length(tree)-1 }
      
      # compute risk neutral probabilities
      u = tree[[2]][1]/tree[[1]]
      d = tree[[2]][2]/tree[[1]]
      q = (1+ir-d)/(u-d)
      omq = 1-q
      df = 1/(1+ir)
      # print(c(q,df))
      otree = list()
      y = pmax(0.5*(1+sign(tree[[as.character(per)]]-K)),0)
      otree[[as.character(per)]] = y
      for (i in per:1) {
        x = df*(q * y[-(i+1)] + omq * y[-1])
        otree[[as.character(i-1)]] = x
        y = x
      }
      otree
    }
    option.call.american = function(tree,ir,K,per){
      if (missing(per)) { per = length(tree)-1 }
      
      # compute risk neutral probabilities
      u = tree[[2]][1]/tree[[1]]
      d = tree[[2]][2]/tree[[1]]
      q = (1+ir-d)/(u-d)
      omq = 1-q
      df = 1/(1+ir)
      # print(c(q,df))
      otree = list()
      y = pmax(tree[[as.character(per)]]-K,0)
      otree[[as.character(per)]] = y
      for (i in per:1) {
        x = df*(q * y[-(i+1)] + omq * y[-1])
        x = pmax(x,tree[[ as.character(i-1)]] -K)
        otree[[as.character(i-1)]] = x
        y = x
      }
      otree
    }
    option.call.asian <- function(tree, K, ir, dT){ 
      u = tree[[2]][1]/tree[[1]]
      d = tree[[2]][2]/tree[[1]]
      per <- ncol(tree)-1
      up = u
      down = 1/u
      R = 1+ir*dT
      q = (R-down)/(up-down)
      
      nPath <- 2^(per)
      
      payout.expected <- 0
      for (i in 0:(nPath-1)) { #i<-0
        tmp <- as.numeric(digitsBase(i,base=2,per))
        n.down <- sum(tmp)
        n.up <- per - n.down
        idx <- cbind(
          cumsum(c(1,tmp)),
          1:(length(tmp)+1)
        )
        S.mean <- mean(tree[idx[-1,]])
        payout <- max(S.mean-K,0)
        payout.expected <- payout.expected + payout * q^n.up * (1-q)^n.down
      }
      return(payout.expected/R^per)
    }
    
    option.put = function(tree,ir,K,per){
      if (missing(per)) { per = length(tree)-1 }
      
      # compute risk neutral probabilities
      u = tree[[2]][1]/tree[[1]]
      d = tree[[2]][2]/tree[[1]]
      q = (1+ir-d)/(u-d)
      omq = 1-q
      df = 1/(1+ir)
      #print(c(u,d,q,df))
      otree = list()
      y = pmax(K-tree[[as.character(per)]],0)
      otree[[as.character(per)]] = y
      for (i in per:1) {
        x = df*(q * y[-(i+1)] + omq * y[-1])
        otree[[as.character(i-1)]] = x
        y = x
      }
      otree
    }
    option.put.digital = function(tree,ir,K,per){
      if (missing(per)) { per = length(tree)-1 }
      
      # compute risk neutral probabilities
      u = tree[[2]][1]/tree[[1]]
      d = tree[[2]][2]/tree[[1]]
      q = (1+ir-d)/(u-d)
      omq = 1-q
      df = 1/(1+ir)
      #print(c(u,d,q,df))
      otree = list()
      y = pmax(0.5*(1+sign(K-tree[[as.character(per)]])),0)
      otree[[as.character(per)]] = y
      for (i in per:1) {
        x = df*(q * y[-(i+1)] + omq * y[-1])
        otree[[as.character(i-1)]] = x
        y = x
      }
      otree
    }
    option.put.american = function(tree,ir,K,per){
      if (missing(per)) { per = length(tree)-1 }
      
      # compute risk neutral probabilities
      u = tree[[2]][1]/tree[[1]]
      d = tree[[2]][2]/tree[[1]]
      q = (1+ir-d)/(u-d)
      omq = 1-q
      df = 1/(1+ir)
      # print(c(q,df))
      otree = list()
      y = pmax(K-tree[[as.character(per)]],0)
      otree[[as.character(per)]] = y
      for (i in per:1) {
        x = df*(q * y[-(i+1)] + omq * y[-1])
        x = pmax(x,K-tree[[ as.character(i-1) ]])
        otree[[as.character(i-1)]] = x
        y = x
      }
      otree
    }
    
    # turns a tree into a matrix or back again
    tree2matrix <- function(tree, type="base", layout="upper"){
      d = length(tree)
      m = matrix(nrow=d,ncol=d)
      if (layout == "upper") {
        for (i in 1:d)
        {
          if (type=="base") {
            m[1:i,i] = tree[[i]]
          } else if (type=="valuation") {
            m[1:i,i] = tree[[d-i+1]]
          }
        }
        dimnames(m)[[1]] <- rep(" ",d)
        
      } else if (layout == "lower" ) {
        for (i in 0:(d-1))
        {
          if (type=="base") {
            m[(d-i):d,i+1] = tree[[i+1]]
          } else if (type=="valuation") {
            m[(d-i):d,i+1] = tree[[d-i]]
          }
        }
        dimnames(m)[[1]] <- (d-1):0    
      } else
        stop("Unkown layout")
      dimnames(m)[[2]] <- as.character(0:(d-1))
      m
    }
    tree2matrix(option.call,type="valuation")
    matrix2tree <- function(mat, type="base", layout="upper"){
      dimnames(mat) = list(NULL,NULL)
      tree <- list()
      d <- dim(mat)[1]
      if (layout=="upper") {
        for (i in 1:d) {
          tree[[as.character(i-1)]] <- mat[1:i,i] 
        }
      } else
        stop ("Not yet implemented.")
      tree
    }
    
    bewertung.mit.pacht <- function(payout,q,ir,per,pacht=0,dig=3,emph=c("\fbox{","}")){
      # compute risk neutral probabilities
      omq <- 1-q
      df <- 1/(1+ir)
      
      d <- dim(payout)[1]
      wert.mat <- matrix(nrow=d, ncol=d)
      exercise.mat <- matrix(nrow=d,ncol=d)
      exercise.mat[,] <- FALSE
      
      exercise.mat[,d] <- TRUE
      wert.mat[,d] <- payout[,d]
      for (i in ((d-1):1)) {
        w.vec <-df*(q * wert.mat[1:i,i+1] + omq * wert.mat[2:(i+1),i+1]) - pacht
        wert.mat[1:i,i] <- pmax(w.vec, payout[1:i,i])
        exercise.mat[1:i,i] <- w.vec < payout[1:i,i]
      }
      out = as.character(round(wert.mat,dig))
      dim(exercise.mat) = NULL
      for (i in 1:length(out)) {
        if (exercise.mat[i]) {
          out[i] = paste(emph[1],out[i],emph[2],sep="")
          out[is.na(wert.mat)] == ""
        }
      }
      dim(out) = dim(wert.mat)
      #    list(value=wert.mat, exercise=exercise.mat)
      out
    }
  }
  # Options Vinc
  function(){
    ## Binom Funktion
    # R = 1+ir
    
    future.prices <- function(tree, q){
      mat <- matrix(data = NA,nrow = nrow(tree),ncol = ncol(tree))
      mat[,ncol(tree)] <- tree[,ncol(tree)]
      for(i in ncol(tree) : 2){
        for( j in nrow(tree):2 ){
          mat[j,i-1] <- q*mat[j-1,i]+(1-q)*mat[j,i]
        }
      }
      return(mat)  
    }
    future.prices(prices$baum.mat,2/3)
    
    options.prices <- function(S,K,u,d=1/u,ir,per,div=0){
      # Functions
      createBTree = function(u,d=1/u,per=1,S=1){
        tree = list()
        tree[["0"]] = S
        for (i in 1:per) {
          x = tree[[as.character(i-1)]]
          x = c(u*x,d*x[i])
          tree[[as.character(i)]] = x
        }
        tree
      }
      createBTreeSym = function(sig,dT,per,S){
        u = exp(sig*sqrt(dT))
        d = 1/u
        tree = list()
        tree[["0"]] = S
        for (i in 1:per) {
          x = tree[[as.character(i-1)]]
          x = c(u*x,d*x[i])
          tree[[as.character(i)]] = x
        }
        tree
      }
      tree2matrix <- function(tree, type="base", layout="upper"){
        d = length(tree)
        m = matrix(nrow=d,ncol=d)
        if (layout == "upper") {
          for (i in 1:d)
          {
            if (type=="base") {
              m[1:i,i] = tree[[i]]
            } else if (type=="valuation") {
              m[1:i,i] = tree[[d-i+1]]
            }
          }
          dimnames(m)[[1]] <- rep(" ",d)
          
        } else if (layout == "lower" ) {
          for (i in 0:(d-1))
          {
            if (type=="base") {
              m[(d-i):d,i+1] = tree[[i+1]]
            } else if (type=="valuation") {
              m[(d-i):d,i+1] = tree[[d-i]]
            }
          }
          dimnames(m)[[1]] <- (d-1):0    
        } else
          stop("Unkown layout")
        dimnames(m)[[2]] <- as.character(0:(d-1))
        m
      }
      option.call = function(tree,ir,K,per){
        if (missing(per)) { per = length(tree)-1 }
        
        # compute risk neutral probabilities
        u = tree[[2]][1]/tree[[1]]
        d = tree[[2]][2]/tree[[1]]
        q = (1+ir-d)/(u-d)
        omq = 1-q
        df = 1/(1+ir)
        # print(c(q,df))
        otree = list()
        y = pmax(tree[[as.character(per)]]-K,0)
        otree[[as.character(per)]] = y
        for (i in per:1) {
          x = df*(q * y[-(i+1)] + omq * y[-1])
          otree[[as.character(i-1)]] = x
          y = x
        }
        otree
      }
      option.call.american = function(tree,ir,K,per){
        if (missing(per)) { per = length(tree)-1 }
        
        # compute risk neutral probabilities
        u = tree[[2]][1]/tree[[1]]
        d = tree[[2]][2]/tree[[1]]
        q = (1+ir-d)/(u-d)
        omq = 1-q
        df = 1/(1+ir)
        # print(c(q,df))
        otree = list()
        y = pmax(tree[[as.character(per)]]-K,0)
        otree[[as.character(per)]] = y
        for (i in per:1) {
          x = df*(q * y[-(i+1)] + omq * y[-1])
          x = pmax(x,tree[[ as.character(i-1)]] -K)
          otree[[as.character(i-1)]] = x
          y = x
        }
        otree
      }
      option.put = function(tree,ir,K,per){
        if (missing(per)) { per = length(tree)-1 }
        
        # compute risk neutral probabilities
        u = tree[[2]][1]/tree[[1]]
        d = tree[[2]][2]/tree[[1]]
        q = (1+ir-d)/(u-d)
        omq = 1-q
        df = 1/(1+ir)
        #print(c(u,d,q,df))
        otree = list()
        y = pmax(K-tree[[as.character(per)]],0)
        otree[[as.character(per)]] = y
        for (i in per:1) {
          x = df*(q * y[-(i+1)] + omq * y[-1])
          otree[[as.character(i-1)]] = x
          y = x
        }
        otree
      }
      option.put.american = function(tree,ir,K,per,div=div){
        if (missing(per)) { per = length(tree)-1 }
        
        # compute risk neutral probabilities
        u = tree[[2]][1]/tree[[1]]
        d = tree[[2]][2]/tree[[1]]
        q = (1+ir-div-d)/(u-d)
        omq = 1-q
        df = 1/(1+ir)
        # print(c(q,df))
        otree = list()
        y = pmax(K-tree[[as.character(per)]],0)
        otree[[as.character(per)]] = y
        for (i in per:1) {
          x = df*(q * y[-(i+1)] + omq * y[-1])
          x = pmax(x,K-tree[[ as.character(i-1) ]])
          otree[[as.character(i-1)]] = x
          y = x
        }
        otree
      }
      
      
      # Prices
      baum <- createBTree(u = u,d = d,per = per,S = S)
      baum.mat <- tree2matrix(baum)
      
      # Calloption
      o.c.e <- option.call(baum,ir=ir,K=K,per=per)
      o.c.e.mat <- tree2matrix(o.c.e,type="valuation") # Als Matrix
      bin.c.e <- baum.mat-K-o.c.e.mat
      bin.c.e[which(bin.c.e<0)] <- NA
      o.c.a <- option.call.american(baum,ir=ir,K=K,per=per)
      o.c.a.mat <- tree2matrix(o.c.a,type="valuation") #  Als Matrix
      bin.c.a <- baum.mat-K-o.c.a.mat
      bin.c.a[which(bin.c.a<0)] <- NA
      
      # Putoption
      o.p.e <- option.put(baum,ir=ir,K=K,per=per)
      o.p.e.mat <- tree2matrix(o.p.e,type="valuation")
      bin.p.e <- K-baum.mat-o.p.e.mat
      bin.p.e[which(bin.p.e<0)] <- NA
      o.p.a <- option.put.american(baum,ir=ir,K=K,per=per,div=0)
      o.p.a.mat <- tree2matrix(o.p.a,type="valuation")
      bin.p.a <- K-baum.mat-o.p.a.mat
      bin.p.a[which(bin.p.a<0)] <- NA
      
      #Output
      return(list(baum.mat=baum.mat,o.p.e.mat=o.p.e.mat,bin.p.e=bin.p.e,o.p.a.mat=o.p.a.mat,bin.p.a=bin.p.a
                  ,o.c.e.mat=o.c.e.mat,bin.c.e=bin.c.e,o.c.a.mat=o.c.a.mat,bin.c.a=bin.c.a,K=K))
      
    }
    options.prices <- function(S,K,sigma,ir,per,dT,div=0){
      # Functions
      createBTree = function(u,per,S){
        d <- 1/u
        tree = list()
        tree[["0"]] = S
        for (i in 1:per) {
          x = tree[[as.character(i-1)]]
          x = c(u*x,d*x[i])
          tree[[as.character(i)]] = x
        }
        tree
      }
      createBTreeSym = function(sig,dT,per,S){
        u = exp(sig*sqrt(dT))
        d = 1/u
        tree = list()
        tree[["0"]] = S
        for (i in 1:per) {
          x = tree[[as.character(i-1)]]
          x = c(u*x,d*x[i])
          tree[[as.character(i)]] = x
        }
        tree
      }
      tree2matrix <- function(tree, type="base", layout="upper"){
        d = length(tree)
        m = matrix(nrow=d,ncol=d)
        if (layout == "upper") {
          for (i in 1:d)
          {
            if (type=="base") {
              m[1:i,i] = tree[[i]]
            } else if (type=="valuation") {
              m[1:i,i] = tree[[d-i+1]]
            }
          }
          dimnames(m)[[1]] <- rep(" ",d)
          
        } else if (layout == "lower" ) {
          for (i in 0:(d-1))
          {
            if (type=="base") {
              m[(d-i):d,i+1] = tree[[i+1]]
            } else if (type=="valuation") {
              m[(d-i):d,i+1] = tree[[d-i]]
            }
          }
          dimnames(m)[[1]] <- (d-1):0    
        } else
          stop("Unkown layout")
        dimnames(m)[[2]] <- as.character(0:(d-1))
        m
      }
      option.call = function(tree,ir,K,per){
        if (missing(per)) { per = length(tree)-1 }
        
        # compute risk neutral probabilities
        u = tree[[2]][1]/tree[[1]]
        d = tree[[2]][2]/tree[[1]]
        q = (1+ir-d)/(u-d)
        omq = 1-q
        df = 1/(1+ir)
        # print(c(q,df))
        otree = list()
        y = pmax(tree[[as.character(per)]]-K,0)
        otree[[as.character(per)]] = y
        for (i in per:1) {
          x = df*(q * y[-(i+1)] + omq * y[-1])
          otree[[as.character(i-1)]] = x
          y = x
        }
        otree
      }
      option.call.american = function(tree,ir,K,per){
        if (missing(per)) { per = length(tree)-1 }
        
        # compute risk neutral probabilities
        u = tree[[2]][1]/tree[[1]]
        d = tree[[2]][2]/tree[[1]]
        q = (1+ir-d)/(u-d)
        omq = 1-q
        df = 1/(1+ir)
        # print(c(q,df))
        otree = list()
        y = pmax(tree[[as.character(per)]]-K,0)
        otree[[as.character(per)]] = y
        for (i in per:1) {
          x = df*(q * y[-(i+1)] + omq * y[-1])
          x = pmax(x,tree[[ as.character(i-1)]] -K)
          otree[[as.character(i-1)]] = x
          y = x
        }
        otree
      }
      option.put = function(tree,ir,K,per){
        if (missing(per)) { per = length(tree)-1 }
        
        # compute risk neutral probabilities
        u = tree[[2]][1]/tree[[1]]
        d = tree[[2]][2]/tree[[1]]
        q = (1+ir-d)/(u-d)
        omq = 1-q
        df = 1/(1+ir)
        #print(c(u,d,q,df))
        otree = list()
        y = pmax(K-tree[[as.character(per)]],0)
        otree[[as.character(per)]] = y
        for (i in per:1) {
          x = df*(q * y[-(i+1)] + omq * y[-1])
          otree[[as.character(i-1)]] = x
          y = x
        }
        otree
      }
      option.put.american = function(tree,ir,K,per,div=div){
        if (missing(per)) { per = length(tree)-1 }
        
        # compute risk neutral probabilities
        u = tree[[2]][1]/tree[[1]]
        d = tree[[2]][2]/tree[[1]]
        q = (1+ir-div-d)/(u-d)
        omq = 1-q
        df = 1/(1+ir)
        # print(c(q,df))
        otree = list()
        y = pmax(K-tree[[as.character(per)]],0)
        otree[[as.character(per)]] = y
        for (i in per:1) {
          x = df*(q * y[-(i+1)] + omq * y[-1])
          x = pmax(x,K-tree[[ as.character(i-1) ]])
          otree[[as.character(i-1)]] = x
          y = x
        }
        otree
      }
      
      
      # Prices
      baum <- createBTreeSym(sig=sigma,dT=dT,per=per,S=S)
      baum.mat <- tree2matrix(baum)
      
      # Calloption
      o.c.e <- option.call(baum,ir=ir,K=K,per=per)
      o.c.e.mat <- tree2matrix(o.c.e,type="valuation") # Als Matrix
      bin.c.e <- baum.mat-K-o.c.e.mat
      bin.c.e[which(bin.c.e<0)] <- NA
      o.c.a <- option.call.american(baum,ir=ir,K=K,per=per)
      o.c.a.mat <- tree2matrix(o.c.a,type="valuation") #  Als Matrix
      bin.c.a <- baum.mat-K-o.c.a.mat
      bin.c.a[which(bin.c.a<0)] <- NA
      
      # Putoption
      o.p.e <- option.put(baum,ir=ir,K=K,per=per)
      o.p.e.mat <- tree2matrix(o.p.e,type="valuation")
      bin.p.e <- K-baum.mat-o.p.e.mat
      bin.p.e[which(bin.p.e<0)] <- NA
      o.p.a <- option.put.american(baum,ir=ir,K=K,per=per,div=0)
      o.p.a.mat <- tree2matrix(o.p.a,type="valuation")
      bin.p.a <- K-baum.mat-o.p.a.mat
      bin.p.a[which(bin.p.a<0)] <- NA
      
      #Output
      return(list(baum.mat=baum.mat,o.p.e.mat=o.p.e.mat,bin.p.e=bin.p.e,o.p.a.mat=o.p.a.mat,bin.p.a=bin.p.a
                  ,o.c.e.mat=o.c.e.mat,bin.c.e=bin.c.e,o.c.a.mat=o.c.a.mat,bin.c.a=bin.c.a,K=K))
      
    }
    
    # S = So Preis der Aktie zum zeitpunkt 0
    # K = Ausübungspreis (Preis der mit der Option abgemacht wird)
    # sigma = Volatilitäd des (logarythmischen) Preises
    # R = Risikofreie Zins pro Jahr = 1+ir
    # ir = Zins auf Periode gerechnet
    # per = Perioden (anz Schritte)
    # dT = Periodenlänge im Jahr (1/12 = 1 Schritt ist 1 Monat)
    
    # Einlesen der Werte
    S <- 53
    K <- 50
    sigma <- 0.25
    per <- 5
    dT <- 1/12
    r <- 0.12
    ir <- r * dT
    
    
    prices <- options.prices(S=S,K=K,sigma=sigma,ir=ir,per=per,dT=dT)
    
    # Baum mit den Preisen der Basisaktie (S)
    prices$baum.mat
    
    # Optionspreise für Call und Put und amerikanische und europäische 
    prices$o.p.e.mat
    prices$o.p.a.mat
    prices$o.c.e.mat # bei Call Option ist amerikanische und europ?ische das Selbe !
    
    # Wie viel mehr verdient man beim Auslösen der Option, zu wenn man sie behält
    # 
    prices$bin.p.e
    prices$bin.p.a
    prices$bin.c.e
    prices$bin.c.a
    
    # Wo lohnt es sich?
    (K - prices$baum.mat - prices$o.p.a.mat)==prices$bin.p.a
    
  }
  # BinomialTreeOption???
  function(){
    # basic concept created by Diethelm Wuertz (of Package: fOptions)
    # customized by Spasojevic Aleksandar (IDP)
    BinomialTreeOption = function (TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, sigma, n){
      # TODO: Exception-Handling
      
      # loading/installing package "fOptions"
      if(!is.element("package:fOptions", search())){
        if(is.element("fOptions",installed.packages())){
          library("fOptions")
        }else{
          install.packages("fOptions")
        }
      } 
      
      # Parameter
      dt = Time/n
      u = exp(sigma*sqrt(dt))
      d = 1/u
      q = (1+(r*dt)-d)/(u-d)
      p = 0.5*(1+(r/sigma)*sqrt(dt))
      Df = 1/(1+r*dt)
      
      # Declare if Price-Dynamic or OptionValue-Dynamic
      if(missing(TypeFlag)){
        Value = S*u^(0:n)*d^(n:0)
      }else{
        TypeFlag = TypeFlag[1]
        if (TypeFlag == "ce" || TypeFlag == "ca") 
          z = +1
        if (TypeFlag == "pe" || TypeFlag == "pa") 
          z = -1
        Value = z*(S*u^(0:n)*d^(n:0)-X)
      }
      
      # Calculating the Tree
      offset = 1
      Tree = Value = (abs(Value) + Value)/2
      if (TypeFlag == "ce" || TypeFlag == "pe") {
        for (j in (n - 1):0) {
          Tree <- c(Tree, rep(0, times = n-j))
          for (i in 0:j) {
            Value[i + offset] = (q*Value[i+1+offset]+
                                   (1-q)*Value[i+offset])*
              Df
            Tree = c(Tree, Value[i+offset])
          }
        }
      }
      if (TypeFlag == "ca" || TypeFlag == "pa") {
        for (j in (n-1):0) {
          Tree <- c(Tree, rep(0, times = n-j))
          for (i in 0:j) {
            Value[i+offset] = max((z*(S*u^i*d^(abs(i-j))-X)), 
                                  (q*Value[i+1+offset]+(1-q)*Value[i+offset])* 
                                    Df)
            Tree = c(Tree, Value[i + offset])
          }
        }
      }
      Tree = matrix(rev(Tree), byrow = FALSE, ncol = n + 1)
      invisible(Tree)
    }
    
    # Tree-Plot / von library(fOptions)
    # BinomialTreePlot(BinomialTreeOptionObject)
  }
  
  # TREE Short-Rates und Kassazinsen (vorher tree erstellen)
  function(){
    # Skript Beispiel
    tree <- createBTree(u = 1.3,d = 0.9,per = 5,S = 0.07)
    (tree.mat <- tree2matrix(tree))
    
    ### Baum wird zurück gerechnet
    # Anfangs wird ein kleiner Baum zurückgerechnit danach ein immer grösserer
    # Bis der ganze angegebene Baum mit allen Möglichkeiten durchlaufen wurde
    # und von jedem der "Preis" d.h. der Abzinsfaktor berechnet wurde.
    q <- 0.5 #
    l = 6 # per+1
    P = numeric(l)
    A = list()
    for (n in 1:l)
    {
      V = rep(1,n+1)
      A[[n]] = mat.or.vec(n+1,n+1)
      A[[n]][,n+1] <- rep(1,n+1)
      for (t in (n-1):0)
      {
        V = (q*V[-1]+(1-q)*V[-(t+2)])/(1+tree.mat[1:(t+1),t+1])
        if(is.vector(A[[n]])) A[[n]] = V
        if(is.matrix(A[[n]])) A[[n]][1:(t+1),t+1] = V
      }
      P[n] = V[1]
    }
    A # A = Bondvalue Verlauf
    c(1,P) # P = Bond Prices (letzes Bondvalue in Matrix)
    (s = P^(-1/1:l) - 1) # s = "Kassazinsen" / Spot rates
    
    plot(100*s, xlab="Laufzeit", ylab = "Zinssatz in %",type="b")
    
  }
  # Ho-Lee modell
  function(){
    ### 1. Short-Rate-Gitter
    srLattice = function(a,b,n=length(a)){
      lat = matrix(rep(0,n*n),nrow=n,ncol=n)
      for (t in 0:(n-1)) {
        for (i in 0:t) {
          lat[i+1,t+1] = a[t+1] + b * i
        }
      }  
      lat
    }
    par.a = c(2.50, 3.50, 4.36, 5.10, 5.67, 6.15, 6.66, 7.02, 6.47, 7.68)
    par.b = 1
    srLat = srLattice(par.a,par.b)
    srLat
    dim(srLat)
    
    ### 2. Elementarpreise Elementarpreisgitter
    ep = function(srLat){
      d = dim(srLat)+1
      ep.lat = matrix(rep(0,prod(d)),nrow=d[1],ncol=d[2])
      ep.lat[1,1] = 1
      
      # Vorwärtsgleichung
      for (t in 1:1) {
        ep.lat[1,2] = 0.5 * 1/(1+0.01*srLat[1,t]) * ep.lat[1,1]
        ep.lat[2,2] = 0.5 * 1/(1+0.01*srLat[t,t]) * ep.lat[1,1]
      }
      for (t in 2:(d[1]-1)) {
        ep.lat[1,t+1] = 0.5 * 1/(1+0.01*srLat[1,t]) * ep.lat[1,t]
        ep.lat[2:t,t+1] = 0.5 * (
          1/(1+0.01*srLat[1:(t-1),t]) * ep.lat[1:(t-1),t] +
            1/(1+0.01*srLat[2:t,t]) * ep.lat[2:t,t]
        )
        ep.lat[t+1,t+1] = 0.5 * 1/(1+0.01*srLat[t,t]) * ep.lat[t,t]
        
      }
      ep.lat
    }
    epLat = ep(srLat)
    epLat
    
    ### 3. die Preise der Nullcouponanleihen (P beim Short-Rates Gitter)
    function("Unnötige Funktion"){
      ep2zcb = function(ep.lat){
        colSums(ep.lat)[-1]
      }
      zcb.prices = ep2zcb(epLat) 
    }
    zcb.prices = colSums(ep.lat)[-1]
    zcb.prices
    
    ### 4. Berechnung der implizierten Kassazinskurve (Zinsen von 0 bis t, s beim Short-Rates Gitter)
    zcb2irSpot = function(zcb){
      ir.Spot = (zcb^(-1/(1:length(zcb))) - 1)
    }
    irSpot.HL = zcb2irSpot(zcb.prices)
    irSpot.HL
    ## Bondpreise berechnen
    cf.bond = c(rep(5,9),105)
    sum((1+irSpot.HL)^(-(1:10)) * cf.bond)
    zcb.prices%*%cf.bond # das selbe wie oben
    
    ### 5. Quadratischer Fehler zwischen empirischer Kassazinskurve und "Ho-Lee-Kurve"
    irSpot.emp = c(7.67, 8.27, 8.81, 9.31, 9.75, 10.16, 10.52, 10.85, 11.15, 11.42, 11.67, 11.89, 12.09, 12.27) # Kassazins
    irSpot.emp
    mse.hoLee = function(par,irSpot){
      b = par[1]
      a = par[-1]
      srLat = srLattice(a,b)
      #  print(srLat)
      epLat = ep(srLat)
      zcb.prices = ep2zcb(epLat)
      irSpot.hl = zcb2irSpot(zcb.prices)
      sum((irSpot-irSpot.hl)^2)
    }
    par.init = c(1, irSpot.emp[1], diff(irSpot.emp))
    par.init
    mse.hoLee(par.init, irSpot.emp) # Quadratischer Fehler
    
    ### 6. Anpassung
    res = optim(par.init, mse.hoLee, irSpot=irSpot.emp)
    res
    res$par
  }
  
  
  ### Futures Prices (future.prices) ???
  function(){
    
    future.prices <- function(tree, q){
      mat <- matrix(data = NA,nrow = nrow(tree),ncol = ncol(tree))
      mat[,ncol(tree)] <- tree[,ncol(tree)]
      for(i in ncol(tree) : 2){
        for( j in nrow(tree):2 ){
          mat[j,i-1] <- q*mat[j-1,i]+(1-q)*mat[j,i]
        }
      }
      return(mat)  
    }
    
    future.prices(baum.mat,2/3)
  }
  
  # BS-Formeln
  function(){
    
    # Mat = Time T (5 Monate = 5/12)
    # t = Startzeitpunkt
    # ir = JAHRESZINS
    # verbose = TRUE für zusatzinformationen d1, d2, Nd1, Nd2
    
    # Call
    BS.call = function(S,t=0,K=S,Mat=1,sig,ir, verbose=T){
      #  print(c(S,t,K,Mat,sig,ir))
      t2Mat = Mat - t
      fct = sig*sqrt(t2Mat)
      d1 = ( log(S/K) + (ir+sig^2/2)*t2Mat )/fct
      d2 = d1 - fct
      if(verbose) {
        print(c("d1=",d1,"d2=",d2))
        print(c("N(d1)=",pnorm(d1),"N(d2)=",pnorm(d2)))
      }
      S*pnorm(d1)-K*exp(-ir*t2Mat)*pnorm(d2)
    }
    
    BS.call.Delta = function(S,t=0,K=S,Mat=1,sig,ir, verbose=T){
      t2Mat = Mat - t
      fct = sig*sqrt(t2Mat)
      d1 = ( log(S/K) + (ir+sig^2/2)*t2Mat )/fct
      if(verbose) {
        print(c("d1=",d1))
        print(c("N(d1)=",pnorm(d1)))
      }
      pnorm(d1)
    }
    
    BS.call.Gamma = function(S,t=0,K=S,Mat=1,sig,ir, verbose=T){
      t2Mat = Mat - t
      fct = sig*sqrt(t2Mat)
      d1 = ( log(S/K) + (ir+sig^2/2)*t2Mat )/fct
      
      if(verbose) {
        print(c("d1=",d1))
        print(c("N(d1)=",pnorm(d1)))
      }
      exp(-0.5*d1^2) / (sig*S*sqrt(2*pi*Mat))
    }
    
    BS.call.Theta = function(S,t=0,K=S,Mat=1,sig,ir, verbose=T){
      t2Mat = Mat - t
      fct = sig*sqrt(t2Mat)
      d1 = ( log(S/K) + (ir+sig^2/2)*t2Mat )/fct
      d2 = d1 - fct
      if(verbose) {
        print(c("d1=",d1))
        print(c("N(d1)=",pnorm(d1)))
      }
      -ir*K*exp(-ir*Mat)*pnorm(d2) -
        sig*S*exp(-0.5*d1^2)/(2*sqrt(2*pi*Mat))
    }
    
    # Put
    BS.put = function(S,t=0,K=S,Mat=1,sig,ir,verbose=T){BS.put
      t2Mat = Mat - t
      fct = sig*sqrt(t2Mat)
      d1 = ( log(S/K) + (ir+sig^2/2)*t2Mat )/fct
      d2 = d1 - fct
      if(verbose) {
        print(c("d1=",d1,"d2=",d2))
        print(c("N(-d1)=",pnorm(-d1),"N(-d2)=",pnorm(-d2)))
      }
      -S*pnorm(-d1)+K*exp(-ir*t2Mat)*pnorm(-d2)
    }
    
    BS.put.Delta = function(S,t=0,K=S,Mat=1,sig,ir, verbose=T){
      t2Mat = Mat - t
      fct = sig*sqrt(t2Mat)
      d1 = ( log(S/K) + (ir+sig^2/2)*t2Mat )/fct
      if(verbose) {
        print(c("d1=",d1))
        print(c("N(d1)=",pnorm(d1)))
      }
      -pnorm(-d1)
    }
    
    BS.put.Gamma = function(S,t=0,K=S,Mat=1,sig,ir, verbose=T){
      t2Mat = Mat - t
      fct = sig*sqrt(t2Mat)
      d1 = ( log(S/K) + (ir+sig^2/2)*t2Mat )/fct
      if(verbose) {
        print(c("d1=",d1))
        print(c("N(d1)=",pnorm(d1)))
      }
      exp(-0.5*d1^2)/(sig*S*sqrt(2*pi*Mat))
    }
    
    BS.put.Theta = function(S,t=0,K=S,Mat=1,sig,ir, verbose=T){
      t2Mat = Mat - t
      fct = sig*sqrt(t2Mat)
      d1 = ( log(S/K) + (ir+sig^2/2)*t2Mat )/fct
      d2 = d1 - fct
      if(verbose) {
        print(c("d1=",d1))
        print(c("N(d1)=",pnorm(d1)))
      }
      ir*K*exp(-ir*Mat)*pnorm(-d2) -
        sig*S*exp(-0.5*d1^2)/(2*sqrt(2*pi*Mat))
    }
    
    
    # Digital
    BS.call.Digital = function(S, t = 0, K = S, Mat = 1, sig, ir, verbose = T) {
      # print(c(S,t,K,Mat,sig,ir))
      t2Mat = Mat - t
      fct = sig * sqrt(t2Mat)
      d1 = (log(S/K) + (ir + sig^2/2) * t2Mat)/fct
      d2 = d1 - fct
      if (verbose) {
        print(c("d1=", d1, "d2=", d2))
        print(c("N(d1)=", pnorm(d1), "N(d2)=", pnorm(d2)))
      }
      exp(-ir * t2Mat) * pnorm(d2)
    }
    BS.put.Digital = function(S, t = 0, K = S, Mat = 1, sig, ir, verbose = T) {
      t2Mat = Mat - t
      fct = sig * sqrt(t2Mat)
      d1 = (log(S/K) + (ir + sig^2/2) * t2Mat)/fct
      d2 = d1 - fct
      if (verbose) {
        print(c("d1=", d1, "d2=", d2))
        print(c("N(d1)=", pnorm(d1), "N(d2)=", pnorm(d2)))
      }
      exp(-ir * t2Mat) * pnorm(-d2)
    }
    
    
    ### vom Optionspreis implizierte Volatilität des Basisinstruments
    # allgemeine Werte werden verwendet in Funktion! müssen auch definiert werden.
    # S0 = 
    # K =
    # r = 
    # dT = 
    # C = vorher berechnet
    S0 <- 36.12;K <- 35;r <- 0.07;dT <- 7/52;sigma <- 0.3;C <- 2.4
    BSsigma <- function(sigma, C){
      integrand <- function(x) {1 / sqrt(2 * pi) * exp(-x ^ 2 / 2)}
      d1  <- (log(S0 / K) + (r + sigma ^ 2 / 2) * dT) / (sigma * sqrt(dT))
      d2  <- (log(S0 / K) + (r - sigma ^ 2 / 2) * dT) / (sigma * sqrt(dT))
      Nd1 <- integrate(integrand, lower = -Inf, upper = d1)$value
      Nd2 <- integrate(integrand, lower = -Inf, upper = d2)$value
      C1  <- S0 * Nd1 - K * exp(-r * dT) * Nd2
      C - C1
    }
    
    # Price und sigma berechnet von vorher
    price <- 0.24; sigma <- 0.3
    uniroot(BSsigma, c(0, 1), C = price)$root
    
  }
  # Eigenschaften von BS berechnen
  # auch mit BS-Formeln berechenbar
  function(){
    # über alles gültig
    K     <- 35
    price <- 2.15
    r     <- 0.07
    S0    <- 36.12
    sigma <- 0.05
    dT    <- 7/52 # T // Fälligkeitsdatum
    Nd1   <- NULL
    Nd2   <- NULL
    
    # d berechnen
    d1    <- (log(S0 / K) + (r + sigma ^ 2 / 2) * dT) / (sigma * sqrt(dT))
    d2    <- (log(S0 / K) + (r - sigma ^ 2 / 2) * dT) / (sigma * sqrt(dT))
    
    # merhere sigma ablaufen
    function(){
      sigma <- seq(0.01, 0.3, 0.01)
      integrand <- function(x) {1 / sqrt(2 * pi) * exp(-x ^ 2 / 2)}
      for(i in 1 : length(sigma))
      {
        Nd1[i] <- integrate(integrand, lower = -Inf, upper = d1[i])$value
        Nd2[i] <- integrate(integrand, lower = -Inf, upper = d2[i])$value
      }
      (C <- S0 * Nd1 - K * exp(-r * dT) * Nd2)
    }
    
    ## Call Option
    integrand <- function(x) {1 / sqrt(2 * pi) * exp(-x ^ 2 / 2)}
    (CNd1 <- integrate(integrand, lower = -Inf, upper = d1)$value)
    (CNd2 <- integrate(integrand, lower = -Inf, upper = d2)$value)
    ## Put Option
    integrand <- function(x) {1 / sqrt(2 * pi) * exp(-x ^ 2 / 2)}
    (PNd1 <- integrate(integrand, lower = -Inf, upper = -d1)$value)
    (PNd2 <- integrate(integrand, lower = -Inf, upper = -d2)$value)
    
    # BS-Formel anwenden
    (C <- S0 * CNd1 - K * exp(-r * dT) * CNd2)
    (P <- -S0 * PNd1 + K * exp(-r * dT) * PNd2)
  }
  
  # Monte Carlo Methode
  function(){
    
    # Durchläufe erstellen
    # Spalten = Realisationen
    # Zeilen = S(t) anz Schritte
    GBM1 =
      # Log-normalverteilt
      # Geometric Brownian motion
      # Algorithmus gemäss Gl. (22) in Kapitel 11.
      function(len, mu=sig^2/2, sig=1, dT=1, x0=1, nPath=1)
      {
        # 1. Term: nu * Delta t
        drift = (mu-sig^2/2)*dT
        # 2. Term: Faktor sigma * sqrt(Delta t)
        # (Standardabweichung von Zufallszahl)
        diff = sig*sqrt(dT)
        # Matrix mit Anz. Zeilen "len" und Anz. Spalten "nPath"
        # (Jede Spalte wird einen simulierten Pfad enthalten)
        # In jedes Matrixelement wird eine standard-normalverteilte 
        # Zufallszahl geschrieben (eps)
        res = matrix(rnorm(len*nPath), nrow=len, ncol=nPath)
        # Jedes Matrixelement wird so transformiert, dass sich 
        # der Term in der geschweiften Klammer in Gl. (22) ergibt.
        res = drift+diff*res
        # Am unteren Rand der Matrix wird eine Zeile hinzugefügt,
        # die überall den log. Anfangswert enthält
        res = rbind(rep(log(x0),nPath), res)
        # Die Werte werden spaltenweise aufsummiert und am Schluss wird
        # die ganze Matrix exponenziert.
        res = apply(res, 2, cumsum)
        exp(res)
      }
    # len:      Anzahl Zeitschritte in Pfad (länge)
    # mu:       Momentane Rendite
    # sig:      Volatilität
    # dT:       Länge des Zeitschritts
    # x0:       Anfangswert
    # nPath:    Anzahl Pfade
    dT <- 1/52    # dt = Zeitschritte (Bsp. Wochen)
    steps <- 21   # steps = Anazahl Zeitschritte (Bsp. 21)
    # Zeilen = S(t) anz Schritte, Spalten = Realisationen
    S.end = GBM1(steps, mu=0.12, sig=0.25, dT=dT,x0=53,nPath=10000)
    hist(S.end[(steps+1),])
    
    GBM2 =
      # Normalverteilt
      # Geometric Brownian motion,
      # Algorithmus gemäss Gl. (21) in Kapitel 11.
      function(len, mu=sig^2/2, sig=1, dT=1, x0=1, nPath=1)
      {
        # deterministischer Term in Gl (21)
        determ = 1+mu*dT
        #print(f1)
        # Faktor für stochastischen Term sigma * sqrt(Delta t)
        diff = sig*sqrt(dT)
        
        # Matrix mit Anz. Zeilen "len" und Anz. Spalten "nPath"
        # (Jede Spalte wird einen simulierten Pfad enthalten)
        # In jedes Matrixelement wird eine standard-normalverteilte 
        # Zufallszahl geschrieben (eps)
        res = matrix(rnorm(len*nPath), nrow=len, ncol=nPath)
        # Jedes Matrixelement wird so transformiert, dass sich 
        # der Term in der eckigen Klammer in Gl. (21) ergibt.
        res = determ + diff*res
        # Am unteren Rand der Matrix wird eine Zeile hinzugefügt,
        # die überall den log. Anfangswert enthält
        res = rbind(x0,res)
        # Die Werte werden spaltenweise aufmultipliziert
        res = apply(res, 2, cumprod)
      }
    
    option.call.mc =
      # Calloption mit Monte Carlo Methode berechnet
      function(n, S0, K, Tend, sig, r, dT, process=GBM1)
      {
        pLen = floor(Tend/dT)
        dFact = exp(-r*Tend)
        # Simulation der risikoneutralen Version der GBB:
        paths = process(pLen, mu=r, sig=sig, dT=dT, x0=S0, nPath=n)
        #print( dim(paths) )
        # Berechnung der risikoneutralen, abgezinsten erwarteten Auszahlung
        C = dFact * pmax( paths[pLen+1,]-K,0)
        #  c( mean(C), sqrt( var(C)/n ) )
        list( mean=mean(C), fluct=sqrt(var(C)/n), values=C )
      }
    # n:        Anzahl durchlaufene Realisationen
    # dT:       Zeitschrittgrösse (in Wochen => 1/52)
    # Tend:     anzahl Schritte in Zeitschritten (5 Wochen => 5/52)!!!
    # process:  GMB1 oder GMB2? Auswahl Namen eingeben
    opt.res = option.call.mc(n = 2000000,S0=53,K=50,Tend=5/12,sig=0.25,r=0.12, dT=1/52)
    c(opt.res$mean, opt.res$fluct)
    
    # Varianz berechnen und prozentzahl zu der Messung berechnen.
    # 10 mal besser => 100 mal mehr pfade     verbesserung^2 = multipliziert pfade
  }
  
  
}







