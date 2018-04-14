

function("Diskrete Verteilungsfamilien"){
  help(dbern)
  
  # Binominal Verteilung
  function(){
    dbinom(x, size, prob)
    pbinom(q, size, prob)
    qbinom(p, size, prob)
    rbinom(n, size, prob)
    
    prob = wahrscheinlichkeit eines einzelnen Ereignisses
    size = wie viele Versuche gibt es
    
    q = wie viele Erfogle soll es geben
    p = wie wahrscheinlich soll ein erfolg sein
  }
  
  # Poisson Verteilung
  function(){
    dpois(x = , lambda = )
    ppois(q = , lambda = )
    qpois(p = , lambda = )
    rpois(n = , lambda = )
    
    Lambda = durchschnittlicher Wert (vermutlich)
    
    q = wahrscheinlichkeit das das Ereigniss häufger als "q"mal vorkomt
  }
  
  # Geometrische Verteilung
  function(){
    dgeom(x = , prob = )
    pgeom(q = , prob = )
    qgeom(p =, prob = )
    rgeom(n = , prob = )
    
    prob = wahrscheinlichkeit bis zum Erfolg
    
    q = wie viele Misserfolge bis zum Erfolg
  }
  
  # Negative Binominal Verteilung
  function(){
    dnbinom(x = , size = ,prob = )
    pnbinom(q = , size = ,prob = )
    qnbinom(p = , size = ,prob = ) #gibt einem anzahl Fehlversuche, Anzahl Erfolge müssen addiert werden
    rnbinom(n = , size = ,prob = )
    
    prob = wahrscheinlichkeits eiens einzelnen Ereignisses
    size = wie viele Erfolge sind nötig
    
    q = Anzahl wie viele versuche es gibt
    p = wahrscheinlichkeit diese anzahl in Size zu erfüllen
  }
  
  # Hypergeometrische Verteilung
  function(){
    dhyper(x= ,m= ,n= ,k= )
    phyper(q= ,m= ,n= ,k= )
    qhyper(p= ,m= ,n= ,k= )
    rhyper(nn= ,m= ,n= ,k= )
    
    m = anz defekter                                     1. Stichprobe anz
    n = anz intakter                                     nicht gefundene Fehler anz
    k = Stichprobengrösse                                2. Stichprobe anz
    
    x = wie viele dürfen erwischt werden (exakte Zahl)   von beiden gefundene Bugs anz
    q = anz BIS wie viele erwischt werden dürfen!
      p = wahrscheinlichkeit das es erwischt wird
  }
}

function("Stetige Verteilungsfamilien"){
  
  # Exponentialverteilung
  function(){
    dexp(x= , rate= ) Dichte (keine W‘keit !!!)
    pexp(q= , rate= ) Verteilungsfunktion
    qexp(p= , rate= ) Quantilfunktion
    rexp(n= , rate= ) Simulation von Zufallszahlen
    
    hist(rexp(n= 100000, rate= 1))
    
    rate = "Lambda"
    
    Erklärung:
      Wir wollen nun den Fokus auf die Wartezeit zwischen zwei ankommenden Kunden richten.
    Da es sich hierbei um eine Zeit handelt,
    die mit beliebiger Präzision gemessen werden kann,
    ist es eine stetige Zufallsgrösse, die eine stetige Verteilung hat.
  }
  
  # Poisson-Verteilung
  function(){
    dpois()
    ppois()
    qpois()
    rpois()
  }
  
  # Stetige Gleichverteilung
  function(){
    dunif(x= , min = a, max = b) Dichte (keine W‘keit !!!)
    punif(q= , min = a, max = b) Verteilungsfunktion
    qunif(p= , min = a, max = b) Quantilfunktion
    runif(n= , min = a, max = b) Simulation von Zufallszahlen
    
    hist(runif(n=100000 , min = 0, max = 10))
    
    a = wie lange dauert es min
    b = wie lange dauert es max
    q = zb. zu wartende Min
    
    Erklärung:
      Die stetige Uniformverteilung(auch stetige Gleichverteilung) tritt dann auf,
    wenn alle Werte gleich oft eintreffen.
    Ein Beispiel für eine Zufallsvariable X mit dieser Verteilung ist die x-Koordinate des ersten Regentropfens,
    der bei einsetzendem Regen auf einem Gartentisch auftrifft.
    Sie kann Werte zwischen a=0 und der Länge des Tisches in Metern, z.B. b=2, annehmen.
  }
  
  # Normalverteilung
  function(){
    dnorm(x= , mean= , sd= ) Dichte (keine W‘keit !!!)
    pnorm(q= , mean= , sd= ) Verteilungsfunktion
    qnorm(p= , mean= , sd= ) Quantilfunktion
    rnorm(n= , mean= , sd= ) Simulation von Zufallszahlen
    
    hist(rnorm(n= 100000, mean= 0, sd= 10))
    
    X ~ N(a,b)
    Wurzel(b)= sd !!!!!!!!
      zusammenzähenl nur der b`s möglich
    a = mean
    
    Erklärung:
      Die Normalverteilung ist die wichtigste stetige Verteilung.
    Sie wurde als Modell für Messfehler entwickelt und passt für Messdaten meist sehr gut.
    Die Normalverteilung mit ihrer Glockenform ist also so etwas wie eine „Naturkonstante“ in der Welt der Zufallsvariablen.
  }
  
  # Multivariate Normalverteilung
  function(){
    dmvnorm(x, mean, sigma) Dichte
    pmvnorm(lower, upper, mean, sigma) Verteilungsfunktion
    qmvnorm(p, mean, sigma) Quantilfunktion
    rmvnorm(n, mean, sigma) Simulation von Zufallszahlen
  }
  
  # Chi-Quadrat-Verteilung (Symobl ist ein x^2)
  function(){
    Zusammenhang zwischen zwei Eigenschaften (wie häufig wurde welche erfüllt)
    aus der Matrix! => (anz Zeilen-1)*(anz Spalten-1) = Freiheitsgrad
    Formel für schätzwerte => sum(Spalte)*sum(Zeile)= Schätzung ein Wert
    Formel Teststatistik sum(sum((Beobachtet-Erwartet)/Erwartet))
    dchisq(x, df = Freiheitsgrade)
    pchisq
    qchisq
    rchisq
  }
  
  #wilcox
  function(){
    #falls nur wenige Beobachtungen zur Verfügung stehen
    dwilcox(x, m (anzahl 1. beobachtungen), n (anzahl 2. beobachtungen))
    pwilcox(p,m,n)
    qwilcox(q,m,n)
    rwilcox(nn,m,n)
  }
  
  # t-Verteilung Students-distribution
  function(){
    # df = Freiheitsgrade => Matrix (1.Dim - 1)*(2.Dim - 1)
    dt(x, df)
    pt(q, df)
    qt(p, df)
    rt(n, df)
  }
  
  # Weibul-Verteilung
  function(){
    dweibull()
    pweibull()
    qweibull()
    rweibull()
  }
  
  # Gamma
  function(){
    dgamma(x, shape, rate = 1, scale = 1/rate, log = FALSE)
    pgamma(q, shape, rate = 1, scale = 1/rate, lower.tail = TRUE,log.p = FALSE)
    qgamma(p, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
    rgamma(n, shape, rate = 1, scale = 1/rate)
  }
}




function("ordnen von Vektoren Reihenfolge"){
  q <- c(0,3,2,6,4,5)
  q
  sort(q,decreasing=T)     #Reihenfolge verdrehen
  sort(q)                  #nach Grösse ordnen
  
  q[order(q)]
  order(q) # An welcher stelle müsste die erste Zahl sein um sortiert zu sein. und ide 2. 3. usw. Zahl
  
  rank(q) # die wie vielt grösste Zahl ist an dieser Stelle?
  
  rev(q) #dreht die Zahlenreihe um
}

function("Vectorize"){
  f <- function(x,p,y) x*y+p
  
  x <- 1:5
  p <- 1:5
  y <- 1:5
  
  f(x,p,y)
  
  ?Vectorize
  f2 <- Vectorize(f,vectorize.arg="p")
  f2(x,p,y)
}

function("Schlaufen und Bedingungen"){
  # If else Bedingungen
  ifelse(Bedingung, wenn ja , wenn nein)
  
  # Bessere Schlaufe
  replicate(N,Function)
  
  # Schlaufe und Funktions
  for ( i in 1:try){
    a <- sample(1:365, k, replace = T)
    if(max(table(a))>1){
      counter <- counter + 1
    }
  }
  
  fak <- function(x){
    count <- 1
    for(i in 1:x){
      count <- count*x
      x <- x-i
      count
    }
  }
  fak(1)
}

function("Informationen aus der System lesen (Datum, Zeit)"){
  Sys.Date()
  Sys.time()
}

function("QQ-Plot"){
  install.packages("car")
  library(car)
  
  # QQ-Plot
  qqPlot(A, distribution="norm")
  # "norm" "exp" "unif" "t"
}

function("Vertrauensintervalle  Konfidenzintervalle Testen"){
    
    # T-Test  (Diverse Vektoren Testen)
    function(){
      t.test(x, y, conf.level = 0.95)
      t.test(rnorm(10))$p.value
      #vergelich von Werten (falls Paired T => kein var-equal!!!)
      t.test(x,y,alternative="greater",paired=F,var.equal=T)
    }
    
    # Z-Test wenn n<30 Messungen oder sigma nicht bekannt
    function(){
      library(BSDA)
      z.test(iq, sigma.x=15 , mu=100 ,conf.level=0.975)
    }
    
    # Wilcox-Test (falls zweifel an t-Test)
    function(){
      #  boxplot(Wert1, Wert2)          Verzogen?
      #  qqPlot(differenz der Werte)    Ausreisser?
      #falls nur wenige Beobachtungen zur Verfügung stehen
      #und falls Boxplots verzogen aussehen, nicht parametrisch
      wilcox.test(x, y, paired=T)
      # "greater" => x<y rechtsverschiebung     "less"=> x>y linksverschiebung
    }
    
    # Binominal-Test
    function(){
      binom.test(anz.Erfolge , anz.Tests , p=getestete Wahrscheinlichkeit)
      binom.test(502, 970, p=0.5, conf.level = 0.95) # Wert wird bis und mit 502 gerechnet
      
      poisson.test(x=27,T=3,r=referenz Anzhal, conf.level = 0.95, alternative = "less/greater")
      # 27 Störungen in 3 Monaten (entweder conf.level oder referenz Anzahl)
    }
    
    #Chisq-Test (Symobl ist ein x^2)
    function(){
      # Erwartete Häufigketit Tabelle (Vergleich von allen Werten zueinander)
      beob <- matrix(c(68,130,84,137,72,26), 2,3,byrow=TRUE)
      erw <- as.matrix(rowSums(beob)) %*% t(as.matrix(colSums(beob)))/sum(beob)
      erw
      
      chisq.test(beob)$expected
      chisq.test(beob) # nur Tabelle!  (Symobl ist ein x^2)
      chisq.test(beob)$residuals #(Abweichung zwischen beobachteten und erwarteten)
      
      # alternative = "two.sided", "greater" oder "less"
      # paired = TRUE/FALSE
    }
    
    #--------------------------------------------------------------------------------------
    # Paarweise Testen
    function(){
      pairwise.t.test(Messung, Messperson, p.adjust.method="none", pool.sd=F)
      pairwise.t.test(Messung, Messperson, p.adjust.method="bonf", pool.sd=F)
      # alle Daten in eine Zeile schreiben und mit einem 2. Vektor mit Bennenungen unterteilen
      # in Messung und Messpersonen. Bsp:
      Messung <- matrix(rnorm(n= 250, mean= 1, sd= 1),ncol=1)
      Messperson <- c(rep("V1",50),rep("V2",50),rep("V3",50),rep("V4",50),rep("V5",50))
      
      #zusammenzählen der werte die Bedingung erfüllen
      sum(pairwise.t.test(z,z1,p.adjust.method="bonf")$p.value<0.05,na.rm=T)
      # weitere Methode "none"
      
      tapply(Bereich, funktions Wert , funktion)
      tapply(isolier$Fehler, isolier$Hersteller,sum)
    }
  }

function("Transformationen"){
  Transformationen: Wurzel sqrt(), für Zähldaten => log(), für Beträge
  
  W <- lm(M[,1]~M[,2],data=M)
  summary(W)
  par(mfrow=c(1,6))
  plot(W, which=1:6)
  
  plot(M[,1]~M[,2],data=M)
  abline(W, col=2)
}






