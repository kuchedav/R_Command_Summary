
# DATEN VERWALTEN--------------STATISTIKEN--------------TABELLEN


#Konsole Befehle
rm(list = ls())
rm(Name)
ls()
setwd("~")
setwd("/Users/davidkuchelmeister/")
getwd()
load("Users/davidkuchelmeister/Downloads/xtsdata.Rda")
head(dir())
find.file("functions.R", dir = "Benchmark/")
list.files(pattern = "Portfolio")

# geht nur in einer Funktion (parentframe anzeigen)
dirname(parent.frame(2)$ofile)


### Nicht fertige Arbeiten
function("Verteilung anpassen???"){
  fitdistr(x, densfun, start, ...)
  
  library(fBasics) ## package loading
  skewness(x.norm) ## skewness of a normal distribution
  kurtosis(x.norm) ## kurtosis of a normal distribution
  skewness(x.wei) ## skewness of a Weibull distribution
  kurtosis(x.wei) ## kurtosis of a Weibull distribution
  
  library(MASS) ## loading package MASS
  x.gam <- rnorm(1000)
  truehist(x.gam)
  plot(table(round(x.gam,digits = 2)))
  fit <- fitdistr(x.gam,"t")
  plot(diff(pt(q = seq(-10,10,0.1),df = fit$estimate[3]),1),type="l")
  
  x <- diff(pt(q = seq(-10,10,0.1),df = fit$estimate[3]),1)
  lines(seq(-4,4,length.out = length(x)),x,col=2)
  
  plot(diff(pgamma(q = seq(-10,10,0.1),shape = fit$estimate[1],rate = fit$estimate[2]),1),type="l")
}
function("kaufen/verkaufen über Plot einzeichnen"){
  weight_one <- weights[,1]
  
  for(i in 1:length(weight_one)){
    change <- as.numeric(weight_one[1])-as.numeric(weight_one[i])
    if(change > 0.04) break
  }
  i
  for(j in i:length(weight_one)){
    change <- as.numeric(weight_one[i])-as.numeric(weight_one[j])
    if(change < -0.03) break
  }
  j
  
  
  plot(NULL,ylim=c(-1,1),xlim=c(0,length(weight_one)),
       yaxt="n",ylab="",xlab="",xaxt="n",frame=F)
  rect(0, 0, i, 1,col = "blue",density = 30)
  rect(i, -1, j, 0,col = "red",density = 30)
  par(new=T)
  plot.xts(weight_one)
  
  
  plot(sign(weight_one))
  as.numeric(weight_one>=0)
  
  
  plot(NULL,ylim=c(-1,1),xlim=c(0,length(weight_one)),
       yaxt="n",ylab="",xlab="",xaxt="n",frame=F)
  rect(0, 0, i, 1,col = "blue",density = 30)
  rect(i, -1, j, 0,col = "red",density = 30)
  par(new=T)
  plot.xts(weight_one)
  
  
  plot(sign(weight_one),col="lightskyblue2",type="h")
  plot(diff(sign(weight_one)))
  par(new=T)
  barplot(sign(weight_one),col="red",ylim=c(-1,1),border=NA,axes=F,axisnames=F)
  plot.xts(weight_one)
  colors()
  rgb()
}


function("Interessantes"){
  
  # Array
  array(dim=c(2,2,2))
  
  
  any() # ist eines der Statements TRUE?
  any(c(F,F,F,T,F,F))
  
  all.equal(1:10,1:10) # Wie 1:10==1:10
  
  # nimmt x als ausgabe element in Funktion, gibt aber kein Putin der Konsole.
  invisible(x)
  
  # make a text into a function
  eval(parse(text=paste0("hal","lo")))
  
  # erinnerung an which
  which(wg$Type == "turtles" & wg$Alive == "Y")
  # finden wo sich ein text in dem Vektor befindet
  grep("AA",names(DJ))
  match(c("hallo","ja","nein"),c("hallo","nein"),nomatch = 0) != 0
  
  
  install.packages("beepr")
  library(beepr)
  beep()
  
  
  # Funktionen
  ... # gibt Attribute weiter
  %in% # sucht ob Element in Objekt enthalten
    
    # Funktionen aus Package betrachten
    ls("package:QRM")
  ls() # Global Environment
  
  ### Diverses
  library(zoo)
  x <- c(1,NA,10,NA)
  x1 <- na.exclude(x)
  na.locf(x)
  na.omit(c(NA,1:10,NA,11:12))
  na.aggregate(c(NA,1:10,NA,11:12))
  plot(round(na.spline(c(NA,1:10,rep(NA,10),11:16)),2),type="l")
  
  ts.plot(na.spline(x))
  # der letzte wird widerholt
  ts.plot(na.locf(c(1:3,NA,5:8)))
  
  library(gmodels)
  CrossTable()
  
  ### Table plotten
  library(PerformanceAnalytics)
  PerformanceAnalytics:::textplot(halign = "center", valign = "center",object = T1)
  
  ### korrekt aus str() auslesen
  str(cmddata)
  attr(cmddata$points, "dimnames")[[1]] # f??r die zweite Stufe
  
  # optimierungs Aufgaben
  # funktion hat bsp. 2 Attribute a und b, die erstere wird ??berlaufen von 0 bis 1
  uniroot(funktion, a = c(0, 1), b = 1)
  
  ## Falls etwas fehlt in einer Funktion
  if (missing(per)){}
  
  # Return -1,0,1 ob die Zahl positiv negativ oder 0 ist.
  sign(-10:10)
  
  # attribut
  attr()
  attr(plo,"daves") <- c(1,2)
  attr(plo,"daves")
  
  # zahlenreihe umdrehen
  q <- c(0,3,2,6,4,5)
  rev(q)
  
  # Bessere Schlaufen
  vapply(X, FUN, FUN.VALUE, ..., USE.NAMES = TRUE)
  replicate(n, expr, simplify = "array")
  simplify2array(x, higher = TRUE)
  
  # HTML Matrix erstellen. (im Internet explorer anschaubar)
  HTML(dailyTests, file = paste(getwd(), "/dailyTests.htm", sep = ""),innerBorder = 1)
  
  # Args
  formalArgs(matrix)
  
  # Ausgabe der Version des R
  print(version)
  
  # Switch Command
  AA = 'foo'
  switch(AA, 
         foo={
           # case 'foo' here...
           print('foo')
         },
         bar={
           # case 'bar' here...
           print('bar')    
         },
         {
           print('default')
         }
  )
  
  
  ##### 2D Plots ####.
  persp(kde2d(runif(10000),runif(10000)),theta = -40, phi = 30, col = "steelblue")
  # erstellen eins Grids
  uv = grid2d(x = (0:25)/25) # grid für fCopulae
  plot(x=uv$x,y=uv$y)
  # draufsicht kontur plot
  contour(u, v, Verteilung, main="unabhängige Copula F")
}

function("Trading"){
  # Rolling functions
  library(RcppRoll)
  RcppRoll::
    
    # Returns
    (x[1:(length(x)-1)]-x[2:(length(x))])/x[2:(length(x))]
  # log-Returns
  getReturns(x)
  diff(log(x),differences=1)
  
  
  
  # Kreuz Korrelation
  ccf()
  ccf(x, y, lag.max = NULL, type = c("correlation", "covariance"), plot = TRUE, na.action = na.fail, ...)
  
  grep(pattern = ,x = ) # Strukturen in Zahlenreihen erkennen
  
  .hiddenfunction <- function(){cat("hidden function")} # wird nich angezeigt
}

function("Functions"){
  try( function...,silent = T) # code l??uft weiter, auch wenn funktion fehl schl??gt
  
  hallo <- function(){}
  eval(parse(text = "hallo")) # kann Text zu Command umwandeln
}

function("XTS"){
  xts(x = Zeitreihe,order.by = Datum der Zeitreihe) # von numeric
  as.xts() # NUR von einer Zeitreihe
  
  na.fill()
  na.exclude()
  
  xts::ndays()
  xts::nhours()
  xts::nquarters()
  
  xts::lag.xts()
  
  xts::dimnames.xts(Data)
  xts::first(Data); xts::last(Data)
  xts::plot.xts(na.exclude(Data[,1]))
  xts::shift.time()
  xts::to.quarterly(Data)
  xts::merge.xts()
  
  xts::to.monthly()
  
  xts.data["2010-01-01/2016-01-01"] # von bis Daten auslesen
  
  window(Data, window.start, window.end)
  
  # schöner XTS matplot
  matplot(return.mat,type="l",axes=F)
  axis(1, 1:nrow(return.mat), index(return.mat))
  axis(2)
}

function("Datum"){
  # umwandeln von Charaktern zu Datum
  strptime()
}

function("Knombinatorik"){
  ### Bin??res aufz??hlen (m??glichkeiten durchlaufen) Zahlen in Bin??r umwandeln
  library(sfsmisc)
  digitsBase(1:16)
  digitsBase(1:10,base = 2,ndigits = 2) # in base=2 system transorfmiereny
  as.intBase(digitsBase(100), base = 2) # zur??ck transformieren in 10er System von base=2
  # M??glichkeiten durchlaufen
  digitsBase(1:4,base = 2,ndigits = 2)
  
  # Kombinationen durchgehen
  (a <- expand.grid(1:10,1:10))
  
  # m Elemente werden aus Menge mit x Elementen gezogen
  combn(x = 10,m = 2,simplify = T)
}

function("neue Environments bauen"){
  # neu erstellen
  data <- new.env()
  
  # anschauen
  ls(data)
  class(data)
  
  a <- 100
  assign("a", 999, envir=data)
  ls(data)
  get('a', envir=data)
  a
  
  env.profile(data)
  environmentName(data)
  parent.env(data)
}

function("Text"){
  # ausschneiden
  strsplit("/tangency","/")[[1]][2]
  
  # ersetzen
  gsub("e", "", group)
  
  # doppelte löschen
  unique()
  
  # vereiningt die Werte des Vektors x mit den Werten des Vektors y,
  undion(x,y) # die nicht bereits in x enthalten sind. Umgekehrt funktioniert auch union(y, x)
  
  # Korrekturlesen
  library(hunspell)
  words <- c("Numbre", "Rate", "Year")
  hunspell_check(words)
}

function("Plots"){
  
  function("leerer Plot"){
    plot(NULL,xlim=c(-1,1),ylim=c(-1,1), yaxt="n",ylab="",xlab="",xaxt="n",frame=F)
  }
  
  function("Par / Layout"){
    # Par
    par(mfrow=c(1,1))
    par(mfrow=c(1,1),oma = c(1,1,1,1))
    par(mar=c(0,0,0,0)) # unten - links - oben - rechts
    # mar = Margin: freier Raum zwischen Rand und Grafik, die Angabe erfolgt z.B. in inches  
    # oma = Breite des a ??u??eren Rands in inches
    
    par(usr = c(0, 1, 0, 1)) # Alle Positionen können relativ zum Plot angegebwen werden
    par(new=T) # ein neuer Plot kan ÜBER den alten Plot gelegt werden
    
    par(mfrow=c(1,1),xaxs="i",yaxs="i",mar=c(4.5,4.5,2,1.5)+0.1)
    #Anzahl und Format der Abbildungen (mehrere Plots)
    # mar = Position und gr??sse des Plots im Fenster
    
    # Layout
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
    hist(rnorm(1000))
  }
  
  function("Details / Text"){
    #Unterschied der Symbole auf: http://de.wikibooks.org/wiki/GNU_R:_plot
    points(12, 5, pch = 3, cex = 2, lwd= 2)
    
    lines(density(wg$Miete), col="red")
    abline(v=10); abline(h=10)
    
    text(x1,x2,labels = "min. Var",pos = 4,col=2)
    
    # Plotten von mathematischen Symbolen
    ?plotmath
  }
  
  function("legenden für cex,lwd,usw..."){
    # cex = Skalierung der Symbolgr????e
    # lty = typ von Strich
    # lwd = Gr??sse der Striche
    # pch = Art der Punkte
    # pos = position
    # bty = ???
  }
  
  function("legende"){
    # (Plot fenster in gew??nschte h??he einstellen und plottden, damit Abst??nde korrekt sind)
    legend(0,27000,
           lty=c(1,1),
           lwd=1,col=1:ncol(Data),
           legend = colnames(Data),
           cex = 0.5,
           text.width = 50
    )
  }
  # Datenpunkt auslesen
  identify(x=Forbes$x, y=Forbes$y, cex= 0.7)
  
  function("altes zeugs"){
    # Pie Charts Kuchendiagramme
    pie(table(Erhebung$Augenfarbe)) # Erstellen eines Kuchediagramms mit Kategoriellen Werten
    
    # Stripchart (Boxplot, bei wenigen Werten)
    stripchart(dat, ylab = "dat", main = "horizontal") #mehrfach vorkommende Werte werden ??bereinander geplottet!
    
    # Mosaikplot (z.B. f??r zwei kategorielle Merkmale)
    mosaicplot(table(Erhebung$Augenfarbe,Erhebung$Nat??rliche.Harrfarbe))
    
    # alle Beziehungen
    pairs(data)
    
    #Histogramm Abbildung in der Konsole
    stem(data)
    
    # Viereck in Plot einzeichnen
    rect(xleft, ybottom, xright, ytop, density = NULL, angle = 45,
         col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"))
    rect(links, unten, rechts, oben,col = "red",density = 20)
  }
}

function("Ordner auf Pc verwalten"){
  # standart workin direktory
  setwd("~")
  setwd("/Users/davidkuchelmeister/")
  
  # current directory
  getwd()
  
  # which files are in this directory?
  head(dir())
}

function("Listen"){
  unlist(data)
  
  
  reduce(merge.xts, data)
}

function("Matrizen Rechnungen"){
  
  # Matrix mit Nullen (Dimensionen)
  mat.or.vec(10,5)
  
  install.packages("expm")
  library(expm)
  P%^%t
  
  t()                       # Transponieren
  upper.tri()  lower.tri()  # Untere Obere Dreicks Matrix
  %*%                       # Matrixmultiplikation
  det()                     # Determinante
  
  eigen(Mat)$values       # Eigenwerte
  abs(eigen(Mat)$values)
  
  # Die Matrix B wird mit jedem Element der Matrix A einzeln multipliziert
  # die daraus entstehenden Matrizen werden zu einer neuen grossen MAtrix angeordnet
  # in der selben ordnung wie A (siehe kronecker auf Wikipedia)
  kronecker(A,B)
  
  # Inverse einer nicht invertierbaren Matrix
  # (annäherung)
  inverse <- function(x) {
    # require(MASS)
    sing <- tryCatch(!is.matrix(solve(x)), error=function(e){ return(TRUE)}, finally={})
    if (sing==TRUE) {
      xinv <- MASS::ginv(x)
    } else {
      xinv <- solve(x)
    }
    return(xinv)
  }
}

function("Klassen Analyse"){
  class("hello")
  inherits("hello",what = "character")
}

function("Excel schreiben und lesen"){
  library(xlsx)
  Book <- createWorkbook()
  cars1 <- createSheet(wb=Book, sheetName="a")
  addDataFrame(x=a, sheet=cars1)
  saveWorkbook(data, "Studium/GitHub/Datenverwaltung/data.xlsx")
  
  read.xlsx(file = "Studium/GitHub/Datenverwaltung/data.xlsx",sheetIndex = 1)
}

function("Faster Code"){
  
  # user Matrixes instead of data frames. the simples the form the better
  
  
  # small calculations sometimes make a lot of error checking
  # turn this off with:
  .Internal(mean(x))
  
  
  
  # multicore calculation
  require(parallel)
  require(doMC)
  registerDoMC(cores=detectCores())
  level.premium.list <- foreach(i = 1:100) %dopar% chosen.function()
  
  
  # compiling functions
  library(compiler)
  chosen.function.new <- cmpfun(chosen.function)
  level.premium.list <- foreach(i = 1:100) %dopar% chosen.function.new
  
  
  
  
  
  # comaprison
  library(microbenchmark)
  compare <- microbenchmark(
    function1,
    function2,
    times = 20
  )
  
  
  
  function("finding Bottlenecks"){
    Rprof("out.out")
    
    level.premium.list <- foreach(i = Policy.Records$Policy.Number) %dopar% 
      optimise(f = level_of_premium.compiler,optimising.range,tol = tolerance, i)$minimum
    
    Rprof(NULL)
    summaryRprof("out.out")
    
    
    # different analysation
    proftable <- function(file, lines = 10) {
      profdata <- readLines(file)
      interval <- as.numeric(strsplit(profdata[1L], "=")[[1L]][2L]) / 1e+06
      filelines <- grep("#File", profdata)
      files <- profdata[filelines]
      profdata <- profdata[-c(1, filelines)]
      total.time <- interval * length(profdata)
      ncalls <- length(profdata)
      profdata <- gsub("\\\"| $", "", profdata)
      calls <- lapply(profdata, function(x) rev(unlist(strsplit(x, " "))))
      stacktable <- as.data.frame(table(sapply(calls, function(x) paste(x, collapse = " > "))) / ncalls * 100, stringsAsFactors = FALSE)
      stacktable <- stacktable[order(stacktable$Freq[], decreasing = TRUE), 2:1]
      colnames(stacktable) <- c("PctTime", "Call")
      stacktable <- head(stacktable, lines)
      shortcalls = strsplit(stacktable$Call, " > ")
      shortcalls.len <- range(sapply(shortcalls, length))
      parent.call <- unlist(lapply(seq(shortcalls.len[1]), function(i) Reduce(intersect, lapply(shortcalls,"[[", i))))
      shortcalls <- lapply(shortcalls, function(x) setdiff(x, parent.call))
      stacktable$Call = sapply(shortcalls, function(x) paste(x, collapse = " > "))
      if (length(parent.call) > 0) {
        parent.call <- paste(paste(parent.call, collapse = " > "), "> ...")
      } else {
        parent.call <- "None"
      }
      frac <- sum(stacktable$PctTime)
      attr(stacktable, "total.time") <- total.time
      attr(stacktable, "parent.call") <- parent.call
      attr(stacktable, "files") <- files
      attr(stacktable, "total.pct.time") <- frac
      print(stacktable, row.names=FALSE, right=FALSE, digits=3)
      if(length(files) > 0) {
        cat("\n")
        cat(paste(files, collapse="\n"))
        cat("\n")
      }
      cat(paste("\nParent Call:", parent.call))
      cat(paste("\n\nTotal Time:", total.time, "seconds\n"))
      cat(paste0("Percent of run time represented: ", format(frac, digits=3)), "%")
      
      invisible(stacktable)
    }
    proftable("out.out")
    
    
  }
}




