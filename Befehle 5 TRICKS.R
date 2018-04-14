# R-Commander ausf??hren
library(Rcmdr)

# OOP Programmieren (OESSY2/WMAR/Folien)
function(){
  # Eine Klasse wird erstellt
  # --------------------------------------------------------------------------------
  function(){
    
    # Eine Normale Klasse erstellen
    function(){
      # Es ben??tigt mind. 2 Angaben:
      #  Name der Klasse: "Bond"
      #  Informationen die gespeichert werden sollen: "representation"
      
      setClass(Class="Bond",
               representation=representation(issueDate = "character",
                                             maturityDate = "character",
                                             coupon = "numeric",
                                             nominal = "numeric"))
      
      
      
      
      
      # Die Klasse kann somit erstellt und einem Namen zugewiesen werden
      myBond <- new("Bond", issueDate="2013-01-01",
                    maturityDate="2023-01-01", coupon=5, nominal=100)
      
      
      # Werden nicht alle Angaben eingegeben, k??nnen sie im nachhinein eingetragen werden
      myBond2 <- new("Bond", issueDate="2013-01-01",
                    maturityDate="2023-01-01",nominal=100)
      myBond2@coupon <- 5
    }
  
    # Klasse mit eingebauten Grundwerten erstellen
    function(){
      # Die Klasse wird mit Grundwerten erstellt
      # d.h. wird die Klasse ohne vollst??ndige Angaben erstellt,
      # dann nehmen die prototypangeben ihren Platz ein
      setClass(Class="Bond",
               representation=representation(issueDate = "character",
                                             maturityDate = "character",
                                             coupon = "numeric",
                                             nominal = "numeric"
               ),
               prototype=prototype(issueDate = as.character(Sys.Date()),
                                   maturityDate = "2099-12-31",
                                   coupon = 0,
                                   nominal = 100
               )
      )
      
      # Die nicht eingegeben Daten werden automatisch angenommen
      myBond <- new("Bond", issueDate="2013-01-01",
                    maturityDate="2023-01-01", nominal=100)
     }
  
    # Klasse, welche Bedingungen ??berpr??ft erstellen
    function(){
  
      # Validity ??berpr??ft Bedingungnen
      
      setClass(Class="Bond",
               representation=representation(issueDate = "character",
                                             maturityDate = "character",
                                             coupon = "numeric",
                                             nominal = "numeric"),
               validity=function(object){
               if(object@coupon<0 || object@nominal<0)
                 stop("Coupon and Nominal should be positive !")
               }
      )
      new("Bond", issueDate="2013-01-01",
          maturityDate="2023-01-01", coupon=-5, nominal=100)
    
      
      #Validity nach der Erstellung der Klasse noch einf??gen
      function(){
          # Bedingung erstellen (Beispiel nicht von Bond Variablen k??nnen variieren)
          validTrackObject <- function(object) {
            if(length(object@x) == length(object@y)) TRUE
            else paste("Unequal x,y lengths: ", length(object@x), ", ",
                       length(object@y), sep="")
          }
          
          setValidity("Bond", validTrackObject)
        }
      }
  }
  
  # Erstellen einer Methode
  # --------------------------------------------------------------------------------
  function(){
    
    # Normale Methode erstellen (muss bereits existieren)
    function(){
      # (Die Methode muss bereits existieren!!! zB Plot, falls nicht siehe weiter unten)
      # eine Methode ist eine Funktion, welche der Klasse zugewiesen wird
      # Sie ben??titg mind. 3 Angaben
      #  f = "Name"
      #  signature = "zugeh??rige Klasse"
      #  definition= "Funktion die ausgef??hrt werden soll"
      
      setMethod(
                f="plot",
                signature="Bond",
                definition=function(x, y){
                n <- floor(as.numeric(as.Date(x@maturityDate) -
                                        as.Date(x@issueDate))/360)
                cfs <- rep(x@coupon, n)
                cfs[n] <- cfs[n] + x@nominal
                mat <- as.POSIXlt(as.Date(x@maturityDate))
                dates <- as.Date(paste(mat$mday, mat$mon+1,
                                       mat$year+1900-seq(0,n-1),
                                       sep="."), format="%d.%m.%Y")
                barplot(cfs, names.arg=sort(dates), main = "Cash Flows")
                }
      )
      
      # Ausf??hren der Methtode
      plot(myBond)
    }
  
    # Methode, welche Daten nur Anzeigt erstellen
    function(){
      # Abgek??rzte Version einer Mehtode
      # Gibt nur Informationen aus mit dem Befehl "cat"
      
      setMethod("print", "Bond",
                function(x){
                cat("********** Bond ***********")
                cat("\nIssue Date: ", x@issueDate)
                cat("\nMaturity Date: ", x@maturityDate)
                cat("\nCoupon: ", x@coupon)
                cat("\nNominal: ", x@nominal)
                cat("\n***************************")
                }
      )
      
      # Ausf??hren der Methtode
      print(myBond)
    }
  
    # Methode die noch nicht existiert erstellen
    function(){
      # Falls die Methode noch nicht existiert
      # muss sie zuerst erzeugt werden mit set Generic
      
      setGeneric(name="cashFlows",
                 def = function(object){
                 standardGeneric("cashFlows")
                 }
      )
      cashFlows(myBond)
      
      # Falls diese Funktion bereits exisiter wird sie einfach ??berschrieben!
      # Um ein ??berschrieben einer erstellten Generischen funktion zu verhindern wird dies benutzt:
      # (vorsichtig sein mit dieser Funktion)
      
      lockBinding("cashFlows", .GlobalEnv )
      
      
      
      # Nicht ganz sicher f??r was das genau ist, aber man braucht es einfach xD
      library(timeSeries, verbose=FALSE)
      
      # festlegen einer Methode, was passiert, wenn sie aufgerufen wird
      # Diese Methode existierte zuvor noch nicht und musste zuerst mit setGeneric erzeugt werden
      setMethod(f="cashFlows",
                signature="Bond",
                definition=function(object){
                n <- floor(as.numeric(as.Date(object@maturityDate) -
                                        as.Date(object@issueDate))/360)
                cfs <- rep(object@coupon, n)
                cfs[n] <- cfs[n] + object@nominal
                mat <- as.POSIXlt(as.Date(object@maturityDate))
                dates <- as.Date(paste(mat$mday, mat$mon+1,
                                       mat$year+1900-seq(0,n-1),
                                       sep="."), format="%d.%m.%Y")
                cash.flows <- timeSeries(cfs, charvec=sort(dates))
                colnames(cash.flows) <- "Cash_Flows"
                return(cash.flows)
                }
      )
      
      # Ausf??hren der eigenen Methode
      cashFlows(myBond)
    }
  
    # Methode zum der Klasse neue Werte zuweisen <-
    function(){
      setGeneric("setCoupon<-", function(object, value) { standardGeneric("setCoupon<-") }
                 
      setReplaceMethod(f="setCoupon", signature="Bond",
                                 definition=function(object, value){
                                 object@coupon <- value
                                 return(object)
                                 }
                )}
    
      setCoupon(myBond) <- 100
    
    # Methode zur eingeben von Neuen Werten in die Klasse (nicht so wichtig)
      function(){
      # wird nur ausgef??hrt, wenn issueDate & maturityDate vorhanden
      
      setMethod(f="initialize",
                signature="Bond",
                definition=function(.Object, issueDate, maturityDate, coupon, nominal){
                  if (!missing(issueDate) && !missing(maturityDate)){
                    issueDate <- format(as.Date(issueDate), "%d.%m.%Y")
                    maturityDate <- format(as.Date(maturityDate), "%d.%m.%Y")
                    ## Slots Values
                    .Object@issueDate <- issueDate
                    .Object@maturityDate <- maturityDate
                    .Object@coupon <- coupon
                    .Object@nominal <- nominal
                    ## Return Object
                    validObject(.Object)
                   }
                return(.Object)
                }
      )
      
      initialize
      
      new("Bond", issueDate="2013-01-01", maturityDate="2023-01-01",
          coupon=-5, nominal=-100)
    }
  }
  
  # ??bersicht ??ber die Klasse erhalten! (Wichtig)
  # --------------------------------------------------------------------------------
  function(){
    
    getMethod()
    existsMethod()
    hasMethod()
    selectMethod()
    setMethod()
    
    # Father/Son OOP
    callNextMethod()
    print(as(myFRN2, "Bond"))
    
    
    # Zugriff auf einzelne Informationen aus der Klasse
    myBond@coupon
    
    # Ab??ndern der Informationen
    myBond@coupon <- 4 ## NOT RECOMMENDED
    
    # Eigenschaften der Klasse abrufen
    slotNames("Bond")
    getSlots("Bond")
    
    
    # Welche Methoden enth??lt "Bond"?
    showMethods(class="Bond", where = ".GlobalEnv")
    
    # Wie sieht diese Methode aus?
    getMethod(f="print", signature="Bond")
    
    # Existiert die methode "..." TRUE/FALSE
    existsMethod(f="summary", signature="Bond")
  }
  
  # Benutzerfreundliches Konstruieren
  # --------------------------------------------------------------------------------
  function(){
    # Erstellen der Klasse ??ber eine Funktion, mit Test zur Angabe wa passiert
    function(){
      
      # Erstellen der Funktion
      bond <- function(issueDate, maturityDate, coupon, nominal){
        
        # Bedingung, falls der Benutzer die maturity im falschen Format eingibt
        maturityDate <- ifelse(!is.numeric(maturity), as.character(maturity),
                               as.character(nYears(issueDate, maturity)))
        
      cat("**** BOND : CONSTRUCTOR ****nn")
      new(Class="Bond", issueDate=issueDate, maturityDate=maturityDate,
          coupon=coupon, nominal=nominal)
      }
      
      # Informationen k??nnen bequem in Funktion eingegeben werden
      bond(issueDate="2013-01-01", maturityDate="2023-01-01", coupon=5, nominal=100)
    }
  }

  # Begriffe
  # --------------------------------------------------------------------------------
  function(){
      redemptionValue = Nominaler Wert (meist 100)
      issue Date = beginn der Zeit
      Value Date = Zu diesem Zeitpunkt wird der Preis berechnet
      Maturity Date = Ende der Zeitspanne
      initial Date = 
      spot Date = ausgelesenes Datum aus Curve (meist mit Value Date gleichgesetzt)
      
      finCalcOessy:::plot.Curve(myCurve) #zwingen die Funktion aus finCalcOessy zu nehmen
    }
}

# Server verbindung mit R-Studio erstellen
function(){
  # Internetseite f??r Hilfe:
  # http://stackoverflow.com/questions/26210317/installation-of-rodbc-roracle-packages-on-os-x-mavericks
  # ben??tigt unixodbc, muss gegoogelt und instelliert werden
  # Im Dropdown Menu (Tools > Shell...) kann eine Konsole ge??ffnet werde, wo die befehle geschrieben werden
  
  install.packages("RODBC",type="source")
  library(RODBC)
  
  odbcDriverConnect("driver=SQL Server;server=01wh155073")
  
  
  # Anleitung von Kevin
  function(){
  1) Rmysql herunterladen!
  2) 1 mal verbinden
  3) functions einlesen
  4) functions bedienen
  5) Allgemein: SQL befehl werdet als String mit send query a de datebank geschickt ( FROM * SELECT table ~ sql)
  
  
  
  # Input--------------------------------------------------------------
  
  # Package herunterladen
  require(RMySQL) 
  library(RMySQL)
  library(DBI)
  #verbinden (1-mal,keine Fehlermeldung?,Warten bis "stop" Symbol verschwindet, i.o)
  con<-dbConnect(dbDriver("MySQL"),user='sql373333',password='dA8%qJ4*',host='sql3.freesqldatabase.com',dbname='sql373333')
  
  send("lol","troll",con) 
  chat<-get(con)
  dbDisconnect(con) #Verbindung trennen , nicht n??tig
  
  
  #--------------------------------------------------------------------
  
  
  #---functions-Nur einlesen-------------------------------------------
  send<-function(from,msg,con){
    a<-"INSERT INTO `table`(`indexx`, `from`, `Msg`) VALUES ('','"               
    comma<-"','"                     
    need<-"')"
    sendq<-paste(a,from,comma,msg,need)
    dbSendQuery(con, sendq)       
  }
  get<-function(con){
    res<-dbSendQuery(con, "SELECT * FROM `table`")
    allMsg<- fetch(res, n = -1)
    View(allMsg)
    return(allMsg)
    dbClearResult(res)
  }
  #---functions-Nur einlesen-------------------------------------------
  
  
  abc<-dbConnect(dbDriver("MySQL"),user='sql374253',password='xW6*vJ2!',host='sql3.freemysqlhosting.net',dbname='sql374253')
  
  send("lol",abc)
  dbDisconnect(con)
  
  
  
  
  
    }
}

# shiny App
function(){
  # Anmelden
  function(){
    # geheim
    function(){
      shinyapps::setAccountInfo(name='davidch12345',
                              token='0BB2CC7CD728C436D704D6E474485B97',
                              secret='nonXWdJ0mCl0nub75/3g8vEV+gUGRvEt9iaUiqJk')
    }
    
    library(shinyapps)
    # shinyapps::deployApp('path/to/your/app')
    
    install.packages('devtools')
    devtools::install_github('rstudio/shinyapps')
    library(shinyapps)
    install.packages(c('ggplot2', 'shiny'))
  }
  
  # f??r jedes App muss ein Ordner erstellt werden der die Datei ui.R und server.R enth??lt (exact diese Namen)
  # am besten Tutorials auf Youtube anschauen
  
  setwd("/Users/davidkuchelmeister/Desktop/Alt Desktop/Studium/2. Jahr/4. Semester/")
  library(shiny)
  runApp("Shiny/App1/")  # funktioniert ??hnlich wie read.table oder so...
  
  library(shinyapps)
  deployApp("Shiny/App1/")
  
  # Tutorial
  http://shiny.rstudio.com/tutorial/lesson1/
  
  
  # Problembehandlung
  shinyapps::showLogs("Shiny/App1/")
  terminateApp("Shiny/App1/") # App sofort beenden
  
  # Beispiele
  function(){
    runExample("01_hello") # a histogram
    runExample("02_text") # tables and data frames
    runExample("03_reactivity") # a reactive expression
    runExample("04_mpg") # global variables
    runExample("05_sliders") # slider bars
    runExample("06_tabsets") # tabbed panels
    runExample("07_widgets") # help text and submit buttons
    runExample("08_html") # Shiny app built from HTML
    runExample("09_upload") # file upload wizard
    runExample("10_download") # file download wizard
    runExample("11_timer") # an automated timer
  }
}
