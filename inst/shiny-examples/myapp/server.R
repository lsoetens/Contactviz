library(shiny)
#load packages
library(xlsx)
library(igraph)
library(ggplot2)
library(plyr)
library(grid)
library(scales)
library(gtable)
library(gridExtra)
library(RGraphics)
library(binom)
library(Hmisc)
library(boot)
library(RColorBrewer)

#read data
#setwd("R:/EPI/MOD/Projects/EW_GIS_CIb_Datavisualisatie/Work-data_scripts_programs/5_Complexe_time_type_place_person_data/Contact_tracing_information/")


# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
        
        output$plotcontacttrace <- renderPlot({
        
          
        inFile <- input$file1
                
        if (is.null(inFile))
        return(NULL)
        
                
        data<- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                        stringsAsFactors=F)
                        
        #set the current date
        currentdate<- reactive({ as.Date(input$currentdate, origin = "1970-01-01")})
        currentdate<- currentdate()
        #set dates for annotation
        #Outbreakdetected
        outbreakdetect<- reactive({ as.Date(input$date1, origin = "1970-01-01")})
        outbreakdetect<- outbreakdetect()
        date1label<- reactive({input$date1label})
        date1label<- date1label()
        #start control measures
        controlstart<- reactive({ as.Date(input$date2, origin = "1970-01-01")})
        controlstart<- controlstart()
        date2label<- reactive({input$date2label})
        date2label<- date2label()
        #set mean generation time in days
        generationtime<- reactive({ as.numeric(input$generationinterval)})
        generationtime<- generationtime()
        #set mean incubation period in days
        incubationperiod<- reactive({as.numeric(input$incubationperiod)})
        incubationperiod<- incubationperiod()
        #set sd incubation period in days
        incubationsd<- reactive({as.numeric(input$sdincubationperiod)})
        incubationsd<- incubationsd()
        
        
        #country and disease of outbreak for title plot
        country<- reactive({input$country})
        country<- country()
        disease<- reactive({input$disease})
        disease<- disease()
        
        #width of surveillance interval
        surveillance<- reactive({input$surveillance})
        surveillance<- surveillance()
        surveillance<- as.numeric(surveillance)
        problow<- ((100-surveillance)/2)/100
        probhigh<- 1-(((100-surveillance)/2)/100)
        
        
        #create some new variables and do some formatting
        data<- within (data, {
                date_onset_disease<- as.Date(DOD, format="%Y-%m-%d")
                date_exposure<- as.Date(Exposure_date, format="%Y-%m-%d")
                Gender<- factor(Gender, levels= c("man", "woman", "unknown"), ordered=TRUE)
                Exposure_type<- factor(Exposure_type, levels=c("Family/Household", "Friend/neighbor", "Medical/Religious","Work", "Bar/Shop", ""))
        }
        )
        
        
        #Remove cases and their contacts if dod of case is after currentdate and exposure type is unknown
        idselect<- data$ID[data$date_onset_disease> currentdate & data$Exposure_type==""]
        data<- subset(data, !(data$ID %in% idselect) & !(data$IDSource %in% idselect))
        
        #set cases diagnosed after the currentdate as controls
        idselect<- data$ID[data$Type=="Case" & data$date_onset_disease > currentdate]
        data<- subset(data, !(data$IDSource %in% idselect))
        data$Type<- ifelse(data$ID %in% idselect, "Contact", data$Type)
        
        
        #########################
        #calculate Feff for every contact
        
        #first select contacts
        data2<- subset(data, (data$Type== "Contact")) 
        
        
        #create function for creating convolution of exposure interval and incubation period distribution
        convolution<- function(expdate, expduration, nsim){
                set.seed(200)
                x<- expdate + runif(nsim, min=0, max= expduration)+ rlnorm(nsim, log(incubationperiod), log(incubationsd))
                return(x)
        }
        
        #returns a matrix, with in every column the simulated dates for one contact
        simdates<- mapply(convolution, expdate= data2[,"date_exposure"], expduration= data2[,"Exposure_duration"], nsim=10000 )
        
        
        #apply ecdf function to every column in simdates
        cumsim<- apply(simdates, 2, ecdf)
        
        
        #calculate cum. probability (and 1- prob) for certain date for every column
        arglist<- as.list(c(rep(as.numeric(currentdate), length(data2[, "ID"]))))
        prob<- mapply(function(x,y) x(y), x= cumsim, y=arglist) #Feff
        probsymp<- round(1-prob,6) #Seff
        date05<- as.Date(round(as.numeric(apply(simdates, 2, quantile, prob=problow))),origin = "1970-01-01")
        date95<- as.Date(round(as.numeric(apply(simdates, 2, quantile, prob=probhigh))),origin = "1970-01-01")
        
        
        #attach the calculated probabilities to original dataframe 
        data3<- data.frame(ID=data2$ID, prob= prob, probsymp=probsymp, date05= date05, date95= date95)
        data<- merge(data, data3, by="ID", all=T)
        
        
        # make a graphdate, which corresponds to the 1st day of symptom onset for cases,
        # and to the lower 5% of the symptom onset interval for contacts
        data<- within (data, {
          graphdate<- ifelse(Type=="Case", date_onset_disease, date05)
          graphdate<- as.Date(graphdate, origin = "1970-01-01")
        })
        
        
        
        #########################################
        #MLE for probability on infection for every contact type 
        
        LE<- function(ns, na, Feff, pi){
          x<- ns*log10(pi)+ na*log10(1-pi) + sum(log10(1-pi*Feff))
          return(x)
        }
        
        AR<- c()
        ARlo<- c()
        ARhi<- c()
        Exposure_type<- c()
        ntotal<- c()
        for(i in 1: length(unique(data$Exposure_type))){
          test<- subset(data, Exposure_type== unique(data$Exposure_type)[i])
          Exposure_type<- c(append(Exposure_type, as.character(unique(test$Exposure_type)), after=length(Exposure_type)))
          ntotal<- c(append(ntotal, length(test$Exposure_type), after=length(ntotal)))
          if (length(test$Exposure_type)> 5){
            
            ns<- length(test$Type[test$Type== "Case"]) #number of cases
            na<- length(test$Type[test$Type=="Contact" & test$probsymp< 0.05]) #number of contacts no symptoms in surveillance period
            pi<- c()
            ll<- c()
            for (j in 1:1000){
              PI<- j/1000
              LL<- LE(ns=ns, na=na, Feff=test$prob[test$Type=="Contact" & test$probsymp> 0.049], pi=PI)
              pi<- c(append(pi, PI, after=length(pi)))
              ll<- c(append(ll, LL, after=length(ll)))
            }
            MLE.df<- data.frame(pi= pi, ll=ll)
            maxll<- MLE.df$pi[MLE.df$ll == max(MLE.df$ll)] #maxLL
            MLE.df$dev<- 2*(max(MLE.df$ll)-MLE.df$ll)
            lllow<- MLE.df$pi[which(MLE.df$dev < 3.841)[1]] #lower CI MAXLL
            llhi<- MLE.df$pi[which(MLE.df$dev < 3.841)[length(which(MLE.df$dev < 3.841))]] #higher CI MaxLL
            AR<- c(append(AR, maxll, after=length(AR)))
            ARlo<- c(append(ARlo, lllow, after=length(ARlo)))
            ARhi<- c(append(ARhi, llhi, after=length(ARhi)))
          }
          else {
            
            ns<- length(data$Type[data$Type== "Case"]) #number of cases
            na<- length(data$Type[data$Type=="Contact" & data$probsymp< 0.05]) #number of contacts no symptoms in surveillance period
            pi<- c()
            ll<- c()
            for (j in 1:1000){
              PI<- j/1000
              LL<- LE(ns=ns, na=na, Feff=data$prob[data$Type=="Contact" & data$probsymp> 0.049], pi=PI)
              pi<- c(append(pi, PI, after=length(pi)))
              ll<- c(append(ll, LL, after=length(ll)))
            }
            MLE.df<- data.frame(pi= pi, ll=ll)
            maxll<- MLE.df$pi[MLE.df$ll == max(MLE.df$ll)] #maxLL
            MLE.df$dev<- 2*(max(MLE.df$ll)-MLE.df$ll)
            lllow<- MLE.df$pi[which(MLE.df$dev < 3.841)[1]] #lower CI MAXLL
            llhi<- MLE.df$pi[which(MLE.df$dev < 3.841)[length(which(MLE.df$dev < 3.841))]] #higher CI MaxLL
            AR<- c(append(AR, maxll, after=length(AR)))
            ARlo<- c(append(ARlo, lllow, after=length(ARlo)))
            ARhi<- c(append(ARhi, llhi, after=length(ARhi)))        
          }
        }
        
        ARdf<- data.frame(Exposure_type=Exposure_type, AR=AR, ARlo=ARlo, ARhi=ARhi, ntotal=ntotal)
        ARdf_5total<- subset(ARdf, ntotal>= 5)
        
        #################################################
        #Calculate probability on infection for a contact of type x that has not developed symptoms up to time t
        
        data<-merge(data, ARdf, by="Exposure_type", all.x=T) 
        
        data$probinf<- (data$AR*data$probsymp)/(1-data$AR*data$prob)
        data$probinf<- ifelse(data$Type=="Case", 1, data$probinf)
        
        
        
        #################################################
        # Calculate reproduction number
        
        
        #first calculate the first component for all cases: add over the number of contacts of that case and wigh each contact by their probability of being infected
        
        rcase<- tapply(X=data$probinf, INDEX=data$IDSource, FUN=sum, simplify=T)
        rcase<- data.frame(ID=rownames(rcase), rcase=rcase)
        row.names(rcase) <- NULL 
        data<- merge(data, rcase, by="ID", all.x=T)
        
        #next calculate the second component for orphan cases: add over the number of orphans and weigh each by their probability of being infected with case j
        
        orphans<-subset(data, is.na(data$IDSource))
        
        cases<- subset(data, Type=="Case")
        for (i in 1: length(orphans$ID)){
          days<- as.numeric(as.Date(orphans$DOD[i])- as.Date(cases$DOD))
          orphan<- round(dgamma(days, generationtime, 1), 2)
          dftemp<- data.frame(ID= cases$ID, days=days, orphan=orphan)
          dftemp$proborphan<- round(dftemp$orphan/sum(dftemp$orphan),2)
          dftemp$proborphan<- ifelse(is.na(dftemp$proborphan), 0, dftemp$proborphan)
          dftemp$proborphan<- ifelse(is.nan(dftemp$proborphan), 0, dftemp$proborphan)
          colnames(dftemp)<- c("ID", paste("days",i, sep="_"), paste("orphan", i, sep="_"), paste("proborphan", i, sep="_"))
          data<- merge(data, dftemp, by="ID", all.x=T)
        }
        
        data$nsecondary<- round(data$rcase + data$proborphan_1 + data$proborphan_2 + data$proborphan_3, 2)
        data$nsecondary<- ifelse(is.na(data$nsecondary) & data$Type=="Case", 0, data$nsecondary)
        
        #apply moving average for plotting
        ma_dates<- data.frame(graphdate=data$graphdate[data$Type=="Case"], nsecondary= data$nsecondary[data$Type=="Case"]) 
        alldates<- data.frame(graphdate=c(seq(min(cases$graphdate),  max(data$date95, na.rm=T), by=1)))
        ma.dataframe<- merge(alldates, ma_dates, by="graphdate", all=T)
        alldates2<- data.frame(graphdate=c(seq(min(cases$graphdate),  currentdate, by=1)))
        
        mean<- c()
        sd<- c()
        nsample<- c()
        for (i in 1: length(alldates2$graphdate)){
          if (i<= generationtime) {
            temp.dataframe<- subset(ma.dataframe, subset= (graphdate >= alldates2$graphdate[1] & graphdate <= alldates2$graphdate[i]))
          } else {
            temp.dataframe<- subset(ma.dataframe, subset= (graphdate >= alldates2$graphdate[i-generationtime] & graphdate <= alldates2$graphdate[i-1])) 
          }        
          mean<- c(append(mean, mean(temp.dataframe$nsecondary, na.rm=T), after=length(mean)))
          sd<- c(append(sd, sd(temp.dataframe$nsecondary, na.rm=T), after=length(sd)))
          nsample<- c(append(nsample, sum(complete.cases(temp.dataframe)), after=length(nsample)))
        }
        
        ma.dataframe2<- data.frame(graphdate= alldates2$graphdate, ma=mean, sd=sd, n=nsample)
        ma.dataframe2<- merge(alldates, ma.dataframe2, by="graphdate", all.x=T)
        ma.dataframe2$lowerexact<- (qchisq(0.025, 2*ma.dataframe2$ma* ma.dataframe2$n)/2)/ma.dataframe2$n  #exact ma
        ma.dataframe2$upperexact<- (qchisq(0.975, 2*(ma.dataframe2$ma* ma.dataframe2$n+1))/2)/ma.dataframe2$n  #exact ma
        maxexact2<- max(ma.dataframe2$upperexact, na.rm=T) #for plotting, max y-axis
        ma.dataframe2$upperexact[ma.dataframe2$n <= 1]<- Inf
        ma.dataframe2$lowerexact[ma.dataframe2$n <= 1]<- -Inf
        ma.dataframe2$SEM<- ma.dataframe2$sd/sqrt(ma.dataframe2$n)
        ma.dataframe2$lowerSEM<- ma.dataframe2$ma - ma.dataframe2$SEM
        ma.dataframe2$higherSEM<- ma.dataframe2$ma + ma.dataframe2$SEM
        ma.dataframe2$lowerSEM[ma.dataframe2$n <= 1]<- Inf
        ma.dataframe2$higherSEM[ma.dataframe2$n <= 1]<- -Inf
        maxsEM<- max(ma.dataframe2$higherSEM, na.rm=T) #for plotting, max y-axis
        
        
        #Calculate Rt for specific time periods
        #calculate mean R0 before case detected
        casedetect<- subset(data, graphdate <= outbreakdetect & Type=="Case")
        #leave out first case
        casedetect<- arrange(casedetect, DOD)
        casedetect = casedetect[-1,]
        r.casedetect<- mean(casedetect$nsecondary, na.rm=T)
        ndetect<- length(casedetect$nsecondary)
        sddetect<- sd(casedetect$nsecondary, na.rm=T)
        lowerdetect<- (qchisq(0.025, 2*r.casedetect* ndetect)/2)/ndetect  #exact method poisson
        upperdetect<- (qchisq(0.975, 2*(r.casedetect* ndetect+1))/2)/ndetect  #exact method poisson
        lowerdetectSEM<- r.casedetect - (sddetect/sqrt(ndetect))
        upperdetectSEM<- r.casedetect + (sddetect/sqrt(ndetect))
        
        #calculate mean R0 after implementation control measures
        caseafter<- subset(data, graphdate >= controlstart & Type=="Case")
        #remove last case
        caseafter<- arrange(caseafter, desc(DOD))
        caseafter = caseafter[-1,]
        r.caseafter<- mean(caseafter$nsecondary, na.rm=T)
        ncaseafter<- length(caseafter$nsecondary)
        sdcaseafter<- sd(caseafter$nsecondary, na.rm=T)
        lowercaseafter<- (qchisq(0.025, 2*r.caseafter* ncaseafter)/2)/ncaseafter  #exact method poisson
        uppercaseafter<- (qchisq(0.975, 2*(r.caseafter* ncaseafter+1))/2)/ncaseafter  #exact method poisson
        lowercaseafterSEM<- r.caseafter - (sdcaseafter/sqrt(ncaseafter))
        uppercaseafterSEM<- r.caseafter + (sdcaseafter/sqrt(ncaseafter))
        
        #calculate no of secondary cases per case for graph
        count_df<- count(data$IDSource[data$Type=="Case"])
        data<- merge(data, count_df, by.x="ID", by.y="x", all.x=T)
        data$freq<- ifelse(is.na(data$freq)& data$Type=="Case", 0, data$freq)
        count_df2<- aggregate(data[,c("DOD","freq")], by= list(data$DOD, data$freq), FUN=length)
        
        
        
        #########################
        #data processing for visualisation
        
        
        #determine the roots for every cluster
        data$Root<- ifelse(is.na(data$IDSource), 1, 0)
        
        
        #calculate days since first case, to determine xlim of graph
        mindatum<- min(data$graphdate, na.rm=T)
        data$daysdatum<- data$graphdate - mindatum
        
        #For graph, only include the contacts, which still pose a risk, probsymp>0
        
        data4<- subset(data, Type=="Case"| (Type=="Contact" & probsymp >= problow))
        data4$IDSource2<- as.character(data4$IDSource)
        
        
        ####################################
        #use algorithm to layout the links between cases and contacts for network diagram
        
        #maak edges data frame, met edges en de attributes
        edges<- data.frame(IDSource= data4$IDSource2, ID=data4$ID)
        edges2<- na.omit(edges)
        
        #maak vertices data frame, met vertices en de attributes
        vertices<- data.frame(ID= data4$ID, type=data4$Type, exposuretype= data4$Exposure_type, 
                gender= data4$Gender, birthyear= data4$Year_birth, days=data4$daysdatum, 
                root= data4$Root, IDSource= data4$IDSource2, probsymp= data4$probsymp, 
                date05= data4$date05, date95= data4$date95, graphdate= data4$graphdate)
        
        
        infectiongraph <- graph.data.frame(edges2, directed=TRUE, vertices=vertices)
        
        
        #extract the IDs of the root-nodes
        roots<- data4$ID[data4$Root==1]
        roots<- roots[!is.na(roots)]
        roots<- as.character(roots)
        
        #determine which nodes (position) in the vertex dataframe are the roots
        positionroots<- match(roots, V(infectiongraph)$name)
        
        #make new layout for graph (tree layout)
        generationlayout<-layout.reingold.tilford(graph=infectiongraph, 
                root= positionroots, flip.y=T)
        
        #add dummy generations for layout (no crossing lines)
        #first add generation number
        generationdf<- data.frame(ID=V(infectiongraph)$name, generation_rev=generationlayout[,2], generation= 5-generationlayout[,2])
        vertices<- merge(vertices, generationdf, by="ID")
        #create dummy dataframe
        verticedummy<- data.frame(ID= vertices$ID, IDSource= vertices$IDSource, generation=vertices$generation, root=vertices$root)
        #create endpoint variable
        verticedummy$endpoint<- ifelse(verticedummy$ID %in% verticedummy$IDSource, 0, 1)
        
        for (i in 1:(5*sum(verticedummy$endpoint))){
                tryCatch({
                if (verticedummy$endpoint[i]==1 & (verticedummy$generation[i] < max(verticedummy$generation))){
                        verticedummy$endpoint[i]<- c(0)
                        root<- c(0)
                        newrow<- data.frame(ID=paste("dummy",i), IDSource= verticedummy$ID[i], generation= verticedummy$generation[i]+1, endpoint= verticedummy$endpoint[i]+1, root= root)
                        verticedummy<- rbind(verticedummy, newrow)       
                }
                }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                }
        edgesdummy<- data.frame(IDSource= verticedummy$IDSource, ID=verticedummy$ID)
        edgesdummy2<- na.omit(edgesdummy)
        infectiongraphdummy<- graph.data.frame(edgesdummy2, directed=T)
        #extract the IDs of the root-nodes
        rootsdummy<- verticedummy$ID[verticedummy$root==1]
        rootsdummy<- rootsdummy[!is.na(rootsdummy)]
        rootsdummy<- as.character(rootsdummy)
        
        #determine which nodes (position) in the vertex dataframe are the roots
        positionrootsdummy<- match(rootsdummy, V(infectiongraphdummy)$name)
        generationlayoutdummy<-layout.reingold.tilford(graph=infectiongraphdummy, 
                root= positionrootsdummy, flip.y=T)
        
        yposition<- data.frame(ID=V(infectiongraphdummy)$name, yID=generationlayoutdummy[,1])
        
        vertices.dataframe<- merge(vertices, yposition, by="ID", all.x=T)
        
        #involve time, by making a time-layout
        timelayout<- matrix(0, length(generationlayout[,1]), 2)
        timelayout[,1]<- V(infectiongraph)$days
        timelayout[,2]<- generationlayout[,1]
        
        
        #add time-layout parameters to the vertices dataframe
        vertices.dataframe<- cbind(vertices.dataframe, timelayout[,1])
        colnames(vertices.dataframe)[16]<-"xID"
        
        
        
        #create lookup table to add coordinates of source to every vertex
        lookup<- data.frame(IDSource=vertices.dataframe$ID, xsource= vertices.dataframe$xID, ysource= vertices.dataframe$yID)
        
        vertices.dataframe <- join(vertices.dataframe, lookup, by = "IDSource")
        vertices.dataframe$xsourcedate<- min(vertices.dataframe$graphdate) + vertices.dataframe$xsource
        vertices.dataframe$xIDdate<- min(vertices.dataframe$graphdate) + vertices.dataframe$xID
        
        
        
        
        #############################################
        
        #calculate indicators to add to plot
        
        ncase<- length(data$Type[data$Type=="Case"])
        ncontact<- length(data$Type[data$Type=="Contact"])
        ncontactzero<- length(data$Type[data$Type=="Contact" & data$probsymp < problow])
        ncontactfollow<- length(data$Type[data$Type=="Contact" & data$probsymp >= problow])
        nclusters<- sum(data$Root)
        transmission<- subset(vertices.dataframe, vertices.dataframe$exposuretype %in% ARdf_5total$Exposure_type)
        transmission$exposuretype<- factor(transmission$exposuretype, levels(droplevels(transmission$exposuretype)))
        typetransmissiontable<- table(transmission$exposuretype, transmission$type)
        labels<- c()
        for (i in (1:length(typetransmissiontable[,1]))){
                temp<- paste(rownames(typetransmissiontable)[i], " (", typetransmissiontable[i,1], ",", typetransmissiontable[i,2], ")")
                labels<- c(labels, temp)
        }
        gendertable<- table(vertices.dataframe$gender, vertices.dataframe$type)
        
       
        #plot
        
        vertices.dataframe$currentdate<- as.numeric(currentdate)
        vertices.dataframe$currentdate2<- currentdate
        vertices.dataframe$currentdatelabel<- "current date"
        vertices.dataframe$ID <- as.character(vertices.dataframe$ID)
        vertices.dataframe$IDSource <- as.character(vertices.dataframe$IDSource)
        vertices.dataframe$firstdiagnosis<- as.numeric(outbreakdetect)
        vertices.dataframe$firstdiagnosis2<- outbreakdetect
        vertices.dataframe$action<- as.numeric(controlstart)
        vertices.dataframe$action2<- controlstart
        
        
        dfrect<- data.frame(xmin=c(currentdate), xmax= c(max(vertices.dataframe$graphdate, na.rm=T)+3), ymin=c(min(vertices.dataframe$yID, na.rm=T)-3), ymax=c(max(vertices.dataframe$yID, na.rm=T)+5))
        dfrect2<- data.frame(xmin=c(currentdate), xmax= c(currentdate+1), ymin=c(min(vertices.dataframe$yID, na.rm=T)-3), ymax=c(max(vertices.dataframe$yID, na.rm=T)+5))
        textdf<- data.frame(label=c("past", "present", "future"), x= c(currentdate-2, currentdate + 0.5, currentdate+3), y=c(max(vertices.dataframe$yID, na.rm=T)+3.5, max(vertices.dataframe$yID, na.rm=T)+3.5, max(vertices.dataframe$yID, na.rm=T)+3.5), hjust= c(1, 0.5, 0 ))
        
        annotation<- data.frame(xmin=c(outbreakdetect, controlstart), xmax=c(outbreakdetect+1, controlstart+1), ymin=c(min(vertices.dataframe$yID, na.rm=T)-3), ymax=c(max(vertices.dataframe$yID, na.rm=T)+5))
        
        
        
        p<- ggplot(data=vertices.dataframe, aes(x=graphdate))+
          geom_segment(data=annotation, 
                       aes(x = xmin, xend = xmin,  y = ymin, yend = ymax),alpha=0.2, size=1, linetype="dashed")+
          geom_rect(data=dfrect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.1, inherit.aes=FALSE)+
          geom_rect(data=dfrect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.4, inherit.aes=FALSE)+
          geom_text(data=textdf, aes(x=x, y=y, label=label, hjust=hjust), color="black", size=3)+
          geom_segment(data=subset(textdf, label=="past"), aes(x = x, y = y-1, xend = x-3, yend = y-1), colour='black', size=0.4, arrow = arrow(length = unit(0.1, "cm")))+
          geom_segment(data=subset(textdf, label=="future"), aes(x = x, y = y-1, xend = x+3, yend = y-1), colour='black', size=0.4, arrow = arrow(length = unit(0.1, "cm")))+
          coord_cartesian(ylim = c(min(vertices.dataframe$yID, na.rm=T)-3, max(vertices.dataframe$yID, na.rm=T)+5),
                          xlim= c(min(vertices.dataframe$graphdate)-3, max(vertices.dataframe$date95, na.rm=T)+3), expand=F)+
          geom_segment(data= subset(vertices.dataframe, !is.na(IDSource) & IDSource!=ID), 
                       aes(x=xsourcedate, y=ysource, xend=xIDdate, yend=yID, col= exposuretype, linetype=type, drop=F), 
                       size=0.5)+
          scale_linetype_manual(values=c("solid", "dashed"), labels= c("case-case", "case-contact"), name="Type of contact")+
          scale_colour_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
                              name="Exposure type \n(# cases, # contacts in follow-up)",
                              labels= labels)+
          geom_segment(data=subset(vertices.dataframe, type == "Contact"), 
                       aes(x = date05, y = yID, xend = date95, yend = yID), 
                       colour='black', size=0.5)+
          geom_point(data=subset(vertices.dataframe, type == "Contact"), 
                     aes(x=date05, y=yID, shape=gender),
                     col="black", fill= "black", 
                     size=1)+
          geom_point(data=subset(vertices.dataframe, type == "Contact"), 
                     aes(x=date95, y=yID, shape=gender), 
                     col="black", fill= "black",
                     size=1)+
          scale_shape_manual(name= "Gender \n(# cases, # contacts in follow-up)", values=c(15, 23, 19), 
                             labels= c(paste("man (",gendertable[1,1],",", gendertable[1,2], ")"), paste("woman (",gendertable[2,1],",", gendertable[2,2], ")"), paste("unknown (",gendertable[3,1],",", gendertable[3,2], ")")),
                             drop=F)+
          geom_point(data=subset(vertices.dataframe, type=="Case"), 
                     aes(x=graphdate, y=yID, shape=gender), 
                     col="black", fill= "black", 
                     size=2)+
          geom_text(data=subset(vertices.dataframe, type == "Case"), 
                    aes(x=  xIDdate, y= yID+1.5, label= ID),
                    color="black", 
                    size=2)+
          geom_text(data=subset(vertices.dataframe, type == "Contact"), 
                    aes(x= date95+1, y= yID, label= ID),
                    color="black", 
                    size=2)+
                labs(y= "onzin titel", 
                        title= paste("Overview of the ", disease, " outbreak, ", country, "\n", min(vertices.dataframe$graphdate), " to ", currentdate))+
                theme_bw()+
                theme(axis.title.x=element_blank(), 
                        axis.title.y=element_text(size=10,face="bold", colour="white"),
                        axis.text.x=element_blank(),
                        axis.text.y=element_text(size=10, colour="white"),
                        axis.ticks.y= element_blank(),
                        axis.ticks.x= element_blank(),
                        legend.text= element_text(size=10),
                        legend.title= element_text(size=10, face="bold"),
                        legend.key = element_blank(),
                        legend.margin=unit(0, "lines"),
                        panel.grid.minor.y = element_blank(),
                        panel.grid.major.y = element_blank(),
                        panel.margin = unit(c(0,0,0,0), "line"),
                        plot.title= element_text(size=12, face="bold"),
                        plot.margin = unit(c(0,0,0,0.42), "line"))+
                guides(colour = guide_legend(order = 3), 
                        shape = guide_legend(order = 1),
                        linetype= guide_legend(order=2))
        
        
        p1 <- p + theme(legend.position="none") 
        
        p2<- ggplot(data=vertices.dataframe, aes(x=graphdate))+
                geom_segment(data=annotation, 
                        aes(x = xmin, xend = xmin,  y = ymin, yend = ymax),alpha=0.2, size=1, linetype="dashed")+
                geom_rect(data=dfrect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.1, inherit.aes=FALSE)+
                geom_rect(data=dfrect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.4, inherit.aes=FALSE)+
                coord_cartesian(ylim = c(-0.2, 10),
                        xlim= c(min(vertices.dataframe$graphdate)-3, max(vertices.dataframe$date95, na.rm=T)+3), expand=F)+
                geom_histogram(data=subset(vertices.dataframe, type=="Case"), 
                        aes(x=graphdate, ..count..),
                        binwidth= 3,
                        col="black", fill= "black", alpha=0.5 )+
                scale_y_continuous(breaks=c(seq(0,10, by= 2)))+
                labs(x="Date symptom onset",
                        y= "# cases")+
                theme_bw()+
                theme(axis.title.x=element_text(size=10,face="bold"), 
                        axis.title.y=element_text(size=10,face="bold"),
                        axis.text.x=element_text(size=10),
                        axis.text.y=element_text(size=10),
                        axis.ticks.y= element_line(colour="black"),
                        legend.text= element_text(size=10),
                        legend.title= element_text(size=10, face="bold"),
                        legend.key = element_blank(),
                        plot.title= element_text(size=12, face="bold"),
                        plot.margin = unit(c(0,0,0,0.57), "line"))+
                theme(legend.position="none")
        
        
        
        p3<- ggplot(data=ma.dataframe2, aes(x=graphdate))+
          geom_segment(data=annotation, 
                       aes(x = xmin, xend = xmin,  y = ymin, yend = ymax),alpha=0.2, size=1, linetype="dashed")+
          geom_rect(data=dfrect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.1, inherit.aes=FALSE)+
          geom_rect(data=dfrect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.4, inherit.aes=FALSE)+
          coord_cartesian(ylim = c(-0.5, maxexact2+0.5),
                          xlim= c(min(vertices.dataframe$graphdate)-3, max(vertices.dataframe$date95, na.rm=T)+3), expand=F)+
          geom_point(data=count_df2, aes(x=as.Date(Group.1), y=Group.2, size=freq), alpha=0.7, show.legend = FALSE)+
          scale_size_identity()+
          geom_ribbon(data=ma.dataframe2, aes(x=graphdate, ymax=upperexact, ymin=lowerexact), fill="red", alpha=0.2)+
          geom_line(data=ma.dataframe2, aes(x=graphdate, y=ma), col="black", size=0.6)+
          labs(x=" ", y="Rt and # secondary \ncases per case")+
          theme_bw()+
          theme(axis.title.x=element_text(size=7,face="bold"), 
                axis.title.y=element_text(size=7,face="bold"),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=7),
                axis.ticks.y= element_line(colour="black"),
                axis.ticks.x= element_blank(),
                legend.text= element_text(size=7),
                legend.key = element_blank(),
                legend.key.size= unit(c(0.4), "cm"),
                plot.title= element_text(size=10, face="bold"),
                plot.margin = unit(c(0,0,0,0), "line"))
        
        
        #plot attack rates
        ar<- ggplot(data=ARdf_5total, aes(x=Exposure_type, y=AR))+
          geom_bar(aes(x=Exposure_type, ymax = max(AR), fill=Exposure_type, drop=F), stat="identity",size=4, position=position_dodge(width=0.6))+
          geom_errorbar(aes(x=Exposure_type, ymin=ARlo, ymax=ARhi, drop=F), position=position_dodge(width=0.6), width=0.2, show.legend = FALSE)+
          scale_fill_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33"), guide=F)+
          labs(y= "Attack rate", x= "Exposure type") +
          theme_bw()+
          theme(axis.title.y= element_text(size=7,face="bold"),
                axis.text.y=element_text(size=7),
                axis.title.x= element_text(size=7,face="bold"),
                axis.text.x=element_text(size=7, angle=60, hjust=1),
                legend.text= element_text(size=7),
                legend.key = element_blank(),
                legend.key.size= unit(c(0.4), "cm"),
                plot.margin = unit(c(0,1,1,1), "line"))+
          theme(legend.position="none")  
        
        
        
        legend<- gtable_filter(ggplot_gtable(ggplot_build(p)), "guide-box") 
        
        text<- textGrob(x=unit(0.15, "npc"), y=unit(0.5, "npc"),label= paste("Summary statistics:",  "\n\n# cases: ", ncase, "\n# contacts:", ncontact, ", of which: \n",
                ncontactzero, " were no cases, ",ncontactfollow, " are in follow-up","\n# introductions: ", nclusters,
                "\nMax. # secondary cases per case:", round(maxexact2), 
                "\n\nImportant dates (dashed vertical lines):\n",
                outbreakdetect, ": ", input$date1label, "\n", controlstart, ": ", input$date2label
        ), hjust=0, gp= gpar(fontsize= 10))    
        
        
        grid.arrange(arrangeGrob(p1, p3, p2, nrow = 3, heights= c(7.5,1.5,1)), 
                arrangeGrob(text,legend,ar, nrow=3, heights=c(2.5,4,3.5)), 
                widths=unit.c(unit(0.95, "npc") - legend$width, unit(0.05, "npc") + legend$width), 
                nrow=1)
        }, width=1200, height=800)
        
    
        output$tabeldata <- DT::renderDataTable(DT::datatable({
                
                inFile <- input$file1
                
                if (is.null(inFile))
                        return(NULL)
                
                data <- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                        stringsAsFactors=F)
                
                data
        }, options= list(pageLength = 25)))

})
        

        