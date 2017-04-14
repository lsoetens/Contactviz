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
library(TTR)
library(reshape2)



shinyServer(function(input, output) {

        ranges <- reactiveValues(x = NULL, y = NULL)

        output$plotcontacttrace <- renderPlot({


        inFile <- input$file1

        if (is.null(inFile))
        return(NULL)


        data<- read.csv(inFile$datapath, header=input$header, sep=input$sep,
                        stringsAsFactors=F)

        withProgress(message = 'Processing...', value = 0, {

        # Increment the progress bar, and update the detail text.
        incProgress(0.1, detail = paste("Loading data..."))

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

        ############# Start of analysis #################


        #set all character-variables to lower case
        data$Type<- tolower(data$Type)
        data$Exposure_type<- tolower(data$Exposure_type)
        data$Gender<- tolower(data$Gender)
        names(data)<- tolower(names(data))


        #create some new variables and do some formatting
        data<- within (data, {
          id2<- as.character(id)
          idsource2<- as.character(idsource)
          date_onset_disease<- as.Date(dod, format="%Y-%m-%d")
          date_exposure<- as.Date(exposure_date, format="%Y-%m-%d")
          gender<- factor(gender, levels= c("man", "woman", "unknown"), ordered=TRUE)
          exposure_type<- factor(exposure_type, levels= sort(unique(exposure_type)))
        }
        )

        data$id<- data$id2
        data$idsource<- data$idsource2


        #Remove cases and their contacts if dod of case is after currentdate and exposure type is unknown
        idselect<- data$id[data$type == "case" & data$date_onset_disease> currentdate]
        data<- subset(data, !(data$id %in% idselect) & !(data$idsource %in% idselect))


        #set cases diagnosed after the currentdate as contacts and remove their contacts
        idselect<- data$id[data$type=="case" & data$date_onset_disease > currentdate]
        data<- subset(data, !(data$idsource %in% idselect))
        data$type<- ifelse(data$id %in% idselect, "contact", data$type)

        # Increment the progress bar, and update the detail text.
        incProgress(0.3, detail = paste("Estimating probabilities..."))
        #########################
        #calculate Feff for every contact

        #first select contacts
        data2<- subset(data, type== "contact")


        #create function for creating convolution of exposure interval and incubation period distribution
        #function returns distribution for every contact
        convolution<- function(expdate, expduration, nsim){
          set.seed(200)
          x<- expdate + runif(nsim, min=0, max= expduration)+ rlnorm(nsim, log(incubationperiod), log(incubationsd))
          return(x)
        }


        #returns a matrix, with in every column the simulated dates for one contact
        simdates<- mapply(convolution, expdate= data2[,"date_exposure"], expduration= data2[,"exposure_duration"], nsim=10000 )

        #apply ecdf function to every column in simdates
        cumsim<- apply(simdates, 2, ecdf)


        #calculate cum. probability (and 1- prob) for certain date for every column
        arglist<- as.list(c(rep(as.numeric(currentdate), length(data2[, "id"]))))
        prob<- mapply(function(x,y) x(y), x= cumsim, y=arglist) #feff
        probsymp<- round(1-prob,6) #Seff
        date05<- as.Date(round(as.numeric(apply(simdates, 2, quantile, prob=problow))),origin = "1970-01-01")
        date95<- as.Date(round(as.numeric(apply(simdates, 2, quantile, prob=probhigh))),origin = "1970-01-01")


        #attach the calculated probabilities to original dataframe
        data3<- data.frame(id=data2$id, prob= prob, probsymp=probsymp, date05= date05, date95= date95)
        data<- merge(data, data3, by="id", all=T)


        # make a graphdate, which corresponds to the 1st day of symptom onset for cases,
        # and to the lower 5% of the symptom onset interval for contacts
        data<- within (data, {
          graphdate<- ifelse(type=="case", date_onset_disease, date05)
          graphdate<- as.Date(graphdate, origin = "1970-01-01")
        })


        #########################################
        #MLE for probability on infection for every contact type

        # Increment the progress bar, and update the detail text.
        incProgress(0.5, detail = paste("Estimating attack rates..."))

        LE<- function(ns, na, feff, pi){
          x<- ns*log10(pi)+ na*log10(1-pi) + sum(log10(1-pi*feff))
          return(x)
        }



        ar<- c()
        arlo<- c()
        arhi<- c()
        exposure_type<- c()
        ntotal<- c()
        ncase<- c()
        ncontact_no_symp<- c()

        for(i in 1: length(unique(data$exposure_type))){
          test<- subset(data, exposure_type== unique(data$exposure_type)[i])
          exposure_type<- c(append(exposure_type, as.character(unique(test$exposure_type)), after=length(exposure_type)))
          ntotal<- c(append(ntotal, length(test$exposure_type), after=length(ntotal)))
          ns<- length(test$type[test$type== "case"]) #number of cases
          ncase<- c(append(ncase, ns, after= length(ncase)))
          na<- length(test$type[test$type=="contact" & test$probsymp< problow]) #number of contacts no symptoms in surveillance period
          ncontact_no_symp<- c(append(ncontact_no_symp, na, after= length(ncontact_no_symp)))
          pi<- c()
          ll<- c()

          if (length(test$exposure_type)> 5){
            for (j in 1:1000){
              PI<- j/1000
              LL<- LE(ns=ns, na=na, feff=test$prob[test$type=="contact" & test$probsymp>= problow], pi=PI)
              pi<- c(append(pi, PI, after=length(pi)))
              ll<- c(append(ll, LL, after=length(ll)))
            }
            mle.df<- data.frame(pi= pi, ll=ll)
            maxll<- mle.df$pi[which(mle.df$ll == max(mle.df$ll, na.rm=T))] #maxLL
            mle.df$dev<- 2*(max(mle.df$ll, na.rm=T)-mle.df$ll)
            lllow<- mle.df$pi[which(mle.df$dev < 3.841)[1]] #lower CI MAXLL
            llhi<- mle.df$pi[which(mle.df$dev < 3.841)[length(which(mle.df$dev < 3.841))]] #higher CI MaxLL
            ar<- c(append(ar, maxll, after=length(ar)))
            arlo<- c(append(arlo, lllow, after=length(arlo)))
            arhi<- c(append(arhi, llhi, after=length(arhi)))
          }
          else {

            ar<- c(append(ar, NA, after=length(ar)))
            arlo<- c(append(arlo, NA, after=length(arlo)))
            arhi<- c(append(arhi, NA, after=length(arhi)))
          }
        }

        ardf<- data.frame(exposure_type=exposure_type, ar=ar, arlo=arlo, arhi=arhi, ntotal=ntotal, ncase=ncase, ncontact_no_symp= ncontact_no_symp)
        ardf_5total<- subset(ardf, ntotal>= 5)


        #################################################
        #Calculate probability on infection for a contact of type x that has not developed symptoms up to time t

        data<-merge(data, ardf, by="exposure_type", all.x=T)

        data$probinf<- (data$ar*data$probsymp)/(1-data$ar*data$prob)
        data$probinf<- ifelse(data$type=="case", 1, data$probinf)



        #################################################
        # Calculate reproduction number

        # Increment the progress bar, and update the detail text.
        incProgress(0.7, detail = paste("Estimating reproduction number..."))


        #first calculate the first and second component for all cases: add over the number of contacts of that case and weigh each contact by their probability of being infected

        rcase<- tapply(X=data$probinf, INDEX=data$idsource, FUN=sum, simplify=T)
        rcase2<- tapply(X=data$probinf[data$probinf> 0], INDEX=data$idsource[data$probinf>0], FUN= list)
        list_rcrude<- tapply(X=data$probinf[data$probinf== 1], INDEX=data$idsource[data$probinf== 1], FUN= list)
        rcase<- data.frame(id=rownames(rcase), rcase=rcase)
        row.names(rcase) <- NULL
        data<- merge(data, rcase, by="id", all.x=T)

        #next calculate the third component for orphan cases: add over the number of orphans and weigh each by their probability of being infected with case j

        orphans<-subset(data, is.na(data$idsource))

        cases<- subset(data, type=="case")
        for (i in 1: length(orphans$id)){
          days<- as.numeric(as.Date(orphans$dod[i])- as.Date(cases$dod))
          orphan<- round(dgamma(days, generationtime, 1), 2)
          proborphan<- round(orphan/sum(orphan),2)
          dftemp<- data.frame(id= cases$id, proborphan=proborphan)
          dftemp$proborphan<- ifelse(is.na(dftemp$proborphan), 0, dftemp$proborphan)
          dftemp$proborphan<- ifelse(is.nan(dftemp$proborphan), 0, dftemp$proborphan)
          colnames(dftemp)<- c("id", paste("proborphan", i, sep="_"))
          data<- merge(data, dftemp, by="id", all.x=T)
        }

        mean<- c()
        nsample<- c()
        lower.ci<- c()
        upper.ci<- c()

        alldates<- c(seq(min(cases$graphdate),  currentdate, by=1))

        flip.function<- function(k){
          n= length(k)
          p= k
          rbinom(n, 1, p)
        }

        for (i in 1: length(alldates)){
          if (i<= generationtime) {
            temp.dataframe<- subset(data, subset= (graphdate >= alldates[1] & graphdate <= alldates[i]))
          } else {
            temp.dataframe<- subset(data, subset= (graphdate >= alldates[i-generationtime] & graphdate <= alldates[i-1]))
          }

          # estimate R(t)
          id_temp<- temp.dataframe$id
          list_rcase<- rcase2[names(rcase2) %in% id_temp]
          orphandf<- temp.dataframe[, c(which(substr(colnames(temp.dataframe), 1, 10)== "proborphan"),  which(colnames(temp.dataframe)=="id"))]
          orphandf<- melt(orphandf, id.vars=c("id"))

          # remove NA's and zero probabilities
          orphandf$value<- ifelse(orphandf$value == 0, NA, orphandf$value)
          orphandf<- na.omit(orphandf)
          if (nrow(orphandf) >= 1){
            list_orphan <- tapply(X=orphandf$value, INDEX=orphandf$id, FUN=list)
            # merge the two lists
            both <- list(list_rcase, list_orphan)
            n <- unique(unlist(lapply(both, names)))
            names(n) <- n
            list_rcase_orphan<- lapply(n, function(ni) unlist(lapply(both, `[[`, ni)))
          } else {
            list_rcase_orphan<- list_rcase
          }

          # number of secondary cases per case
          list_rcase_orphan_sum<- unlist(lapply(list_rcase_orphan, FUN= sum))

          # add zero secondary cases for those with no secondary cases
          id_nosec<- c(data$id[!(data$id %in% data$idsource) & data$type == "case"])
          id_nosec_temp<- temp.dataframe$id[temp.dataframe$id %in% id_nosec]
          id_nosec_temp<- id_nosec_temp[!(id_nosec_temp %in% names(list_rcase_orphan_sum))]
          list_rcase_orphan_sum<- c(list_rcase_orphan_sum, rep(0, length(id_nosec_temp)))

          # R(t)
          mean<- c(append(mean, mean(list_rcase_orphan_sum, na.rm = T), after = length(mean)))

          # estimating the 95% CI around R(t)
          mean_rcase<- c()
          for (i in 1:1000){
            x<- i
            # step 1: for every possible secondary case, draw if it is going to be a case or not
            sample_rcase<- lapply(list_rcase_orphan, FUN= flip.function)
            # sum total number of secondary cases per case
            sum_rcase<- lapply(sample_rcase, FUN= sum)
            sum_rcase_unlist<- unlist(sum_rcase)
            # add zero secondary cases for those with no secondary cases
            id_nosec_temp<- temp.dataframe$id[temp.dataframe$id %in% id_nosec]
            id_nosec_temp<- id_nosec_temp[!(id_nosec_temp %in% names(sum_rcase_unlist))]
            sum_rcase_unlist<- c(sum_rcase_unlist, rep(0, length(id_nosec_temp)))
            # estimate mean from sample
            mean_rcase<- c(append(mean_rcase, mean(sum_rcase_unlist), after = length(mean_rcase)))
          }

          # determine 95% CI from sample means
          lower.ci<- c(append(lower.ci, quantile(mean_rcase, probs = c(0.025), na.rm = T), after = length(lower.ci)))
          upper.ci<- c(append(upper.ci, quantile(mean_rcase, probs = c(0.975), na.rm = T), after = length(upper.ci)))

        }

        ma.dataframe<- data.frame(graphdate= alldates, rt= mean, lower.ci= lower.ci, upper.ci = upper.ci)
        alldates<- data.frame(graphdate= alldates)
        ma.dataframe2<- merge(alldates, ma.dataframe, by="graphdate", all.x=T)


        #calculate no of secondary cases per case for graph
        count_df<- count(data$idsource[data$type=="case"])
        data<- merge(data, count_df, by.x="id", by.y="x", all.x=T)
        data$freq<- ifelse(is.na(data$freq)& data$type=="case", 0, data$freq)
        count_df2<- aggregate(data[,c("dod","freq")], by= list(data$dod, data$freq), FUN=length)

        maxupper<- max(count_df$freq)



        #########################
        #data processing for visualisation

        # Increment the progress bar, and update the detail text.
        incProgress(0.9, detail = paste("Plotting..."))


        #determine the roots for every cluster
        data$root<- ifelse(is.na(data$idsource), 1, 0)

        #calculate days since first case, to determine xlim of graph
        mindatum<- min(data$graphdate, na.rm=T)
        data$daysdatum<- data$graphdate - mindatum


        #For graph, only include the contacts, which still pose a risk, probsymp>0.05
        data_graph<- subset(data, type=="case"| (type=="contact" & probsymp > problow))
        data_graph$idsource2<- as.character(data_graph$idsource)


        ####################################
        #use algorithm to layout the links between cases and contacts for network diagram
        #maak edges data frame, met edges en de attributes
        edges<- data.frame(idsource= data_graph$idsource2, id=data_graph$id)
        edges2<- na.omit(edges)

        #maak vertices data frame, met vertices en de attributes
        vertices<- data.frame(id= data_graph$id, type=data_graph$type, exposuretype= data_graph$exposure_type,
          gender= data_graph$gender, days=data_graph$daysdatum,
          root= data_graph$root, idsource= data_graph$idsource2, probsymp= data_graph$probsymp,
          date05= data_graph$date05, date95= data_graph$date95, graphdate= data_graph$graphdate)


        suppressWarnings(
          #first network-step
          infectiongraph <- graph.data.frame(edges2, directed=TRUE, vertices=vertices)
        )


        #for tree-layout determine which are the roots
        #extract the ids of the root-nodes
        roots<- data_graph$id[data_graph$root==1]
        roots<- roots[!is.na(roots)]
        roots<- as.character(roots)


        #determine which nodes (position) in the vertex dataframe are the roots
        positionroots<- match(roots, V(infectiongraph)$name)

        #make new layout for graph (tree layout)
        generationlayout<-layout.reingold.tilford(graph=infectiongraph,
          root= positionroots, flip.y=T)


        #add dummy generations for layout (no crossing lines)
        #first add generation number
        generationdf<- data.frame(id=V(infectiongraph)$name, generation_rev=generationlayout[,2], generation= 5-generationlayout[,2])
        vertices<- merge(vertices, generationdf, by="id")
        #create dummy dataframe
        verticedummy<- data.frame(id= vertices$id, idsource= vertices$idsource, generation=vertices$generation, root=vertices$root)
        #create endpoint variable
        verticedummy$endpoint<- ifelse(verticedummy$id %in% verticedummy$idsource, 0, 1)

        #endpointdf<- subset(verticedummy, endpoint==1)
        #endpointdf$nrnewrecords<- max(endpointdf$generation)-endpointdf$generation


        # for (i in 1: (length(verticedummy$id)+ sum(endpointdf$nrnewrecords))){
        #
        #         tryCatch({
        #         if (verticedummy$endpoint[i]==1 & (verticedummy$generation[i] < max(verticedummy$generation))){
        #                 verticedummy$endpoint[i]<- c(0)
        #                 root<- c(0)
        #                 newrow<- data.frame(id=paste("dummy",i), idsource= verticedummy$id[i], generation= verticedummy$generation[i]+1, endpoint= verticedummy$endpoint[i]+1, root= root)
        #                 verticedummy<- rbind(verticedummy, newrow)
        #         }
        #         }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        # }
        #
        for (i in 1:(5*sum(verticedummy$endpoint))){
          tryCatch({
          if (verticedummy$endpoint[i]==1 & (verticedummy$generation[i] < max(verticedummy$generation))){
            verticedummy$endpoint[i]<- c(0)
            root<- c(0)
            newrow<- data.frame(id=paste("dummy",i), idsource= verticedummy$id[i], generation= verticedummy$generation[i]+1, endpoint= verticedummy$endpoint[i]+1, root= root)
            verticedummy<- rbind(verticedummy, newrow)
          }
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }


        edgesdummy<- data.frame(idsource= verticedummy$idsource, id=verticedummy$id)
        edgesdummy2<- na.omit(edgesdummy)
        infectiongraphdummy<- graph.data.frame(edgesdummy2, directed=T)
        #extract the ids of the root-nodes
        rootsdummy<- verticedummy$id[verticedummy$root==1]
        rootsdummy<- rootsdummy[!is.na(rootsdummy)]
        rootsdummy<- as.character(rootsdummy)

        #determine which nodes (position) in the vertex dataframe are the roots
        positionrootsdummy<- match(rootsdummy, V(infectiongraphdummy)$name)
        generationlayoutdummy<-layout.reingold.tilford(graph=infectiongraphdummy,
          root= positionrootsdummy, flip.y=T)

        yposition<- data.frame(id=V(infectiongraphdummy)$name, yid=generationlayoutdummy[,1])

        vertices.dataframe<- merge(vertices, yposition, by="id", all.x=T)


        #involve time, by making a time-layout
        timelayout<- matrix(0, length(generationlayout[,1]), 2)
        timelayout[,1]<- V(infectiongraph)$days
        timelayout[,2]<- generationlayout[,1]


        #add time-layout parameters to the vertices dataframe
        vertices.dataframe<- cbind(vertices.dataframe, timelayout[,1])
        colnames(vertices.dataframe)[ncol(vertices.dataframe)]<-"xid"


        #create lookup table to add coordinates of source to every vertex
        lookup<- data.frame(idsource=vertices.dataframe$id, xsource= vertices.dataframe$xid, ysource= vertices.dataframe$yid)

        vertices.dataframe <- join(vertices.dataframe, lookup, by = "idsource")
        vertices.dataframe$xsourcedate<- min(vertices.dataframe$graphdate) + vertices.dataframe$xsource
        vertices.dataframe$xiddate<- min(vertices.dataframe$graphdate) + vertices.dataframe$xid



        #############################################

        #calculate indicators to add to plot

        ncase<- length(data$type[data$type=="case"])
        ncontact<- length(data$type[data$type=="contact"])
        ncontactzero<- length(data$type[data$type=="contact" & data$probsymp < problow])
        ncontactfollow<- length(data$type[data$type=="contact" & data$probsymp >= problow])
        nclusters<- sum(data$root)
        transmission<- subset(vertices.dataframe, vertices.dataframe$exposuretype %in% ardf_5total$exposure_type)
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
        vertices.dataframe$id <- as.character(vertices.dataframe$id)
        vertices.dataframe$idsource <- as.character(vertices.dataframe$idsource)
        vertices.dataframe$firstdiagnosis<- as.numeric(outbreakdetect)
        vertices.dataframe$firstdiagnosis2<- outbreakdetect
        vertices.dataframe$action<- as.numeric(controlstart)
        vertices.dataframe$action2<- controlstart


        dfrect<- data.frame(xmin=c(currentdate), xmax= c(max(vertices.dataframe$date95, na.rm=T)+3), ymin=c(-Inf), ymax=c(Inf))
        dfrect2<- data.frame(xmin=c(currentdate), xmax= c(currentdate+1), ymin=c(-Inf), ymax=c(Inf))
        textdf<- data.frame(label=c("past", "present", "future"), x= c(currentdate-2, currentdate + 0.5, currentdate+3), y=c(max(vertices.dataframe$yid, na.rm=T)+4, max(vertices.dataframe$yid, na.rm=T)+4, max(vertices.dataframe$yid, na.rm=T)+4), hjust= c(1, 0.5, 0 ))

        annotation<- data.frame(xmin=c(outbreakdetect, controlstart), xmax=c(outbreakdetect+1, controlstart+1), ymin=c(min(vertices.dataframe$yid, na.rm=T)-3), ymax=c(max(vertices.dataframe$yid, na.rm=T)+5))


        p<- ggplot(data=vertices.dataframe, aes(x=graphdate))+
          geom_segment(data=annotation,
                       aes(x = xmin, xend = xmin,  y = ymin, yend = ymax),alpha=0.2, size=1, linetype="dashed")+
          geom_rect(data=dfrect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.1, inherit.aes=FALSE)+
          geom_rect(data=dfrect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.3, inherit.aes=FALSE)+
          geom_text(data=textdf, aes(x=x, y=y, label=label, hjust=hjust), color="black", size=3)+
          geom_segment(data=subset(textdf, label=="past"), aes(x = x, y = y-1, xend = x-3, yend = y-1), colour='black', size=0.4, arrow = arrow(length = unit(0.1, "cm")))+
          geom_segment(data=subset(textdf, label=="future"), aes(x = x, y = y-1, xend = x+3, yend = y-1), colour='black', size=0.4, arrow = arrow(length = unit(0.1, "cm")))+
          coord_cartesian(ylim = c(min(vertices.dataframe$yid, na.rm=T)-3, max(vertices.dataframe$yid, na.rm=T)+5),
                          xlim= c(min(vertices.dataframe$graphdate)-3, max(vertices.dataframe$date95, na.rm=T)+3), expand=F)+
          geom_segment(data= subset(vertices.dataframe, !is.na(idsource) & idsource!=id),
                       aes(x=xsourcedate, y=ysource, xend=xiddate, yend=yid, col= exposuretype, linetype=type, drop=F))+
          scale_linetype_manual(values=c("solid", "dashed"), labels= c("case-case", "case-contact"), name="type of contact")+
          scale_colour_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
                              name="Exposure type \n(# cases, # contacts in follow-up)",
                              labels= labels)+
          geom_rect(data=subset(vertices.dataframe, type == "contact"),
                    aes(xmin=date05, xmax=date95, ymin=yid-0.2, ymax=yid+0.2, fill = type), inherit.aes=FALSE)+
          scale_fill_manual(values = c("#99000d"), labels= c("Monitoring period of \ncontacts in follow-up"), name= " ")+
          # geom_segment(data=subset(vertices.dataframe, type == "contact"),
          #         aes(x = date05, y = yid, xend = date95, yend = yid),
          #         colour='black', size=0.5)+
          # geom_point(data=subset(vertices.dataframe, type == "contact"),
          #         aes(x=date05, y=yid, shape=gender),
          #         col="black", fill= "black",
          #         size=1)+
          # geom_point(data=subset(vertices.dataframe, type == "contact"),
          #         aes(x=date95, y=yid, shape=gender),
          #         col="black", fill= "black",
          #         size=1)+
          geom_point(data=subset(vertices.dataframe, type=="case"),
                   aes(x=graphdate, y=yid, shape=gender),
                   col="black", fill= "black",
                   size=2)+
          scale_shape_manual(name= "gender \n(# cases, # contacts in follow-up)", values=c(15, 23, 19),
                    labels= c(paste(rownames(gendertable)," (", gendertable[,1],", ", gendertable[,2], ")", sep="")))+
          geom_text(data=subset(vertices.dataframe, type == "case"),
                    aes(x=  xiddate, y= yid+1.5, label= id),
                    color="black",
                    size=2)+
          geom_text(data=subset(vertices.dataframe, type == "contact"),
                    aes(x= date95+1, y= yid, label= id),
                    color="black",
                    size=2)+
          labs(y= "onzin titel",
                    title= paste("Overview of the ", disease, " outbreak, ", country, "\n", min(vertices.dataframe$graphdate), " to ", currentdate))+
          theme_bw()+
          theme(axis.title.x=element_blank(),
                axis.title.y=element_text(size=7,face="bold", colour="white"),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=7, colour="white"),
                axis.ticks.y= element_blank(),
                axis.ticks.x= element_blank(),
                legend.text= element_text(size=7),
                legend.title= element_text(size=7, face="bold"),
                legend.key = element_blank(),
                legend.margin=unit(0, "lines"),
                panel.grid.minor.y = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.margin = unit(c(0,0,0,0), "line"),
                plot.title= element_text(size=10, face="bold"),
                plot.margin = unit(c(0,0,0,0.2), "line"))+
          guides(colour = guide_legend(order = 3),
                 shape = guide_legend(order = 1),
                 linetype= guide_legend(order=2))

        p1 <- p + theme(legend.position="none")

        p2<- ggplot(data=vertices.dataframe, aes(x=graphdate))+
          geom_segment(data=annotation,
            aes(x = xmin, xend = xmin,  y = ymin, yend = ymax),alpha=0.2, size=1, linetype="dashed")+
          geom_rect(data=dfrect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.1, inherit.aes=FALSE)+
          geom_rect(data=dfrect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.4, inherit.aes=FALSE)+
          coord_cartesian(ylim = c(-0.2, 20),
            xlim= c(min(vertices.dataframe$graphdate)-3, max(vertices.dataframe$date95, na.rm=T)+3), expand=F)+
          geom_histogram(data=subset(vertices.dataframe, type=="case"),
            aes(x=graphdate, ..count..),
            binwidth= 3,
            col="black", fill= "black", alpha=0.5 )+
          #scale_y_continuous(breaks=c(seq(0,20, by= 2)))+
          scale_y_continuous(breaks= pretty_breaks())+
          labs(x="Date symptom onset",
            y= "# cases")+
          theme_bw()+
          theme(axis.title.x=element_text(size=7,face="bold"),
            axis.title.y=element_text(size=7,face="bold"),
            axis.text.x=element_text(size=7),
            axis.text.y=element_text(size=7),
            axis.ticks.y= element_line(colour="black"),
            legend.text= element_text(size=7),
            legend.title= element_text(size=7, face="bold"),
            legend.key = element_blank(),
            plot.title= element_text(size=10, face="bold"),
            plot.margin = unit(c(0,0,0,0.6), "line"))+
          theme(legend.position="none")


        p3<- ggplot(data=ma.dataframe2, aes(x=graphdate))+
          geom_segment(data=annotation,
                       aes(x = xmin, xend = xmin,  y = ymin, yend = ymax),alpha=0.2, size=1, linetype="dashed")+
          geom_rect(data=dfrect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.1, inherit.aes=FALSE)+
          geom_rect(data=dfrect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.4, inherit.aes=FALSE)+
          coord_cartesian(ylim = c(-0.5, maxupper+0.5),
                          xlim= c(min(vertices.dataframe$graphdate)-3, max(vertices.dataframe$date95, na.rm=T)+3), expand=F)+
          geom_point(data=count_df2, aes(x=as.Date(Group.1), y=Group.2, size=freq), alpha=0.7, show.legend = FALSE)+
          scale_size_identity()+
          geom_ribbon(data=ma.dataframe2, aes(x=graphdate, ymax=upper.ci, ymin=lower.ci), fill="red", alpha=0.2)+
          geom_line(data=ma.dataframe2, aes(x=graphdate, y=rt), col="black", size=0.6)+
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
                plot.margin = unit(c(0.5,0,0,0.35), "line"))


        #plot attack rates
        ar<- ggplot(data=ardf_5total, aes(x=exposure_type, y=ar))+
          geom_bar(aes(x=exposure_type, ymax = max(ar), fill=exposure_type, drop=F), stat="identity",size=4, position=position_dodge(width=0.6))+
          geom_errorbar(aes(x=exposure_type, ymin=arlo, ymax=arhi, drop=F), position=position_dodge(width=0.6), width=0.2, show.legend = FALSE)+
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

        text<- textGrob(x=unit(0.2, "npc"), y=unit(0.41, "npc"),label= paste("Summary statistics:",  "\n\n# cases: ", ncase, "\n# contacts:", ncontact, ", of which: \n",
                ncontactzero, " were no cases, ",ncontactfollow, " are in follow-up",
                "\nMax. # secondary cases per case:", round(maxupper),
                "\n\nImportant dates (dashed vertical lines):\n",
                outbreakdetect, ": ", input$date1label, "\n", controlstart, ": ", input$date2label
        ), hjust=0, gp= gpar(fontsize= 7))

        })


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


