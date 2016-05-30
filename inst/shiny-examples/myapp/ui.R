library(shiny)

# Define UI for application that plots random distributions 
shinyUI(fluidPage(
        
        titlePanel(h1("Contact tracing overview tool")),
        br(),
        br(),
        
        tabsetPanel(type=c("tabs"),
                tabPanel(h4("Data upload"),
                        fluidRow(
                                column(6,
                                        br(),
                                        p("Please upload your csv file containing contact tracing data of the outbreak you want to study.
                                        In the file every row should contain a record of a case or contact."),
                                        br(),
                                        p("The csv file should contain at least the following columns:"),
                                        p(strong("- ID"),": the unique identifier of the case/contract"),
                                        p(strong("- IDSource"),": the unique identifier of the source of the case/contact"),
                                        p(strong("- Gender"),": the gender of the case/contact (man, women or unknown)"),
                                        p(strong("- Type"),": if it is a case or contact (Case or Contact)"),
                                        p(strong("- DOD"),": date of onset of disease (format: yyyy-mm-dd)"),
                                        p(strong("- Exposure_date"),": date of contact between case and contact (format: yyyy-mm-dd)"),
                                        p(strong("- Exposure_duration"),": duration of contact between case and contact (in days)"),
                                        p(strong("- Exposure_type"),": type of contact (in categories, no more than 7)"),
                                        br(),
                                        p("The order of the columns as well as the presence of additional columns is not important.")      
                                ),
                                column(4,
                                        br(),
                                fileInput('file1', 'Choose CSV File',
                                        accept=c('text/csv', 
                                        'text/comma-separated-values,text/plain', 
                                        '.csv')),
                                tags$hr(),
                                checkboxInput('header', 'Header', TRUE),
                                radioButtons('sep', 'Separator',
                                        c(Comma=',', Semicolon=';', Tab='\t'), ';')
                                )
                        )
                        ),
                tabPanel(h4("Disease specific settings"), 
                        fluidRow(
                                column(4,
                                        br(),
                                        textInput("disease", label = "Infectious disease of study", value = "Smallpox")
                                ),
                                column(4,
                                        br(),
                                        textInput("generationinterval", label = "Generation interval in days", value = "18"),
                                        textInput("incubationperiod", label = "Mean incubation period in days", value = "13.0"),
                                        textInput("sdincubationperiod", label = "Sd incubation period in days", value = "1.13"),
                                        p("Note: the incubation period distribution is assumed to be a log-normal distribution")
                                ),
                                column(4,
                                        br(),
                                        sliderInput("surveillance", label = "Percentage of the probability on developing symptoms the surveillance period 
                                                should cover:", min=0, max=100, value= 90)
                                )
                        )
                        ), 
                tabPanel(h4("Outbreak specific settings"), 
                        fluidRow(
                                column(4,
                                        br(),
                                        textInput("country", label = "Country/region of outbreak", value = "The Netherlands"),
                                        dateInput("currentdate", label = "Current date (yyyy-mm-dd)", value = "1951-05-30")
                                ),
                                column(4, 
                                        br(),
                                        h5(strong("Important dates for annotation of the plot:")),
                                        dateInput("date1", label = "Important date 1 (yyyy-mm-dd)", value = "1951/04/24"),
                                        textInput("date1label", label = "Label date 1", value = "Outbreak detected"),
                                        dateInput("date2", label = "Important date 2 (yyyy-mm-dd)", value = "1951/04/29"),
                                        textInput("date2label", label = "Label date 2", value = "Control measures")
                                        )
                        )
                        )
                ),
        
        hr(),
        
        submitButton("Show output"),
        
        hr(),
        
        tabsetPanel(type=c("tabs"),
                tabPanel(h4("Plot"),
                        br(),
                        plotOutput("plotcontacttrace", width=1200, height=800)#,
                        # radioButtons("var3", label= "Select the file type for download", 
                        #         choices= list("png", "pdf")),
                        # downloadButton("down", label = "Download the plot"),
                        # br()
                ),
                tabPanel(h4("Data tabel"),
                        br(),
                        DT::dataTableOutput("tabeldata")
                        )
        )
        
        )
)
