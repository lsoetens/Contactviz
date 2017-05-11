library(shiny)

# Define UI for application that plots random distributions
shinyUI(fluidPage(
        tags$head(
          tags$style(HTML('#plotoutput{color: #fff; background-color: #337ab7; border-color: #2e6da4}')),
          tags$style(HTML('#downloadPlot{color: #fff; background-color: #337ab7; border-color: #2e6da4}'))
          ),

        titlePanel(h1("Contactviz: an analysis and visualization tool for contact tracing data")),
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
                                      radioButtons('sep', 'Separator',
                                                   c(Comma=',', Semicolon=';', Tab='\t'), ','),
                                      fileInput('file1', 'Choose CSV File',
                                                accept=c('text/csv',
                                                         'text/comma-separated-values,text/plain',
                                                         '.csv'))

                               )
                             )
                    ),
                    tabPanel(h4("Disease specific settings"),
                        fluidRow(
                                column(4,
                                        br(),
                                        selectInput("disease", label = "Infectious disease of study",
                                                   choices = list("Measles" = "measles", "Mumps" = "mumps", "RSV" = "RSV", "SARS" = "SARS", "Smallpox" = "smallpox", "Varicella" = "varicella", "Other..." = "other"),
                                                   selected = NULL),
                                        conditionalPanel(
                                         condition = "input.disease == 'other'",
                                         textInput("disease2", label = "Other namely:", value = "")
                                        )
                                ),
                                column(4,
                                        br(),
                                        textInput("generationinterval", label = "Mean generation interval in days", value = ""),
                                        h6(textOutput("generationref")),
                                        textInput("incubationperiod", label = "Mean incubation period in days", value = ""),
                                        h6(textOutput("incubationref")),
                                        textInput("sdincubationperiod", label = "Sd incubation period in days", value = ""),
                                        h6(textOutput("sdincubationref")),
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
                                        dateInput("currentdate", label = "Current date (yyyy-mm-dd)", value = "")
                                )#,
                                # column(4,
                                #         br(),
                                #         h5(strong("Important dates for annotation of the plot (optional):")),
                                #         dateInput("date1", label = "Important date 1 (yyyy-mm-dd)", value = ""),
                                #         textInput("date1label", label = "Label date 1", value = "", placeholder = "E.g. Outbreak detected"),
                                #         dateInput("date2", label = "Important date 2 (yyyy-mm-dd)", value = ""),
                                #         textInput("date2label", label = "Label date 2", value = "", placeholder = "E.g. Control measures")
                                #         )
                        )
                        )

                ),

        hr(),

        p(textOutput("validation")),

        actionButton("plotoutput", "Show output in plot (takes about 20 seconds)"),


        hr(),

        tabsetPanel(type=c("tabs"),
                tabPanel(h4("Data tabel"),
                             br(),
                             DT::dataTableOutput("tabeldata")
                        ),
                tabPanel(h4("Plot"),
                        br(),
                        plotOutput("plotcontacttrace", width=1200, height=800),
                        downloadButton("downloadPlot", "Download")

                        )

        )

        )
)
