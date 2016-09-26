library(shiny)

shinyUI(bootstrapPage(
  
  fluidPage( title= "SEQualyzer",
             theme = "bootstrap.css",
             
             div(class="logo_back", fluidRow(img(id="logo", src="logo.jpg", align = "left", width=150, height=50))),
             
             tagList(
               
               tags$head(
                 
                 
                 tags$style(HTML("
                                 
                                 .logo_back { 
                                 background-color: #90C3D4;
                                 #height: 40px;
                                 #border: 0px;
                                 position: relative;
                                 #top: 40px;
                                 #bottom: 70;
                                 }
                                 
                                 #logo {
                                 margin-left: 40px;
                                 }
                                 
                                 .body {
                                 margin-top: 20px;
                                 }
                                 
                                 #head_plot {
                                 border: 2px;
                                 border-color: black;
                                 }
                                 
                                 #the_plot {
                                 margin-top: 10px;
                                 margin-bottom: 10px;
                                 margin-left: 5%;
                                 border: 5px;
                                 border-color: black;
                                 
                                 }
                                 
                                 #correlation_replicates {
                                 margin-top: 15px;
                                 
                                 }
                                 
                                 #replicate_SNR {
                                 margin-top: 10px;
                                 }

      
                                 "))
               )
               
             ),
             
             div(id="head_plot" , 
#                  uiOutput(outputId="version_warningUI"),
                 uiOutput(outputId= "main_plot")),
             
             div(class="body",
                 sidebarLayout(
                   
                   sidebarPanel(
                     
                     textInput("session_name", label = h3("Session name"), 
                               value = gsub(":", ".", 
                                            paste(c("Untitled1", strsplit(date(), split=" ")[[1]]), 
                                                  collapse="_"))),
                     
                     #                      selectInput(inputId= "probe", 
                     #                                  label= "Probing agent",
                     #                                  choices= c("SHAPE", "DMS"),
                     #                                  selected= "SHAPE"),
                     checkboxInput('sequence_info', 'Sequence', FALSE),
                     
                     uiOutput(outputId= "probed_baseUI"),
                     
                     
                     #                      selectInput(inputId= "priming_strategy",
                     #                                  label= "Priming strategy",
                     #                                  choices= c("Single primer", "Random primer"),
                     #                                  selected= "Random primer"),
                     
                     #                      radioButtons(inputId= "read_type",
                     #                                   label= "Sequenced read type", 
                     #                                   choices= c("Single end", "Paired end"),
                     #                                   selected= "Paired end",
                     #                                   inline= TRUE),
                     
                     checkboxInput("local_coverage", "Local coverage", FALSE),
                     
                     numericInput(inputId= "replicate_count",
                                  label= "Enter number of replicates:",
                                  value= 3, min= 1, step=1),
                     
                     #                      checkboxInput("transcriptome_wide", "Multiple transcripts in dataset", FALSE),
                     
                     tags$hr(),
                     
                     fileInput(inputId= "file",
                               label= "Upload data",
                               multiple= FALSE),
                     
#                      radioButtons("file_rdata", "Import RData format", FALSE),
                     
                     radioButtons("file_rdata", label = "Import format:",
                                  choices = c("RMDB"= "RMDB", 
                                              "RDATA"= "RDATA",
                                              "Default"= "Default"), 
                                  selected = "Default", inline=T),
                     
                     actionButton(inputId= "submit",
                                  label= "Load"),
                     
                     #                      textOutput(outputId= "ff"),
                     
                     uiOutput(outputId= "filter_avg_covUI"),
                     
                     #                      fluidRow(
                     #                        column(5, textOutput(outputId= "f_avg_cov")),
                     #                        
                     #                        column(3, uiOutput(outputId= "filter_avg_covUI"))
                     #                        ),
                     
                     uiOutput(outputId= "filter_all_covUI"),
                     
                     uiOutput(outputId= "filter_lenUI"),
                     
                     # uiOutput(outputId= "filter_meanSNR_UI"), 
                     
                     checkboxInput(inputId= "display_filter_meanSNR", label="Enable filtering by mean SNR", value=FALSE),
                     # 
                     # uiOutput(outputId= "display_filter_meanSNR_min_UI"),
                     # 
                     # uiOutput(outputId= "display_filter_meanSNR_max_UI"),
                     
                     uiOutput(outputId= "sliders_mean_snr"),
                     
                     # uiOutput(outputId="choicesMsgUI"),
                     
                     uiOutput( outputId= "choose_transcript" ),
                     
                     
                     #                      uiOutput(outputId="")
                     
                     uiOutput(outputId="analyzeUI"),
                     
                     tags$hr(),
                     
                     checkboxInput("error_bars", "Show error bars", FALSE),
                     
                     
                     
#                      uiOutput(outputId= "error_methodUI"),
radioButtons(inputId= "error_method",
             label= "Method to calculate error for replicates", 
             choices= "Bootstrap",
             selected= "Bootstrap",
             inline= TRUE),
                     
                     uiOutput(outputId= "sliders"),
                     
                     uiOutput(outputId= "shift_by"),
                     
                     uiOutput(outputId= "shift_buttons"),
                     
                     uiOutput(outputId= "errorBtnUI"),
                     
                     #uiOutput(outputId= "formulaUI"),
                     
                     tags$hr(),
                     
                     actionButton(inputId="save_all", label= "Save analyses to hard disk"),
                     
                     uiOutput(outputId= "dwnPathBtnUI"),
                     
                     uiOutput(outputId= "winSizeUI"),
                     
                     uiOutput(outputId= "dwnldAllUI"),
                     
                     uiOutput(outputId= "createFileMsgUI"),
                     
                     
                     
                     tags$hr()
                   ),
                   
                   mainPanel(
                     
                     tabsetPanel(type = "tabs",
                                 
                                 tabPanel(title= "SEQualyzer summary",
                                          fluidRow(
                                            column(6, selectInput(inputId="estimate.flag",
                                                                  label= "Estimation method",
                                                                  choices= c("Difference of modification rates",
                                                                             "Ratio of modification rates"),
                                                                  selected= "Difference of modification rates")
                                            ),
                                            
                                            column(6, selectInput(inputId= "norm.flag",
                                                                  label= "Normalization strategy",
                                                                  choices= c("Do not normalize",
                                                                             "Boxplot normalization",
                                                                             "2%-8% normalization"),
                                                                  selected= "2%-8% normalization")
                                            )
                                          ),
                                          
                                          fluidRow(
                                            column(6, uiOutput(outputId= "log_transformUI")),
                                            
                                            column(6, uiOutput(outputId= "shift_zeroUI"))
                                          ),
                                          
                                          fluidRow(
                                            column(6, numericInput(inputId= "CQI_confidence",
                                                                   label= "Significance level for 95% CQI estimation",
                                                                   value= 0.95,
                                                                   min= 0.1,
                                                                   max= 0.99)
                                            ),
                                            
                                            column(6, numericInput(inputId= "CQI_percent_variation",
                                                                   label= "Percent fluctuation for 95% CQI estimation",
                                                                   value= 0.4, min= 0.1, max= 3)
                                            )
                                          ),
                                          
                                          uiOutput("CQI_UI"),
                                          
                                          radioButtons("fill_scheme", label = "Fill bars in barplots",
                                                       choices = c("by nucelotide"= "by nucleotide", 
                                                                   "by reactivity"= "by reactivity"), 
                                                       selected = "by reactivity", inline=T),
                                          plotOutput('plot_reactivities'),
                                          #tableOutput('dataset'),
                                          div(id= "replicate_SNR", plotOutput('replicateSNR')),
                                          
                                          div(id= "correlation_replicates",
{fluidRow(
  column(6, 
         textOutput('replicateCorr_title'),
         plotOutput('replicateCorr')),
  
  column(6, 
         #div(       
         
         textOutput('pairwiseSNR_title'),
         plotOutput('pairwiseSNR')
         #                                           tableOutput('pairwise_SNR'))
  )
)}),

downloadButton('downloadReactivityData', 'Save reactivities'),
#downloadButton('downloadBootstrapData', 'Save bootstrap resamples')
plotOutput(outputId= "error_stats"),
checkboxInput("window_with_ind", "Show window level SNR with individual replicate SNR", TRUE),
#tableOutput('mean_SNR'),
tableOutput('med_SNR')

#                                           plotOutput('replicateCorr'),
#                                           textOutput('pairwiseSNR_title'),
#                                           tableOutput('pairwise_SNR')
                                 ),

tabPanel(title= "Data visualization",
         plotOutput('reactivity_dist'),
         #uiOutput('ACGU_dist_UI'),
#          conditionalPanel(condition= "input.sequence_info",
                          plotOutput('ACGU_dist'),#),
         fluidRow(
#            column(6, numericInput(inputId="y_min_counts", 
#                                   label="Minimum for y-axis",
#                                   step=1, value=0)
#            ),
           
           column(6, textInput(inputId="y_max_counts_term_minus", label="Maximum for y-axis",
                               value = "")
#                   numericInput(inputId="y_max_counts", 
#                                   label="Maximum for y-axis",
#                                   step=1, value=NaN, min=10 )
           )
         ),
         plotOutput('terminations_minus'),
fluidRow(
  #            column(6, numericInput(inputId="y_min_counts", 
  #                                   label="Minimum for y-axis",
  #                                   step=1, value=0)
  #            ),
  
  column(6, textInput(inputId="y_max_counts_term_plus", label="Maximum for y-axis",
                      value = "")
         #                   numericInput(inputId="y_max_counts", 
         #                                   label="Maximum for y-axis",
         #                                   step=1, value=NaN, min=10 )
  )
),
         plotOutput('terminations_plus'),
#          conditionalPanel(condition= "input.local_coverage",
                          fluidRow(
                            #            column(6, numericInput(inputId="y_min_counts", 
                            #                                   label="Minimum for y-axis",
                            #                                   step=1, value=0)
                            #            ),
                            
                            column(6, textInput(inputId="y_max_counts_cov_minus", label="Maximum for y-axis",
                                                value = "")
                                   #                   numericInput(inputId="y_max_counts", 
                                   #                                   label="Maximum for y-axis",
                                   #                                   step=1, value=NaN, min=10 )
                            )
                          ),
                          plotOutput('coverage_minus'),
                          fluidRow(
                            #            column(6, numericInput(inputId="y_min_counts", 
                            #                                   label="Minimum for y-axis",
                            #                                   step=1, value=0)
                            #            ),
                            
                            column(6, textInput(inputId="y_max_counts_cov_plus", label="Maximum for y-axis",
                                                value = "")
                                   #                   numericInput(inputId="y_max_counts", 
                                   #                                   label="Maximum for y-axis",
                                   #                                   step=1, value=NaN, min=10 )
                            )
                          ),
                          plotOutput('coverage_plus')#)
),

tabPanel(title= "SNR distributions",
         #                                           div(id= "replicate_SNR", 
         plotOutput('replicate_SNR'),
         plotOutput('replicate_SNR_dist'),
         textOutput('pairwise_SNR_title'),
         plotOutput('pair_SNR'),
         plotOutput('reactivity_overlay'),
         numericInput(inputId="win_size", 
                      label="Enter the window size for window level SNR",
                      max=50, step=1, value=20),
         plotOutput(outputId= "resample_SNR"),
         #tableOutput('mean_SNR'),
         tableOutput('med_resample_SNR'),
         plotOutput('win_replicate_SNR'),
         plotOutput(outputId= "window_SNR"),
         plotOutput(outputId= "resample_SNR_dist")
),

tabPanel(title= "Transcriptome-wide summaries",
         textOutput(outputId="tw_summary"),
         HTML("<br>"),
         textOutput(outputId="selection_tw_summary"),
         HTML("<br>"),
         plotOutput(outputId= "lorenz_curve"),
         HTML("<br>"),
         plotOutput(outputId= "hist_avg_cov")
         
         
)
                     )
                   )

                 )
             )
  ),
singleton(includeScript("active.js"))


)
)
