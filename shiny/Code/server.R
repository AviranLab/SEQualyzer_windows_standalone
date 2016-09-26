if (getRversion() < "3.2.2") {
  message <- paste("You are running R version ", as.character(getRversion()), 
                   ". R version 3.2.2 or more recent is required. The software may not run properly with your version.\n")
  warning(message)
}

options(shiny.maxRequestSize = 100*1024^2)

list.of.packages <- c("ggplot2", "reshape2", "tools", "corrplot", "foreach", "parallel", "doParallel", "plyr", "seqinr", "rmngb", "ineq", "ttutils", "zoo")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
tryCatch({
  if(length(new.packages)) install.packages(new.packages, repos="http://cran.rstudio.com/")
}, error= function(e) {
  print(paste0("Cannot install required packages. Possible cause - No internet connection or CRAN repository inaccessible. Please check your internet connection or try installing following package manually: ", new.packages))
})

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) stop("Required packages missing. Execution aborted.")



source("bootstrap_analyses.R")
source("reactivity_calculator.R")
source("CQI.R")
source("file_reader.R")
source("meanSNR.R")


library(ggplot2)
library(reshape2)
library(tools)
library(corrplot)
library(foreach)
library(parallel)
library(doParallel)
library(seqinr)
library(rmngb)
library(plyr)
library(ineq)
library(ttutils)
library(zoo)
registerDoParallel(cores=detectCores(all.tests=TRUE))


clip_at <- 35
low_t <- 0.3
high_t <- 0.7

shinyServer(function(input, output, session) {
  
  version = as.character(getRversion())
  session$sendCustomMessage(type = "version_checkHandler", version)
  
  output$version_warning <- renderText( paste("Your version of R: ", 
                                              as.character(getRversion()), 
                                              ". R version 3.2.2 or higher is required. The software may not function properly.\n") )
  
  output$probed_baseUI <- renderUI({
    conditionalPanel(condition= "input.sequence_info",
                     checkboxGroupInput(inputId= "probed_base", label= "Bases probed",
                                        choices= c("A","C","G","U"), 
                                        selected = c("A","C","G","U"), inline = T))
  })
  
  observe({
    input$local_coverage
    if (input$local_coverage) updateRadioButtons(session, "error_method", choices = c("Bootstrap", "Use formula"))
  })
  
  output$log_transformUI <- renderUI({
    selectInput(inputId= "log.flag",
                label= "Log transformation",
                choices= if (input$error_method == "Use formula") "Do not transform" else c("Do not transform",
                                                                                            "Transform modification counts and passes",
                                                                                            "Transform reactivities"),
                selected= "Do not transform")
  })
  
  output$shift_zeroUI <- renderUI({
    selectInput(inputId= "zero.flag",
                label= "Handling negatives",
                choices= if (input$error_method == "Use formula") "Shift to zero" else c("Shift to zero",
                                                                                         "Do not shift"),
                selected= "Shift to zero")
  })
  
  est.flag <- reactive({
    if (input$estimate.flag == "Difference of modification rates") return(1)
    if (input$estimate.flag == "Ratio of modification rates") return(2)
  })
  
  nor.flag <- reactive({
    if (input$norm.flag == "Do not normalize") return(0)
    if (input$norm.flag == "Boxplot normalization") return(2)
    if (input$norm.flag == "2%-8% normalization") return(1)
  })
  
  logar.flag <- reactive({
    if (input$log.flag == "Do not transform") return(0)
    if (input$log.flag == "Transform modification counts and passes") return(1)
    if (input$log.flag == "Transform reactivities") return(2)
  })
  
  zeroing.flag <- reactive({
    if (input$zero.flag == "Shift to zero") return(1)
    if (input$zero.flag == "Do not shift") return(0)
  })
  
  oriData <- reactive(
    withProgress(message= "Reading files...", {
      
      if (input$file_rdata == "RDATA") {
        validate(need({!is.null(input[["file"]]) & if (!is.null(input[["file"]])) file_ext(input[["file"]][1, "name"]) == "RData" else FALSE}, "Please upload a valid .RData file. Look at the example datasets or manual for help."))
        load(input[["file"]][1, "datapath"], env = environment())
        validate(need({all(lapply(dat$dat, ncol) == 4*input$replicate_count + input$sequence_info) &identical(names(dat$dat), dat$dat_tw$name)},
                      "Invalid RData."))
        
      } else if (input$file_rdata == "RMDB") {
        readFile_name <- read_all_files(input[["file"]][1, "datapath"],
                                        input$replicate_count, input$sequence_info, 
                                        input$local_coverage, input$session_name, input$probed_base, "rdat")
        
        validate(need(readFile_name$error == "No errors." & (if (input$sequence_info) length(input$probed_base) else TRUE), message= "Please check the inputs."))
        load(readFile_name$file,
             env = environment())
        closeAllConnections()
      } else {
        readFile_name <- read_all_files(readLines(file(input[["file"]][1, "datapath"], open="r")),
                                        input$replicate_count, input$sequence_info, 
                                        input$local_coverage, input$session_name, input$probed_base, "StructureFold")
        
        validate(need(readFile_name$error == "No errors." & (if (input$sequence_info) length(input$probed_base) else TRUE), message= "Please check the inputs."))
        load(readFile_name$file,
             env = environment())
        closeAllConnections()
      }
      get("dat")
    }))
  
  observeEvent(input$submit, {
    
    output$filter_lenUI <- renderUI({
      validate(need(length(oriData()$dat_tw$name) > 1, message= paste0("1 transcript with length = ", 
                                                                       oriData()$dat_tw$length)))
      numericInput(inputId= "filter_len", label= "Length > ", value=0,
                   min = 0, step=1)
    })
    
    output$filter_avg_covUI <- renderUI({
      validate(need(length(oriData()$dat_tw$name) > 1, message= paste0("1 transcript with mean local coverage = ", 
                                                                       oriData()$dat_tw$avg_cov)))
      numericInput(inputId= "avg_cov", label="Mean local coverage > ", value=0,
                   min = 0, step=1)
    })
    
    output$filter_all_covUI <- renderUI({
      validate(need(length(oriData()$dat_tw$name) > 1, message= paste0("1 transcript with overall coverage = ", 
                                                                       oriData()$dat_tw$all_cov)))
      numericInput(inputId= "all_cov", label= "Overall coverage > ", value=0,
                   min = 0, step=1)
    })
    
    filters <- reactive({
      if (input$display_filter_meanSNR) {
        validate(need(!all(is.null(c(est.flag(), nor.flag(), logar.flag(), zeroing.flag(),
                                     input$replicate_count, input$local_coverage))),
                      message= "Wait."))
        list(len = input$filter_len,
             all_cov = input$all_cov,
             avg_cov = input$avg_cov,
             mean_SNR_min = input$range_view_mean_snr[1], #min_mean_SNR,
             mean_SNR_max = if (input$range_view_mean_snr[2] == round(max(mean_snr(), na.rm=T)-0.05, 2)) max(mean_snr(), na.rm=T) else input$range_view_mean_snr[2]) #max_mean_SNR)
      } else {
        list(len = input$filter_len,
             all_cov = input$all_cov,
             avg_cov = input$avg_cov) 
      }
      
    })
    
    mean_snr <- reactive(
      withProgress(message= "Scanning dataset for mean SNR...", {meanSNR(oriData(), input$session_name,
                                                                         est.flag(), 
                                                                         nor.flag(), logar.flag(), 
                                                                         zeroing.flag(),
                                                                         input$replicate_count, 
                                                                         input$local_coverage)}
      )
    )
    
    # create a return sliderInput, which changes based on the range of mean SNR in dataset. 
    output$sliders_mean_snr <- renderUI({
      validate(need(input$replicate_count>1, message="Only one replicate."))
      validate(need(!any(is.null(c(est.flag(), nor.flag(), logar.flag(), zeroing.flag(),
                                   input$replicate_count, input$local_coverage, mean_snr()))),
                    message= "Wait."),
               need(length(oriData()$dat_tw$name) > 1, message= paste0("1 transcript with mean SNR = ", 
                                                                       mean_snr())))
      sliderInput(inputId= "range_view_mean_snr", label="Range of mean SNR", 
                  min= round(min(mean_snr(), na.rm=T), 1), #0,
                  max=round(max(mean_snr(), na.rm=T)-0.05, 2), #clip_at, 
                  value=c(0, clip_at), step=0.05)
      
    })
    
    choiceList <- reactive({
      validate(need(!all(is.null(c(est.flag(), nor.flag(), logar.flag(), zeroing.flag(),
                                   input$replicate_count, input$local_coverage))),
                    message= "Wait."))
      if (input$display_filter_meanSNR & if (input$display_filter_meanSNR) filters()$mean_SNR_min != filters()$mean_SNR_max else TRUE) {
        list1 <- oriData()[["dat_tw"]]$name[which(oriData()$dat_tw$length >= filters()$len &
                                                    oriData()$dat_tw$all_cov >= filters()$all_cov &
                                                    oriData()$dat_tw$avg_cov >= filters()$avg_cov &
                                                    mean_snr() >= filters()$mean_SNR_min &
                                                    mean_snr() <= filters()$mean_SNR_max +0.05)] 
      } else {
        list1 <- oriData()[["dat_tw"]]$name[which(oriData()$dat_tw$length >= filters()$len &
                                                    oriData()$dat_tw$all_cov >= filters()$all_cov &
                                                    oriData()$dat_tw$avg_cov >= filters()$avg_cov)] 
      }
      
      if (!length(list1)) {
        
        updateNumericInput(session, inputId= "filter_len", label = "Length >", value = 0,
                           min = 0, step = 1)
        updateNumericInput(session, inputId= "avg_cov", label = "Mean local coverage > ", value = 0,
                           min = 0, step = 1)
        updateNumericInput(session, inputId= "all_cov", label = "Overall coverage > ", value = 0,
                           min = 0, step = 1)
        
        list1 <- oriData()$dat_tw$name
      }  
      
      list1
    })
    
    output$choose_transcript <- renderUI({
      selectInput(inputId="selected_transcript",
                  label= "Choose a transcript",
                  choices = choiceList())
    })
    
    output$selection_tw_summary <- renderText(paste(length(choiceList()), "selected."))
    
    output$analyzeUI <- renderUI({
      actionButton(inputId="btnAnalyze","Analyze")
    })
    
    #     rData
    data <- reactive({
      validate(need(if (!is.null(input$selected_transcript)) !is.null(oriData()$dat[[input$selected_transcript]]) else FALSE, message="Please check the input.")) 
      oriData()$dat[[input$selected_transcript]] #else oriData()
    })
    
    observeEvent(input$btnAnalyze, {
      
      observe({
        input$selected_transcript
        input$file
        input$estimate.flag
        input$log.flag
        input$zero.flag
        input$norm.flag
        if (exists("reactivity_plot")) reactivity_plot$error_bars <- NULL
        if (exists("getSNR$values")) getSNR$values <- NULL
      })
      
      output$main_plot <- renderUI({
        fluidRow(
          div(id="the_plot", plotOutput("complete_plot", height = "200", width="95%"))
        )
      })
      
      # create a return sliderInput, which changes based on the number of residues of selected dataset. 
      output$sliders <- renderUI({
        sliderInput(inputId= "range_view", label="Select the range of view", min=1, max=nrow(data()), value=c(1, 50), round = TRUE, step=1)
        
      })
      
      output$shift_by <- renderUI({
        numericInput(inputId= "shift_by",
                     label= "Choose distance to shift",
                     value= 50,
                     step=1,
                     min=1, max= nrow(data())
        )
      })
      
      output$shift_buttons <- renderUI({
        p(    actionButton(inputId= "shift_left",
                           label= "<-"),
              actionButton(inputId= "shift_right",
                           label= "->")
        )
        
      })
      
      output$errorBtnUI <- renderUI({
        actionButton(inputId= "errorBtn", label= "Evaluate replicates separately")
      })
      
      # update sliderInput when user choose to shift the plot
      observeEvent(input$shift_right, {
        
        if (!is.null(data())) {
          if(input$range_view[2] + input$shift_by <= nrow(data()) ) {
            updateSliderInput(session, "range_view", value=c(input$range_view[1]+input$shift_by, input$range_view[2]+input$shift_by))
          }
          else if (input$range_view[1] +input$shift_by <= nrow(data())) {
            updateSliderInput(session, "range_view", value=c(input$range_view[1]+input$shift_by, nrow(data())))
          }
          else if (input$shift_by >= nrow(data())) {
            updateSliderInput(session, "range_view", value=c(1, nrow(data())))
          }
        }
        
      })
      
      #shift to left
      observeEvent(input$shift_left, {
        
        if (!is.null(data())) {
          if(input$range_view[1]-input$shift_by > 0){
            updateSliderInput(session, "range_view", value=c(input$range_view[1]-input$shift_by, input$range_view[2]-input$shift_by))
          }
          else if (input$range_view[2] - input$shift_by > 1) {
            updateSliderInput(session, "range_view", value=c(1, input$range_view[2] - input$shift_by))
          } else if(input$shift_by >= nrow(data())) {
            updateSliderInput(session, "range_view", value=c(1, nrow(data())))
          }  
        }
      })
      
      output$CQI_UI <- renderUI({
        validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                      message="Too many zero reactivities. Normalization constant is zero."))
        conditionalPanel(condition= "input.local_coverage",
                         textOutput('CQI_title'),
                         tableOutput('coverage_index'))
      })
      
      # save the range of view
      a <- reactive(return(input$range_view[1]))
      b <- reactive(return(input$range_view[2]))
      
      #break points for plots showing complete range of data
      ticks_complete <- reactive(c(1,floor(nrow(data())/10)*c(1:9)))
      
      #break points for plots showing selected range of data
      ticks_selected <- reactive(c(a(),a() + (floor((b()- a()+1)/10)*c(1:9)) ))
      
      #tc
      output$complete_plot <- renderPlot({
        
        if (input$sequence_info) {
          counts<- data.frame(Nucleotide= 1:nrow(data()), Counts=data()$plus.counts.1, sequence= data()$sequence)
          ggplot(counts, aes(x= Nucleotide, y= Counts, fill= sequence))+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete()) + annotate("rect", xmin=a(), xmax=b(), ymin=-Inf, ymax=Inf, alpha=0.1, fill="red")+
            ggtitle("(+) channel modification counts for Replicate_1") + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position")
        } else {
          counts<- data.frame(Nucleotide= 1:nrow(data()), Counts=data()$plus.counts.1)
          ggplot(counts, aes(x= Nucleotide, y= Counts))+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete()) + annotate("rect", xmin=a(), xmax=b(), ymin=-Inf, ymax=Inf, alpha=0.1, fill="red")+
            ggtitle("(+) channel modification counts for Replicate_1") + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") 
        }
        
      })
      
      #Coverage Quality Index
      observe ({
        output$CQI_title <- renderText("Coverage Quality Indices")
        output$coverage_index <- renderTable({
          if(!is.null(a()) & !is.null(b()) & input$local_coverage) CQI(data()[a():b(),], input$replicate_count, input$CQI_confidence, input$CQI_percent_variation)  
        })
      })
      
      #Reactivity calculation
      reactivityCalc <- reactive({
        validate(need(!is.null(data()), message="Please check input."),
                 need(!is.null(input$replicate_count), message="Please check input."),
                 need(!is.null(input$local_coverage), message="Please check input."))
        reactivities <- matrix(, nrow(data()), 1+ input$replicate_count)
        colnames(reactivities) <- c("Nucleotide", paste("Replicate", c(1: input$replicate_count), sep="_"))
        reactivities[,1] <- c(1:nrow(data()))
        for (j in 1: input$replicate_count) {
          if (input$local_coverage) {
            reactivities[, j+1] <- local.coverage.react.calc(data()[, paste0("minus.counts.", j)],
                                                             data()[, paste0("minus.pass.", j)],
                                                             data()[, paste0("plus.counts.", j)],
                                                             data()[, paste0("plus.pass.", j)],
                                                             est.flag(), nor.flag(), logar.flag(), zeroing.flag())
          } else {
            
            reactivities[, j+1] <- react.calculator(data()[, paste0("minus.counts.", j)],
                                                    data()[, paste0("plus.counts.", j)],
                                                    est.flag(), nor.flag(), logar.flag(), zeroing.flag())
            
          }
        }
        
        return(data.frame(reactivities))
        
      })
      
      #------------------------------------------------------------------
      #Formula calculations
      #------------------------------------------------------------------
      
      #calculate noise for formula based variance calculation
      noiseCalc <- reactive({
        validate(need(!is.null(data()), message="Please check input."),
                 need(!is.null(input$replicate_count), message="Please check input."),
                 need(!is.null(input$local_coverage), message="Please check input."))
        noise <- matrix(, nrow(data()), 1+ input$replicate_count)
        colnames(noise) <- c("Nucleotide", paste("Replicate", c(1: input$replicate_count), sep="_"))
        noise[,1] <- c(1:nrow(data()))
        for (j in 1: input$replicate_count) {
          eval(parse(text=paste("noise[, ", j+1, "] <- data()$minus.counts.", j, "/(data()$minus.counts.", j , " + data()$minus.pass.", j, ")", sep="")))
        }
        return(data.frame(noise))
      })
      
      #calculate raw reactivities for formula based calculation
      rawReactivityCalc <- reactive({
        reactivities <- matrix(, nrow(data()), 1+ input$replicate_count)
        colnames(reactivities) <- c("Nucleotide", paste("Replicate", c(1: input$replicate_count), sep="_"))
        reactivities[,1] <- c(1:nrow(data()))
        for (j in 1: input$replicate_count) {
          eval(parse(text=paste("reactivities[,j+1] <- local.coverage.react.calc(data()$minus.counts.", j, 
                                ", data()$minus.pass.", j, ", data()$plus.counts.", j, ", data()$plus.pass.", j,
                                ", est.flag(), 0, 0, 1)", sep="")))
        }
        return(data.frame(reactivities))
      })
      
      #calculate normalization constants for formula based variance calculation
      normCalc <- reactive({
        norm <- matrix(, 1, input$replicate_count)
        colnames(norm) <- paste("Replicate", c(1: input$replicate_count), sep="_")
        for (j in 1: input$replicate_count) {
          eval(parse(text=paste("norm[1, ", j, "] <- get.normalizer(rawReactivityCalc()[,", j+1,"],", nor.flag(), ")", sep="")))
        }
        return(data.frame(norm))
      })
      
      #calculate sdev from formula
      sdCalc <- reactive({
        formula_sd_calc(rawReactivityCalc(), noiseCalc(), normCalc(), data(), est.flag())
      })
      
      form_SNR_calc <- reactive({
        vals <- (reactivityCalc()/sdCalc())[-1]
        SNR <- cbind(data.frame(1:nrow(data())), replace(vals, which(vals> clip_at, arr.ind=T), clip_at))
        colnames(SNR) <- c("Nucleotide", paste("Replicate", c(1: input$replicate_count), sep="_"))
        return(SNR)
      })
      
      #------------------------------------------------------------------
      
      #------------------------------------------------------------------
      #Plot reactivities with error bars
      #------------------------------------------------------------------
      reactivity_reshaped <- reactiveValues(mean=NA, sdev=data.frame(NA) )
      reactivity_plot <- reactiveValues(means= NULL, error_bars= NULL)
      
      observe({
        updateCheckboxInput(session, inputId="error_bars", value = input$error_bars)
      })
      
      observe({
        validate(need(!is.null(reactivityCalc()), message="Please check input."))
        validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                      message="Too many zero reactivities. Normalization constant is zero."))
        
        reactivity_reshaped$mean = melt(reactivityCalc(), id.vars="Nucleotide", variable.name= "Replicate")
        data_plot = if (input$error_bars & if (nrow(reactivity_reshaped$sdev) >1) nrow(reactivity_reshaped$mean) == nrow(reactivity_reshaped$sdev) else FALSE ) data.frame(reactivity_reshaped$mean, sdev=reactivity_reshaped$sdev) else reactivity_reshaped$mean
        
        #         print(input$fill_scheme)
        if (input$fill_scheme == "by reactivity") {
          #           reac_category <- c()
          # print(cbind(data_plot$value))
          data_plot$sequence <- apply(cbind(data_plot$value), 1, FUN= function(x) {return(if (x< low_t) "L" else if (x > high_t) "H" else "M")})
          #             rep(data()$sequence, input$replicate_count)
          reactivity_plot$means <- ggplot(data_plot, aes(x= Nucleotide, y= value, fill= sequence))+ geom_bar(stat="identity")  +facet_grid( Replicate ~ ., scale="free_y") +
            ggtitle("SHAPE reactivities") + theme(plot.title= element_text(face="bold"))+xlab("Nucleotide position")+ylab("Reactivity") + 
            scale_fill_manual(values = c("red","black", "yellow"), labels = c(paste0(">= ", high_t),
                                                                              paste0("< ", low_t), paste0("< ", high_t, " & ", ">= ", low_t)))
          
          output$plot_reactivities <- renderPlot({ 
            reactivity_plot$means + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()), ylim=c(0, min(4, max(data_plot$value, na.rm=T)))) + if (input$error_bars) reactivity_plot$error_bars
          })
        }
        
        if (input$fill_scheme != "by reactivity" & input$sequence_info) {
          data_plot$sequence <- rep(data()$sequence, input$replicate_count)
          reactivity_plot$means <- ggplot(data_plot, aes(x= Nucleotide, y= value, fill= sequence))+ geom_bar(stat="identity")  +facet_grid( Replicate ~ ., scale="free_y") +
            ggtitle("SHAPE reactivities") + theme(plot.title= element_text(face="bold"))+xlab("Nucleotide position")+ylab("Reactivity")
          
          output$plot_reactivities <- renderPlot({ 
            reactivity_plot$means + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()), ylim=c(0, min(4, max(data_plot$value, na.rm=T)))) + if (input$error_bars) reactivity_plot$error_bars
          })
        } else if (input$fill_scheme != "by reactivity") {
          reactivity_plot$means <- ggplot(data_plot, aes(x= Nucleotide, y= value))+ geom_bar(stat="identity")  +facet_grid( Replicate ~ ., scale="free_y") +
            ggtitle("SHAPE reactivities") + theme(plot.title= element_text(face="bold"))+xlab("Nucleotide position")+ylab("Reactivity")
          
          output$plot_reactivities <- renderPlot({ 
            reactivity_plot$means + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()), ylim=c(0, min(4, max(data_plot$value, na.rm=T)))) +
              if (input$error_bars) reactivity_plot$error_bars
          })
        }
      })
      
      #------------------------------------------------------------------
      
      #------------------------------------------------------------------
      #SEQualyzer summary tab outputs
      #------------------------------------------------------------------
      output$downloadReactivityData <- downloadHandler(
        filename = function() {paste(input$session_name, "_reactivities.tsv", sep="")},
        content = function(file) {
          write.table(reactivityCalc(), file, sep="\t", row.names= FALSE)
        }
      )
      
      output$replicateSNR <- renderPlot({
        validate(need(input$replicate_count>1, message="Only one replicate."))
        validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                      message="Too many zero reactivities. Normalization constant is zero."))
        reactivities <- reactivityCalc()[-1]
        vals <- apply(reactivities, 1, mean, na.rm=TRUE)/apply(reactivities, 1, sd, na.rm=TRUE)
        SNR <- data.frame(Nucleotide= 1:nrow(reactivities), value= replace(vals, which(vals > clip_at), clip_at))
        if (input$sequence_info) {
          SNR$sequence =data()$sequence
          ggplot(SNR, aes(x= Nucleotide, y= value, fill=sequence))+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+
            ggtitle("Signal-to-Noise Ratio (comparing replicates)") + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value")
        } else {
          ggplot(SNR, aes(x= Nucleotide, y= value))+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+
            ggtitle("Signal-to-Noise Ratio (comparing replicates)") + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value")
        }
      })
      
      output$replicateCorr_title <- renderText("Correlation matrix for replicates")
      output$replicateCorr <- renderPlot({
        validate(need(input$replicate_count>1, message="Only one replicate."))
        validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                      message="Too many zero reactivities. Normalization constant is zero."))
        reactivities <- reactivityCalc()[-1]
        corrplot.mixed(cor(reactivities, use = "pairwise.complete.obs"),lower="number", upper="ellipse")
      })
      
      output$pairwiseSNR_title <- renderText("Pairwise SNR for replicates (means)")
      output$pairwiseSNR <- renderPlot({
        validate(need(input$replicate_count>1, message="Only one replicate."))
        validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                      message="Too many zero reactivities. Normalization constant is zero."))
        reactivities <- reactivityCalc()[-1]
        pairwiseSNR <- matrix(, input$replicate_count, input$replicate_count)
        colnames(pairwiseSNR) <- paste("Replicate", 1:input$replicate_count, sep="_")
        rownames(pairwiseSNR) <- paste("Replicate", 1:input$replicate_count, sep="_")
        for (i in 1:input$replicate_count) {
          for (j in 1:input$replicate_count) {
            if(i==j) {
              pairwiseSNR[i,j] <- 1
            } else {
              vals <- apply(reactivities[, c(i,j)], 1, mean, na.rm=TRUE)/apply(reactivities[,c(i,j)], 1, sd, na.rm=TRUE)
              pairwiseSNR[i,j] <- mean(replace(vals, which(vals>clip_at), clip_at), na.rm=TRUE)
            }
          }
        }
        corrplot.mixed(pairwiseSNR,lower="number", upper="ellipse",is.corr=FALSE)
      })
      
      #------------------------------------------------------------------
      
      #------------------------------------------------------------------
      #Bootstrapping/Formula calculations
      #------------------------------------------------------------------
      
      boot <- reactive(
        withProgress(message= "Bootstrapping...", {bootstrap(100, input$replicate_count, data(), 
                                                             !input$local_coverage)}
        )
      )
      
      getSNR <- reactiveValues(values= NULL)
      
      #------------------------------------------------------------------
      
      observeEvent(input$errorBtn, {
        observe({
          input$selected_transcript
          input$file
          input$estimate.flag
          input$log.flag
          input$zero.flag
          input$norm.flag
          
          if (input$error_method == "Bootstrap") {
            booty <- boot()
            
            reactivities_boot <- getReactivities(booty, est.flag(), nor.flag(), logar.flag(), zeroing.flag(), !input$local_coverage)
            getSNR$values <- data.frame(SNRcalc(reactivities_boot))
            
          }  else {
            getSNR$values <- form_SNR_calc()
          }
          
          #------------------------------------------------------------------
          #SEQualyzer summary tab outputs
          #------------------------------------------------------------------
          
          updateCheckboxInput(session, inputId="error_bars", value = TRUE)
          
          output$med_SNR <- renderTable({
            validate(need(input$replicate_count>1, message=paste0("Only one replicate. Mean SNR for this transcript is ", round(mean(getSNR$values[, 2], na.rm=T), 2))))
            validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                          message="Too many zero reactivities. Normalization constant is zero."))
            central_SNR <- rbind(sapply(getSNR$values, median, na.rm=T),sapply(getSNR$values, mean, na.rm=T))[, -1]
            rownames(central_SNR) <- c("Median", "Mean")
            central_SNR
          })
          
          output$error_stats <- renderPlot({
            validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                          message="Too many zero reactivities. Normalization constant is zero."))
            SNR <- melt(getSNR$values, id.vars="Nucleotide", variable.name= "Replicate")
            window_SNR <- melt(SNR.scanner(getSNR$values, input$win_size), id.vars="Nucleotide", variable.name= "Replicate")
            validate(need(!all(is.na(SNR$value)), message= "Poor quality data."))
            if (input$sequence_info) {
              SNR$sequence <- rep(data()$sequence, input$replicate_count)
              window_SNR$sequence <- "A"
              if (input$window_with_ind) {
                ggplot(SNR, aes(x= Nucleotide, y= value, fill= sequence)) + geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+ facet_grid( Replicate ~ ., scale="free_y")+
                  ggtitle(paste("Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value") + geom_line(data=window_SNR, aes(x= Nucleotide, y= value)) + guides(colour=FALSE)
              } else {
                ggplot(SNR, aes(x= Nucleotide, y= value, fill= sequence)) + geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+ facet_grid( Replicate ~ ., scale="free_y")+
                  ggtitle(paste("Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value")
              }
              
            } else {
              if (input$window_with_ind) {
                ggplot(SNR, aes(x= Nucleotide, y= value)) + geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+ facet_grid( Replicate ~ ., scale="free_y")+
                  ggtitle(paste("Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value") + geom_line(data=melt(SNR.scanner(getSNR$values, input$win_size), id.vars="Nucleotide", variable.name= "Replicate"), aes(x= Nucleotide, y= value, colour="Red")) + guides(colour=FALSE)
              } else {
                ggplot(SNR, aes(x= Nucleotide, y= value)) + geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+ facet_grid( Replicate ~ ., scale="free_y")+
                  ggtitle(paste("Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value")
              }
            }
          })
          
          reactivity_reshaped$sdev <- if (input$error_method == "Bootstrap") sd.calc(reactivities_boot) else data.frame(sdev= melt(sdCalc(), id.vars="Nucleotide")$value)
          reactivity_plot$error_bars <- geom_errorbar(aes(ymin= value-sdev, ymax= value+sdev), colour= if (input$sequence_info) "black" else"red", width=0.5)
          
          #------------------------------------------------------------------
          
          #------------------------------------------------------------------
          #SNR distributions tab output (Need bootstrapping/formula)
          #------------------------------------------------------------------
          
          output$med_resample_SNR <- renderTable({
            validate(need(input$replicate_count>1, message="Only one replicate."))
            validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                          message="Too many zero reactivities. Normalization constant is zero."))
            central_SNR <- rbind(sapply(getSNR$values, median, na.rm=T),sapply(getSNR$values, mean, na.rm=T))[, -1]
            rownames(central_SNR) <- c("Median", "Mean")
            central_SNR
          })
          
          output$resample_SNR <- renderPlot({
            validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                          message="Too many zero reactivities. Normalization constant is zero."))
            SNR <- melt(getSNR$values, id.vars="Nucleotide", variable.name= "Replicate")
            if (input$sequence_info) {
              SNR$sequence <- rep(data()$sequence, input$replicate_count)
              ggplot(SNR, aes(x= Nucleotide, y= value, fill= sequence)) + geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+ facet_grid( Replicate ~ ., scale="free_y")+
                ggtitle(paste("Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value")
            } else {
              ggplot(SNR, aes(x= Nucleotide, y= value)) + geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+ facet_grid( Replicate ~ ., scale="free_y")+
                ggtitle(paste("Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value")
            }
          })
          
          output$window_SNR <- renderPlot({
            validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                          message="Too many zero reactivities. Normalization constant is zero."))
            window.SNR <- melt(SNR.scanner(getSNR$values, input$win_size), id.vars="Nucleotide", variable.name= "Replicate")
            ggplot(window.SNR, aes(x= Nucleotide, y= value))+ geom_line(aes(group=Replicate, colour= Replicate)) + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+scale_colour_brewer(palette="Set1")+
              ggtitle(paste("Window level Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value")
          })
          
          output$resample_SNR_dist <- renderPlot({
            validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                          message="Too many zero reactivities. Normalization constant is zero."))
            SNR <- melt(getSNR$values, id.vars="Nucleotide", variable.name= "Replicate")
            medians <-  ddply(SNR, .(Replicate), function(.d) data.frame(x=median(.d$value, na.rm=TRUE)))
            means <-  ddply(SNR, .(Replicate), function(.d) data.frame(x=mean(.d$value, na.rm=TRUE)))
            ggplot(data=SNR, aes(x=value)) + geom_histogram() + facet_grid( Replicate ~ ., scale="free_y") + ggtitle(paste(if (input$error_method== "Bootstrap") "Bootstrap" else "Formula", " SNR distributions", sep="")) + theme(plot.title= element_text(face="bold"))+ xlab("SNR value")+ 
              geom_vline(data=medians, aes(xintercept=x, colour="Median"), show_guide=TRUE) + geom_vline(data=means, aes(xintercept=x, colour="Mean"), show_guide=TRUE) +scale_colour_manual(name="Statistic", values=c("Median"= "red", "Mean"="green"))
          })
        })
        
        
        #------------------------------------------------------------------
      })
      
      #------------------------------------------------------------------
      #SNR distributions tab output (Don't need bootstrapping/formula)
      #------------------------------------------------------------------
      
      output$reactivity_overlay <- renderPlot({
        validate(need(input$replicate_count>1, message="Only one replicate."))
        validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                      message="Too many zero reactivities. Normalization constant is zero."))
        reactivities <- melt(reactivityCalc(), id.vars="Nucleotide", variable.name= "Replicate")
        ggplot(reactivities, aes(x= Nucleotide, y= value))+ geom_line(aes(group=Replicate, colour= Replicate)) + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+scale_colour_brewer(palette="Set1")+
          ggtitle("Reactivity scores from the replicates") + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("Reactivity score")
      })
      
      output$replicate_SNR <- renderPlot({
        validate(need(input$replicate_count>1, message="Only one replicate."))
        validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                      message="Too many zero reactivities. Normalization constant is zero."))
        reactivities <- reactivityCalc()[-1]
        SNR <- data.frame(Nucleotide= 1:nrow(reactivities), value=apply(reactivities, 1, function(x) min(mean(x)/sd(x), clip_at)))
        if (input$sequence_info) {
          SNR$sequence =data()$sequence
          ggplot(SNR, aes(x= Nucleotide, y= value, fill=sequence))+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+
            ggtitle("Signal-to-Noise Ratio (comparing replicates)") + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value")
        } else {
          ggplot(SNR, aes(x= Nucleotide, y= value))+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+
            ggtitle("Signal-to-Noise Ratio (comparing replicates)") + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value")
        }
      })
      
      output$win_replicate_SNR <- renderPlot({
        validate(need(input$replicate_count>1, message="Only one replicate."))
        validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                      message="Too many zero reactivities. Normalization constant is zero."))
        reactivities <- reactivityCalc()[-1]
        SNR <- data.frame(Nucleotide= 1:nrow(reactivities), value=apply(reactivities, 1, function(x) min(mean(x)/sd(x), clip_at)))
        SNR$value <- c(rep(NA, floor(input$win_size/2) ), 
                       rollapply(SNR$value, input$win_size, function(x) min(mean(x, na.rm=T), clip_at)),
                       rep(NA, floor((input$win_size-1)/2)))
        SNR$sequence =data()$sequence
        ggplot(SNR, aes(x= Nucleotide, y= value))+ geom_line() + scale_x_discrete(breaks=ticks_selected())+coord_cartesian(xlim = c(a(), b()))+
          ggtitle("Window level Signal-to-Noise Ratio (comparing replicates)") + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value")
      })
      
      output$replicate_SNR_dist <- renderPlot({
        validate(need(input$replicate_count>1, message="Only one replicate."))
        validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                      message="Too many zero reactivities. Normalization constant is zero."))
        reactivities <- reactivityCalc()[-1]
        SNR <- data.frame(Nucleotide= 1:nrow(reactivities), value=apply(reactivities, 1, mean, na.rm=TRUE)/apply(reactivities, 1, sd, na.rm=TRUE))
        ggplot(data=SNR, aes(x=value)) +geom_histogram()+ggtitle("Replicate SNR distribution") + theme(plot.title= element_text(face="bold"))+ xlab("SNR value")+
          geom_vline(aes(xintercept=median(value, na.rm=TRUE), colour="Median"), show_guide=TRUE)+geom_vline(aes(xintercept=mean(value, na.rm=TRUE), colour="Mean"), show_guide=TRUE) + scale_colour_manual(name="Statistic", values=c("Median"= "red", "Mean"="green")) 
      })
      
      output$pairwise_SNR_title <- renderText("Pairwise SNR for replicates (means)")
      output$pair_SNR <- renderPlot({ #renderTable({
        validate(need(input$replicate_count>1, message="Only one replicate."))
        validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                      message="Too many zero reactivities. Normalization constant is zero."))
        reactivities <- reactivityCalc()[-1]
        pairwiseSNR <- matrix(, input$replicate_count, input$replicate_count)
        colnames(pairwiseSNR) <- paste("Replicate", 1:input$replicate_count, sep="_")
        rownames(pairwiseSNR) <- paste("Replicate", 1:input$replicate_count, sep="_")
        for (i in 1:input$replicate_count) {
          for (j in 1:input$replicate_count) {
            if(i==j) {
              pairwiseSNR[i,j] <- 1
            } else {
              vals <- apply(reactivities[, c(i,j)], 1, mean, na.rm=TRUE)/apply(reactivities[,c(i,j)], 1, sd, na.rm=TRUE)
              pairwiseSNR[i,j] <- mean(replace(vals, which(vals>clip_at), clip_at), na.rm=TRUE)
            }
          }
        }
        corrplot.mixed(pairwiseSNR,lower="number", upper="ellipse",is.corr=FALSE)
      })
      #------------------------------------------------------------------
      
      #------------------------------------------------------------------
      #Data visualizations tab outputs
      #------------------------------------------------------------------
      
      observe({
        if (input$local_coverage) {
          output$coverage_minus <- renderPlot({
            validate(need(input$local_coverage, message="Local coverage information not available."))
            cov_minus <- data.frame(data()[grep(pat="minus.pass", colnames(data()))])
            colnames(cov_minus) <- paste("Replicate", 1:input$replicate_count, sep="_")
            cov_minus$Nucleotide <- c(1:nrow(data()))
            cov_minus <- melt(cov_minus, id.vars="Nucleotide", variable.name= "Replicate")
            
            if (input$sequence_info) {
              cov_minus$sequence <- rep(data()$sequence, input$replicate_count)
              ggplot(cov_minus, aes(x= Nucleotide, y= value, fill= sequence))+ geom_bar(stat="identity") +
                scale_x_discrete(breaks=ticks_selected())+
                coord_cartesian(xlim = c(a(), b()), 
                                ylim=c(0, 
                                       if (is.na(as.numeric(input$y_max_counts_cov_minus))) max(cov_minus$value, na.rm=T) else as.numeric(input$y_max_counts_cov_minus)))+
                facet_grid( Replicate ~ ., scale="free_y")+
                ggtitle("(-) channel local coverage")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Local coverage")
            } else {
              ggplot(cov_minus, aes(x= Nucleotide, y= value))+ geom_line(stat="identity") + 
                scale_x_discrete(breaks=ticks_selected())+
                coord_cartesian(xlim = c(a(), b()),
                                ylim=c(0, 
                                       if (is.na(as.numeric(input$y_max_counts_cov_minus))) max(cov_minus$value, na.rm=T) else as.numeric(input$y_max_counts_cov_minus)))+ 
                facet_grid( Replicate ~ ., scale="free_y")+
                ggtitle("(-) channel local coverage")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Local coverage")
            }
          })
          
          output$coverage_plus <- renderPlot({
            validate(need(input$local_coverage, message="Local coverage information not available."))
            cov_plus <- data.frame(data()[grep(pat="plus.pass", colnames(data()))])
            colnames(cov_plus) <- paste("Replicate", 1:input$replicate_count, sep="_")
            cov_plus$Nucleotide <- c(1:nrow(data()))
            cov_plus <- melt(cov_plus, id.vars="Nucleotide", variable.name= "Replicate")
            
            if (input$sequence_info) {
              cov_plus$sequence <- rep(data()$sequence, input$replicate_count)
              ggplot(cov_plus, aes(x= Nucleotide, y= value, fill= sequence))+ geom_bar(stat="identity") + 
                scale_x_discrete(breaks=ticks_selected())+
                coord_cartesian(xlim = c(a(), b()), 
                                ylim=c(0, 
                                       if (is.na(as.numeric(input$y_max_counts_cov_plus))) max(cov_plus$value, na.rm=T) else as.numeric(input$y_max_counts_cov_plus)))+ 
                facet_grid( Replicate ~ ., scale="free_y")+
                ggtitle("(+) channel local coverage")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Local coverage")
            } else {
              ggplot(cov_plus, aes(x= Nucleotide, y= value))+ geom_line(stat="identity") + 
                scale_x_discrete(breaks=ticks_selected())+
                coord_cartesian(xlim = c(a(), b()), 
                                ylim=c(0, 
                                       if (is.na(as.numeric(input$y_max_counts_cov_plus))) max(cov_plus$value, na.rm=T) else as.numeric(input$y_max_counts_cov_plus)))+ 
                facet_grid( Replicate ~ ., scale="free_y")+
                ggtitle("(+) channel local coverage")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Local coverage")
            }
          })
        }
      })
      
      
      output$terminations_minus <- renderPlot({
        term_minus <- data.frame(data()[grep(pat="minus.counts", colnames(data()))])
        colnames(term_minus) <- paste("Replicate", 1:input$replicate_count, sep="_")
        term_minus$Nucleotide <- c(1:nrow(data()))
        term_minus <- melt(term_minus, id.vars="Nucleotide", variable.name= "Replicate")
        
        if (input$sequence_info) {
          term_minus$sequence <- rep(data()$sequence, input$replicate_count)
          ggplot(term_minus, aes(x= Nucleotide, y= value, fill= sequence))+ geom_bar(stat="identity") + 
            scale_x_discrete(breaks=ticks_selected())+
            coord_cartesian(xlim = c(a(), b()), 
                            ylim=c(0, 
                                   if (is.na(as.numeric(input$y_max_counts_term_minus))) max(term_minus$value, na.rm=T) else as.numeric(input$y_max_counts_term_minus)))+ 
            facet_grid( Replicate ~ ., scale="free_y")+
            ggtitle("(-) channel modification counts")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Number of modifications")
        } else {
          ggplot(term_minus, aes(x= Nucleotide, y= value))+ geom_bar(stat="identity") + 
            scale_x_discrete(breaks=ticks_selected())+
            coord_cartesian(xlim = c(a(), b()), 
                            ylim=c(0, 
                                   if (is.na(as.numeric(input$y_max_counts_term_minus))) max(term_minus$value, na.rm=T) else as.numeric(input$y_max_counts_term_minus)))+ 
            facet_grid( Replicate ~ ., scale="free_y")+
            ggtitle("(-) channel modification counts")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Number of modifications")
        }
        
      })
      
      output$terminations_plus <- renderPlot({
        term_plus <- data.frame(data()[grep(pat="plus.counts", colnames(data()))])
        colnames(term_plus) <- paste("Replicate", 1:input$replicate_count, sep="_")
        term_plus$Nucleotide <- c(1:nrow(data()))
        term_plus <- melt(term_plus, id.vars="Nucleotide", variable.name= "Replicate")
        if (input$sequence_info) {
          term_plus$sequence <- rep(data()$sequence, input$replicate_count)
          ggplot(term_plus, aes(x= Nucleotide, y= value, fill= sequence))+ geom_bar(stat="identity") + 
            scale_x_discrete(breaks=ticks_selected())+
            coord_cartesian(xlim = c(a(), b()),  
                            ylim=c(0, 
                                   if (is.na(as.numeric(input$y_max_counts_term_plus))) max(term_plus$value, na.rm=T) else as.numeric(input$y_max_counts_term_plus)))+ 
            facet_grid( Replicate ~ ., scale="free_y")+
            ggtitle("(+) channel modification counts")+ 
            theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") +
            ylab("Number of modifications")
        } else {
          ggplot(term_plus, aes(x= Nucleotide, y= value))+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_selected())+
            coord_cartesian(xlim = c(a(), b()),  
                            ylim=c(0, 
                                   if (is.na(as.numeric(input$y_max_counts_term_plus))) max(term_plus$value, na.rm=T) else as.numeric(input$y_max_counts_term_plus)))+
            facet_grid( Replicate ~ ., scale="free_y")+
            ggtitle("(+) channel modification counts")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Number of modifications")
        }
      })
      
      output$reactivity_dist <- renderPlot({
        validate(need(if (input$replicate_count == 1) get.normalizer(rawReactivityCalc()[, -1], nor.flag()) != 0 else !any(apply(rawReactivityCalc()[, -1], 2, get.normalizer, flag=nor.flag()) == 0), 
                      message="Too many zero reactivities. Normalization constant is zero."))
        reactivities <- melt(reactivityCalc(), id.vars="Nucleotide", variable.name= "Replicate")
        medians <-  ddply(reactivities, .(Replicate), function(.d) data.frame(x=median(.d$value, na.rm=TRUE)))
        means <-  ddply(reactivities, .(Replicate), function(.d) data.frame(x=mean(.d$value, na.rm=TRUE)))
        ggplot(data=reactivities, aes(x=value)) + geom_histogram() + facet_grid( Replicate ~ ., scale="free_y") + ggtitle("Reactivity score distributions") + theme(plot.title= element_text(face="bold"))+ xlab("Reactivity score")+ 
          geom_vline(data=medians, aes(xintercept=x, colour="Median"), show_guide=TRUE) + geom_vline(data=means, aes(xintercept=x, colour="Mean"), show_guide=TRUE) +scale_colour_manual(name="Statistic", values=c("Median"= "red", "Mean"="green"))
      })
      
      output$ACGU_dist <- renderPlot({
        validate(need(input$sequence_info, message= "Sequence information not available."))
        term_minus <- data.frame(data()[c(grep(pat="minus.counts", colnames(data())), grep(pat="sequence", colnames(data())))], channel= "minus", Nucleotide= 1:nrow(data()) )
        names(term_minus)[grep(pat="minus.counts", colnames(term_minus))] <- paste("Replicate", file_ext(colnames(term_minus)[grep(pat="minus.counts", colnames(term_minus))]), sep="_")
        term_plus <- data.frame(data()[c(grep(pat="plus.counts", colnames(data())), grep(pat="sequence", colnames(data())))], channel= "plus", Nucleotide= 1:nrow(data()))
        names(term_plus)[grep(pat="plus.counts", colnames(term_plus))] <- paste("Replicate", file_ext(colnames(term_plus)[grep(pat="plus.counts", colnames(term_plus))]), sep="_")
        term <- melt(rbind(term_minus, term_plus), id.vars=c("Nucleotide", "sequence", "channel"), variable.name="Replicate")
        sum_term_1 <- data.frame(ddply(term, .(sequence, channel), function(.d) data.frame(Counts=sum(.d$value))), Replicate="All replicates")
        sum_term_2 <- ddply(term, .(sequence, Replicate, channel), function(.d) data.frame(Counts=sum(.d$value)))
        sum_term <- rbind(sum_term_1, sum_term_2)
        ggplot(data=sum_term, aes(x= sequence, y= Counts)) +geom_bar(stat="identity") +facet_grid(Replicate ~ channel, scale="free_y")            
      })
      
      #------------------------------------------------------------------
      
      #------------------------------------------------------------------
      #Transcriptome-wide tab outputs
      #------------------------------------------------------------------
      
      output$tw_summary <- renderText(paste(length(oriData()$dat_tw$avg_cov), "transcripts and", 
                                            round(sum(oriData()$dat_tw$all_cov, na.rm=T)), "reads."))
      
      
      
      output$lorenz_curve <- renderPlot({
        validate(need(length(oriData()$dat) > 1, message="Only one transcript in dataset."))
        distr <- cbind((0:length(oriData()$dat_tw$avg_cov))*(1/length(oriData()$dat_tw$avg_cov)), 
                       as.data.frame(lapply(oriData()$dat_tw[grep("^all_cov.", names(oriData()$dat_tw))],
                                            function(x) Lc(x, n = rep(1,length(x)), plot =F)[2])))
        colnames(distr) <- c("Fraction", paste(paste("Replicate", 1:input$replicate_count), "(Gini index =",
                                               lapply(oriData()$dat_tw[grep("^all_cov.", names(oriData()$dat_tw))],
                                                      function(x) format(round(ineq(x, type="Gini"), 2), nsmall = 2) ), ")"))
        distr <- melt(distr, id.vars= "Fraction", variable.name= "Replicate")
        
        ggplot(distr, aes(x= Fraction, y= value, color= Replicate)) +geom_line()+geom_abline()+
          xlab("Fraction of transcripts") +ylab("Cumulative coverage") + ggtitle("Lorenz curve of read distribution")
        
      })
      
      output$hist_avg_cov <- renderPlot({
        validate(need(length(oriData()$dat) > 1, message="Only one transcript in dataset."))
        qplot(oriData()$dat_tw$avg_cov, geom="histogram") + xlab("Average per-nucleotide coverage") + 
          ylab("Number of transcripts") +ggtitle("Distribution of transcripts by average coverage")
      })
      
      
      
      #------------------------------------------------------------------
      
      #------------------------------------------------------------------
      #Download Analyses
      #------------------------------------------------------------------
      
      observeEvent(input$save_all, {
        
        output$dwnldAllUI <- renderUI({
          actionButton(inputId= "dwnldAll", label= "Begin download")
        })
        
        observe({
          dir <- file.path(getwd(),"Saved_results", file_path_sans_ext(input[["file"]][1, "name"]), input$selected_transcript, "/")
          
          canCreateFile <- FALSE
          if (!file.exists(dir)) {
            canCreateFile <- dir.create(dir, showWarnings= FALSE, recursive =TRUE)
          }
          
          if (canCreateFile | file.exists(dir)) {
            
            output$winSizeUI <- renderUI({
              numericInput(inputId ="winSize", label="Enter length of transcript to visulaize in one panel of plot",
                           value=100, min= 50, max=nrow(data()), step=1)
            })
            
            #Begin download
            
            observeEvent(input$dwnldAll, {
              
              #------------------------------------------------------------------
              #pairwise replicate comparisons
              #------------------------------------------------------------------
              reactivities <- reactivityCalc()[-1]
              if (input$replicate_count > 1) {
                pairwiseSNR <- matrix(, input$replicate_count, input$replicate_count)
                colnames(pairwiseSNR) <- paste("Replicate", 1:input$replicate_count, sep="_")
                rownames(pairwiseSNR) <- paste("Replicate", 1:input$replicate_count, sep="_")
                for (i in 1:input$replicate_count) {
                  for (j in 1:input$replicate_count) {
                    s <- apply(reactivities[, c(i,j)], 1, mean, na.rm=TRUE)/apply(reactivities[,c(i,j)], 1, sd, na.rm=TRUE)
                    s[s > clip_at] <- clip_at
                    pairwiseSNR[i,j] <- if(i==j) 1 else mean(s, na.rm=TRUE)
                  }
                }
                
                setEPS()
                postscript(file= file.path(dir, "pairwise_mean_SNR.eps"))
                corrplot.mixed(pairwiseSNR,lower="number", upper="ellipse",is.corr=FALSE)
                dev.off()
                
                setEPS()
                postscript(file= file.path(dir, "correlation_matrix.eps"))
                corrplot.mixed(cor(reactivities, use = "pairwise.complete.obs"),lower="number", upper="ellipse")
                dev.off()
              }
              
              
              #------------------------------------------------------------------
              #reactivities with error bars
              #------------------------------------------------------------------
              rea_reshaped <- list(mean=NA, sdev=data.frame(NA) )
              rea_plot <- list(means= NULL, error_bars= NULL)
              
              
              rea_reshaped$mean = melt(reactivityCalc(), id.vars="Nucleotide", variable.name= "Replicate")
              
              if (input$error_bars) {
                if (input$error_method == "Bootstrap") {
                  rea_boot <- getReactivities(boot(), est.flag(), nor.flag(), logar.flag(), zeroing.flag(), !input$local_coverage)
                }
                
                rea_reshaped$sdev <- if (input$error_method == "Bootstrap") sd.calc(rea_boot) else data.frame(sdev= melt(sdCalc(), id.vars="Nucleotide")$value)
                rea_plot$error_bars <- geom_errorbar(aes(ymin= value-sdev, ymax= value+sdev), colour= if (input$sequence_info) "black" else"red", width=0.5)
              }
              
              data1_plot = if (input$error_bars & if (nrow(rea_reshaped$sdev) >1) nrow(rea_reshaped$mean) == nrow(rea_reshaped$sdev) else FALSE ) data.frame(rea_reshaped$mean, sdev=rea_reshaped$sdev) else rea_reshaped$mean
              if (input$sequence_info) data1_plot$sequence <- rep(data()$sequence, input$replicate_count)
              start_point <- 1
              end_point <- nrow(data())
              
              ticks <- floor((end_point- start_point+1)/10)
              
              if (input$sequence_info) {
                rea_plot$means <- ggplot(data1_plot, aes(x= Nucleotide, y= value, fill= sequence), environment=environment())+
                  geom_bar(stat="identity")  +facet_grid( Replicate ~ ., scale="free_y") +
                  ggtitle("SHAPE reactivities") + theme(plot.title= element_text(face="bold"))+
                  xlab("Nucleotide position")+ylab("Reactivity")
              } else {
                rea_plot$means <- ggplot(data1_plot, aes(x= Nucleotide, y= value), environment=environment())+ 
                  geom_bar(stat="identity")  +facet_grid( Replicate ~ ., scale="free_y") +
                  ggtitle("SHAPE reactivities") + theme(plot.title= element_text(face="bold"))+ylim(0, 4) +
                  xlab("Nucleotide position")+ylab("Reactivity")
              }
              
              setEPS()
              postscript(file= file.path(dir, "reactivity_complete.eps"))
              print(rea_plot$means + scale_x_discrete(breaks=c(start_point, start_point + (ticks*c(1:9)) ))+ 
                      coord_cartesian(xlim = c(start_point, end_point)) + if (input$error_bars) rea_plot$error_bars)
              dev.off()
              
              end_point <- 0
              
              if (input$sequence_info) {
                rea_plot$means <- ggplot(data1_plot, aes(x= Nucleotide, y= value, fill= sequence), environment=environment())+ geom_bar(stat="identity")  +facet_grid( Replicate ~ ., scale="free_y") +
                  ggtitle("SHAPE reactivities") + theme(plot.title= element_text(face="bold"))+xlab("Nucleotide position")+ylab("Reactivity")
              } else {
                rea_plot$means <- ggplot(data1_plot, aes(x= Nucleotide, y= value), environment=environment())+
                  geom_bar(stat="identity")  +facet_grid( Replicate ~ ., scale="free_y") +
                  ggtitle("SHAPE reactivities") + theme(plot.title= element_text(face="bold"))+xlab("Nucleotide position")+ylab("Reactivity")
              }
              
              for (i in 1:round(nrow(data())/ input$winSize)) {
                start_point <- end_point + 1
                end_point <- if (i*input$winSize < nrow(data()) & i == round(nrow(data())/ input$winSize) )  nrow(data()) else i*input$winSize
                
                ticks <- floor((end_point- start_point+1)/10)
                setEPS()
                postscript(file= file.path(dir, paste("reactivity_", i, ".eps", sep="")))
                print(rea_plot$means + scale_x_discrete(breaks=c(start_point,start_point + (ticks*c(1:9)) ))+coord_cartesian(xlim = c(start_point, end_point)) + if (input$error_bars) rea_plot$error_bars)
                dev.off()
                
              }
              
              observe({
                
                if (input$local_coverage) {
                  
                  #------------------------------------------------------------------
                  #Coverage
                  #------------------------------------------------------------------
                  cov_minus <- data.frame(data()[grep(pat="minus.pass", colnames(data()))])
                  colnames(cov_minus) <- paste("Replicate", 1:input$replicate_count, sep="_")
                  cov_minus$Nucleotide <- c(1:nrow(data()))
                  cov_minus <- melt(cov_minus, id.vars="Nucleotide", variable.name= "Replicate")
                  
                  if (input$sequence_info) {
                    cov_minus$sequence <- rep(data()$sequence, input$replicate_count)
                    setEPS()
                    postscript(file= file.path(dir, "local_coverage_minus.eps"))
                    print(ggplot(cov_minus, aes(x= Nucleotide, y= value, fill= sequence), environment=environment())+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete()) + facet_grid( Replicate ~ ., scale="free_y")+
                            ggtitle("(-) channel local coverage")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Local coverage"))
                    dev.off()
                  } else {
                    setEPS()
                    postscript(file= file.path(dir, "local_coverage_minus.eps"))
                    print(ggplot(cov_minus, aes(x= Nucleotide, y= value), environment=environment())+ geom_line(stat="identity") + scale_x_discrete(breaks=ticks_complete()) + facet_grid( Replicate ~ ., scale="free_y")+
                            ggtitle("(-) channel local coverage")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Local coverage"))
                    dev.off()
                  }
                  
                  cov_plus <- data.frame(data()[grep(pat="plus.pass", colnames(data()))])
                  colnames(cov_plus) <- paste("Replicate", 1:input$replicate_count, sep="_")
                  cov_plus$Nucleotide <- c(1:nrow(data()))
                  cov_plus <- melt(cov_plus, id.vars="Nucleotide", variable.name= "Replicate")
                  
                  if (input$sequence_info) {
                    cov_plus$sequence <- rep(data()$sequence, input$replicate_count)
                    setEPS()
                    postscript(file= file.path(dir, "local_coverage_plus.eps"))
                    print(ggplot(cov_plus, aes(x= Nucleotide, y= value, fill= sequence), environment=environment())+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete()) + facet_grid( Replicate ~ ., scale="free_y")+
                            ggtitle("(+) channel local coverage")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Local coverage"))
                    dev.off()
                  } else {
                    setEPS()
                    postscript(file= file.path(dir, "local_coverage_plus.eps"))
                    print(ggplot(cov_plus, aes(x= Nucleotide, y= value), environment=environment())+ geom_line(stat="identity") + scale_x_discrete(breaks=ticks_complete()) + facet_grid( Replicate ~ ., scale="free_y")+
                            ggtitle("(+) channel local coverage")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Local coverage"))
                    dev.off()
                  }
                }
              })
              
              #------------------------------------------------------------------
              #Terminations
              #------------------------------------------------------------------
              term_minus <- data.frame(data()[grep(pat="minus.counts", colnames(data()))])
              colnames(term_minus) <- paste("Replicate", 1:input$replicate_count, sep="_")
              term_minus$Nucleotide <- c(1:nrow(data()))
              term_minus <- melt(term_minus, id.vars="Nucleotide", variable.name= "Replicate")
              
              if (input$sequence_info) {
                term_minus$sequence <- rep(data()$sequence, input$replicate_count)
                setEPS()
                postscript(file= file.path(dir, "modifications_minus.eps"))
                print(ggplot(term_minus, aes(x= Nucleotide, y= value, fill= sequence), environment=environment())+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete()) + facet_grid( Replicate ~ ., scale="free_y")+
                        ggtitle("(-) channel modification counts")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Number of modifications"))
                dev.off()
              } else {
                setEPS()
                postscript(file= file.path(dir, "modifications_minus.eps"))
                print(ggplot(term_minus, aes(x= Nucleotide, y= value), environment=environment())+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete()) + facet_grid( Replicate ~ ., scale="free_y")+
                        ggtitle("(-) channel modification counts")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Number of modifications"))
                dev.off()
              }
              
              term_plus <- data.frame(data()[grep(pat="plus.counts", colnames(data()))])
              colnames(term_plus) <- paste("Replicate", 1:input$replicate_count, sep="_")
              term_plus$Nucleotide <- c(1:nrow(data()))
              term_plus <- melt(term_plus, id.vars="Nucleotide", variable.name= "Replicate")
              
              if (input$sequence_info) {
                term_plus$sequence <- rep(data()$sequence, input$replicate_count)
                setEPS()
                postscript(file= file.path(dir, "modifications_plus.eps"))
                print(ggplot(term_plus, aes(x= Nucleotide, y= value, fill= sequence), environment=environment())+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete()) + facet_grid( Replicate ~ ., scale="free_y")+
                        ggtitle("(-) channel modification counts")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Number of modifications"))
                dev.off()
              } else {
                setEPS()
                postscript(file= file.path(dir, "modifications_plus.eps"))
                print(ggplot(term_plus, aes(x= Nucleotide, y= value), environment=environment())+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete()) + facet_grid( Replicate ~ ., scale="free_y")+
                        ggtitle("(-) channel modification counts")+ theme(plot.title= element_text(face="bold"))+ xlab("Nucleotide position") + ylab("Number of modifications"))
                dev.off()
              }
              
              #------------------------------------------------------------------
              #Replicate SNR comparison
              #------------------------------------------------------------------
              reactivities <- reactivityCalc()[-1]
              if (input$replicate_count > 1) {
                s <- apply(reactivities, 1, mean, na.rm=TRUE)/apply(reactivities, 1, sd, na.rm=TRUE)
                s[s > clip_at] <- clip_at
                SNR <- data.frame(Nucleotide= 1:nrow(reactivities), value=s)
                if (input$sequence_info) {
                  SNR$sequence <- data()$sequence
                  setEPS()
                  postscript(file= file.path(dir, "SNR_comparing_replicates.eps"))
                  print(ggplot(SNR, aes(x= Nucleotide, y= value, fill=sequence), environment=environment())+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete())+
                          ggtitle("Signal-to-Noise Ratio (comparing replicates)") + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value"))
                  dev.off()
                } else {
                  setEPS()
                  postscript(file= file.path(dir, "SNR_comparing_replicates.eps"))
                  print(ggplot(SNR, aes(x= Nucleotide, y= value), environment=environment())+ geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete())+
                          ggtitle("Signal-to-Noise Ratio (comparing replicates)") + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value"))
                  dev.off()
                }
              }
              
              
              #------------------------------------------------------------------
              #Individual replicate SNR
              #------------------------------------------------------------------
              SNR <- melt(getSNR$values, id.vars="Nucleotide", variable.name= "Replicate")
              window_SNR <- melt(SNR.scanner(getSNR$values, input$win_size), id.vars="Nucleotide", variable.name= "Replicate")
              if (input$sequence_info) {
                SNR$sequence <- rep(data()$sequence, input$replicate_count)
                window_SNR$sequence <- "A"
                if (input$window_with_ind) {
                  setEPS()
                  postscript(file= file.path(dir, "individual_replicate_SNR.eps"))
                  print(ggplot(SNR, aes(x= Nucleotide, y= value, fill= sequence), environment=environment()) + geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete())+ facet_grid( Replicate ~ ., scale="free_y")+
                          ggtitle(paste("Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value") + geom_line(data=window_SNR, aes(x= Nucleotide, y= value)) + guides(colour=FALSE))
                  dev.off()
                } else {
                  setEPS()
                  postscript(file= file.path(dir, "individual_replicate_SNR.eps"))
                  print(ggplot(SNR, aes(x= Nucleotide, y= value, fill= sequence), environment=environment()) + geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete())+ facet_grid( Replicate ~ ., scale="free_y")+
                          ggtitle(paste("Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value"))
                  dev.off()
                }
                
              } else {
                if (input$window_with_ind) {
                  setEPS()
                  postscript(file= file.path(dir, "individual_replicate_SNR.eps"))
                  print(ggplot(SNR, aes(x= Nucleotide, y= value), environment=environment()) + geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete()) + facet_grid( Replicate ~ ., scale="free_y")+
                          ggtitle(paste("Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value") + geom_line(data=melt(SNR.scanner(getSNR$values, input$win_size), id.vars="Nucleotide", variable.name= "Replicate"), aes(x= Nucleotide, y= value, colour="Red")) + guides(colour=FALSE))
                  dev.off()
                } else {
                  setEPS()
                  postscript(file= file.path(dir, "individual_replicate_SNR.eps"))
                  print(ggplot(SNR, aes(x= Nucleotide, y= value), environment=environment()) + geom_bar(stat="identity") + scale_x_discrete(breaks=ticks_complete())+ facet_grid( Replicate ~ ., scale="free_y")+
                          ggtitle(paste("Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value"))
                  dev.off()
                }
              }
              
              #------------------------------------------------------------------
              #Window level SNR
              #------------------------------------------------------------------
              window_SNR <- melt(SNR.scanner(getSNR$values, input$win_size), id.vars="Nucleotide", variable.name= "Replicate")
              setEPS()
              postscript(file= file.path(dir, "window_level_SNR.eps"))
              print(ggplot(window_SNR, aes(x= Nucleotide, y= value), environment=environment())+ geom_line(aes(group=Replicate, colour= Replicate)) + scale_x_discrete(breaks=ticks_complete())+scale_colour_brewer(palette="Set1")+
                      ggtitle(paste("Window level Signal-to-Noise Ratio from ", if (input$error_method == "Bootstrap") "bootstrap" else "formula", sep="")) + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("SNR value"))
              dev.off()
              
              #------------------------------------------------------------------
              #Reactivity distribution
              #------------------------------------------------------------------
              reactivities <- melt(reactivityCalc(), id.vars="Nucleotide", variable.name= "Replicate")
              medians <-  ddply(reactivities, .(Replicate), function(.d) data.frame(x=median(.d$value, na.rm=TRUE)))
              means <-  ddply(reactivities, .(Replicate), function(.d) data.frame(x=mean(.d$value, na.rm=TRUE)))
              setEPS()
              postscript(file= file.path(dir, "reactivity_distribution.eps"))
              print(ggplot(data=reactivities, aes(x=value), environment=environment()) + geom_histogram() + facet_grid( Replicate ~ ., scale="free_y") + ggtitle("Reactivity score distributions") + theme(plot.title= element_text(face="bold"))+ xlab("Reactivity score")+ 
                      geom_vline(data=medians, aes(xintercept=x, colour="Median"), show_guide=TRUE) + geom_vline(data=means, aes(xintercept=x, colour="Mean"), show_guide=TRUE) +scale_colour_manual(name="Statistic", values=c("Median"= "red", "Mean"="green")))
              dev.off()
              
              #------------------------------------------------------------------
              #Modification counts at ACGU
              #------------------------------------------------------------------
              if (input$sequence_info) {
                term_minus <- data.frame(data()[c(grep(pat="minus.counts", colnames(data())), grep(pat="sequence", colnames(data())))], channel= "minus", Nucleotide= 1:nrow(data()) )
                names(term_minus)[grep(pat="minus.counts", colnames(term_minus))] <- paste("Replicate", file_ext(colnames(term_minus)[grep(pat="minus.counts", colnames(term_minus))]), sep="_")
                term_plus <- data.frame(data()[c(grep(pat="plus.counts", colnames(data())), grep(pat="sequence", colnames(data())))], channel= "plus", Nucleotide= 1:nrow(data()))
                names(term_plus)[grep(pat="plus.counts", colnames(term_plus))] <- paste("Replicate", file_ext(colnames(term_plus)[grep(pat="plus.counts", colnames(term_plus))]), sep="_")
                term <- melt(rbind(term_minus, term_plus), id.vars=c("Nucleotide", "sequence", "channel"), variable.name="Replicate")
                sum_term_1 <- data.frame(ddply(term, .(sequence, channel), function(.d) data.frame(Counts=sum(.d$value))), Replicate="All replicates")
                sum_term_2 <- ddply(term, .(sequence, Replicate, channel), function(.d) data.frame(Counts=sum(.d$value)))
                sum_term <- rbind(sum_term_1, sum_term_2)
                setEPS()
                postscript(file= file.path(dir, "ACGU_modification_counts.eps"))
                print(ggplot(data=sum_term, aes(x= sequence, y= Counts), environment=environment()) +geom_bar(stat="identity") +facet_grid(Replicate ~ channel, scale="free_y"))           
                dev.off()
              }
              
              #------------------------------------------------------------------
              #Individual replicate SNR distributions
              #------------------------------------------------------------------
              if (input$replicate_count > 1) {
                SNR <- melt(getSNR$values, id.vars="Nucleotide", variable.name= "Replicate")
                medians <-  ddply(SNR, .(Replicate), function(.d) data.frame(x=median(.d$value, na.rm=TRUE)))
                means <-  ddply(SNR, .(Replicate), function(.d) data.frame(x=mean(.d$value, na.rm=TRUE)))
                setEPS()
                postscript(file= file.path(dir, "replicate_wise_SNR_distributions.eps"))
                print(ggplot(data=SNR, aes(x=value), environment=environment()) + geom_histogram() + facet_grid( Replicate ~ ., scale="free_y") + ggtitle(paste(if (input$error_method == "Bootstrap") "Bootstrap" else "Formula", "SNR distributions", sep="")) + theme(plot.title= element_text(face="bold"))+ xlab("SNR value")+ 
                        geom_vline(data=medians, aes(xintercept=x, colour="Median"), show_guide=TRUE) + geom_vline(data=means, aes(xintercept=x, colour="Mean"), show_guide=TRUE) +scale_colour_manual(name="Statistic", values=c("Median"= "red", "Mean"="green")))
                dev.off()
              }
              
              #------------------------------------------------------------------
              #Reactivities overlay
              #------------------------------------------------------------------
              reactivities <- melt(reactivityCalc(), id.vars="Nucleotide", variable.name= "Replicate")
              if (input$replicate_count > 1) {
                setEPS()
                postscript(file= file.path(dir, "reactivities_overlaid.eps"))
                print(ggplot(reactivities, aes(x= Nucleotide, y= value), environment=environment())+ geom_line(aes(group=Replicate, colour= Replicate)) + scale_x_discrete(breaks=ticks_complete())+scale_colour_brewer(palette="Set1")+
                        ggtitle("Reactivity scores from the replicates") + theme(plot.title= element_text(face="bold")) + xlab("Nucleotide position") + ylab("Reactivity score"))
                dev.off()
              }
              
              #------------------------------------------------------------------
              #Replicate SNR (comparison) distribution
              #------------------------------------------------------------------
              reactivities <- reactivityCalc()[-1]
              if (input$replicate_count > 1) {
                s <- apply(reactivities, 1, mean, na.rm=TRUE)/apply(reactivities, 1, sd, na.rm=TRUE)
                s[s > clip_at] <- clip_at
                SNR <- data.frame(Nucleotide= 1:nrow(reactivities), value=s)
                setEPS()
                postscript(file= file.path(dir, "distribution_of_SNR_from_replicate_comparison.eps"))
                print(ggplot(data=SNR, aes(x=value), environment=environment()) +geom_histogram()+ggtitle("Replicate SNR distribution") + theme(plot.title= element_text(face="bold"))+ xlab("SNR value")+
                        geom_vline(aes(xintercept=median(value, na.rm=TRUE), colour="Median"), show_guide=TRUE)+geom_vline(aes(xintercept=mean(value, na.rm=TRUE), colour="Mean"), show_guide=TRUE) + scale_colour_manual(name="Statistic", values=c("Median"= "red", "Mean"="green")) )
                dev.off()
              }
              
              #dwnldall_closure
            })
            
          } else {
            
            output$createFileMsgUI <- renderUI({
              conditionalPanel(condition= !canCreateFile,
                               textOutput('canCreateFileMsg'))
            })
            
            output$canCreateFileMsg <- renderText("Target directory not creatable.")
          }
        })
      })
      
    })
    
  })
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
  on.exit(registerDoSEQ())
})