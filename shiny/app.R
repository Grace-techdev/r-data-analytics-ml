# --- Libraries ---
options(warn = -1)
library(shiny)
library(shinyWidgets)
library(ggplot2)
library(rpart)
library(rpart.plot)
library('ROCR')
library(ROCit)
library('class')
library(fpc)
library(RColorBrewer)
library('grDevices')
library(dplyr)

# paths
processed_dir    <- file.path("..", "data", "processed")
model_inputs_dir <- file.path("..", "data", "model_inputs")

# processed data
df_shiny <- read.csv(file.path(processed_dir, "dShiny.csv"))
df_shiny$Year <- factor(df_shiny$Year)

df_shiny_long <- read.csv(file.path(processed_dir, "dShiny_long.csv"))
df_shiny_long$Year <- factor(df_shiny_long$Year)

df_wide <- read.csv(
  file.path(processed_dir, "death_original.csv"),
  stringsAsFactors = FALSE
)

# model inputs
dTrain <- read.csv(
  file.path(model_inputs_dir, "dTrain.csv"),
  stringsAsFactors = FALSE
)

dCal <- read.csv(
  file.path(model_inputs_dir, "dCal.csv"),
  stringsAsFactors = FALSE
)

dTest <- read.csv(
  file.path(model_inputs_dir, "dTest.csv"),
  stringsAsFactors = FALSE
)

df_new_2 <- read.csv(
  file.path(model_inputs_dir, "dTransformed.csv"),
  stringsAsFactors = FALSE
)

# --- Functions and variables ---
behavioral_vars <- c("Breastfeeding", "Child.growth", "Low.birth.weight", "Iron.deficiency", 
                     "Vitamin.A.deficiency", "Drug.use", "Alcohol.use", "Unsafe.sex", "Low.physical.activity", 
                     "Tobacco", "Dietary")
environmental_vars <- c("Unsafe.water", "Air.pollution")
metabolic_vars <- c("High.systolic.blood.pressure", "High.fasting.plasma.glucose", "High.body.mass.index", "Low.bone.mineral.density")

numeric_var <- df_shiny %>% select(where(is.numeric)) %>% select(-c(Population, Death_Rate_Normalized, GDP, Behavioral, Environmental, Metabolic))
outcome <- 'Death_Rate_Cat'
pos <- '1'
vars <- setdiff(colnames(dTrain), c(outcome,'rgroup','Country','Year'))
numericVars <- vars[sapply(dTrain[,vars], class) %in% c('numeric','integer')]
factors <- numericVars
d <- subset(df_new_2, select = -c(GDP, Population, Total_Deaths, Death_Rate, Death_Rate_Normalized))
d1 <- subset(df_wide, select = c(numericVars, "Country"))
df_sum <- aggregate(. ~ Country, data = d1, FUN = sum)
log_trans <- function(df, vars) {
  for (var in vars) {
    df[,var]  <- log(df[,var] + 1)
  }
  return(df)
}
df <- log_trans(df_sum, numericVars)
scaled_df <- scale(df[,numericVars])
princ <- prcomp(scaled_df) 
nComp <- 2  
project2D <- as.data.frame(predict(princ, newdata=scaled_df)[,1:nComp])
nK <- 1
calcAUC <- function(ypred, ytrue) {
  perf <- performance(prediction(ypred, ytrue), 'auc')
  as.numeric(perf@y.values)
}
logLikelihood <- function(ypred, ytrue,epsilon=1e-6) {
  sum(ifelse(ytrue, log(ypred+epsilon), log(1-ypred+epsilon)), na.rm=T)
}
logNull <- logLikelihood(sum(dCal[,outcome]==pos)/nrow(dCal), dCal[,outcome]==pos)
plot_roc <- function(predcol, outcol, colour_id, label, overlaid=F) {
  if(length(predcol) == 0 || length(outcol) == 0) {
    return()  
  }
  ROCit_obj <- rocit(score = predcol, class = outcol == "1")
  par(new = overlaid)
  plot(ROCit_obj, col = c(colour_id, "black"), main = "", xlab = "", ylab = "",
       legend = FALSE, YIndex = FALSE, values = FALSE)
  lines(ROCit_obj$fpr, ROCit_obj$tpr, col = colour_id, lwd = 2, type = "b")
}
calcAUD_Deviance <- function(var, dTrain, dCal, dTest, outcome, pos, logNull) {
  pi <- paste('pred', var, sep = '')
  devDrop <- 2 * (logLikelihood(dCal[,pi], dCal[, outcome] == pos) - logNull)
  data.frame(
    Variable = var,
    AUC_Train = round(calcAUC(dTrain[,pi], dTrain$Death_Rate_Cat), 3),
    AUC_Cal = round(calcAUC(dCal[,pi], dCal$Death_Rate_Cat), 3),
    AUC_Test = round(calcAUC(dTest[,pi], dTest$Death_Rate_Cat), 3),
    Deviance_Reduction = round(devDrop, 3)
  )
}
mkPredC <- function(outCol, varCol, appCol) {
  pPos <- sum(outCol==pos)/length(outCol)
  naTab <- table(as.factor(outCol[is.na(varCol)]))
  pPosWna <- (naTab/sum(naTab))[pos]
  vTab <- table(as.factor(outCol), varCol)
  pPosWv <- (vTab[pos,]+1.0e-3*pPos)/(colSums(vTab)+1.0e-3)
  pred <- pPosWv[appCol]
  pred[is.na(appCol)] <- pPosWna
  pred[is.na(pred)] <- pPos
  pred
}
mkPredN <- function(outCol, varCol, appCol) {
  cuts <- unique(as.numeric(
    quantile(varCol, probs=seq(0, 1, 0.1), na.rm=T)))
  varC <- cut(varCol, cuts)
  appC <- cut(appCol, cuts)
  mkPredC(outCol, varC, appC)
}
for (v in numericVars) {
  pi <- paste('pred', v, sep='')
  dTrain[,pi] <- mkPredN(dTrain[,outcome], dTrain[,v], dTrain[,v])
  dTest[,pi] <- mkPredN(dTrain[,outcome], dTrain[,v], dTest[,v])
  dCal[,pi] <- mkPredN(dTrain[,outcome], dTrain[,v], dCal[,v])
}
selNumVars <- c()
minDrop <- 25 
for (v in numericVars) {
  pi <- paste('pred', v, sep='')
  devDrop <- 2*(logLikelihood(dCal[,pi], dCal[,outcome]==pos) - logNull)
  if (devDrop >= minDrop) {
    selNumVars <- c(selNumVars, pi)
  }
}
pca_result <- prcomp(d[, factors], scale. = TRUE)
loadings_matrix <- pca_result$rotation
important_features <- c()
cumulative_threshold <- 0.9
for (i in 1:length(loadings_matrix)) {
  cumulative_variance <- sum((pca_result$sdev[1:i])^2) / sum((pca_result$sdev)^2)
  if (cumulative_variance > cumulative_threshold) break
  pc <- paste0("PC", i)
  important_features <- unique(c(important_features, rownames(loadings_matrix)[abs(loadings_matrix[, pc]) > 0.25]))
}
vars1 <- numericVars
vars2 <- selNumVars
vars3 <- important_features
tmodel <- rpart(as.formula(paste("Death_Rate_Cat >0 ~", paste(vars1, collapse = " + "))), data = dTrain)
tmodel2 <- rpart(as.formula(paste("Death_Rate_Cat >0 ~", paste(vars2, collapse = " + "))), data = dTrain)
tmodel3 <- rpart(as.formula(paste("Death_Rate_Cat >0 ~", paste(vars3, collapse = " + "))), data = dTrain)
knnPredict <- function(df, knnTrain, knnCl, k) {
  knnDecision <- knn(knnTrain, df, knnCl, k=k, prob=TRUE)
  ifelse(knnDecision == TRUE,
         attributes(knnDecision)$prob,
         1 - attributes(knnDecision)$prob)
}
performanceMeasures <- function(ytrue, ypred, model.name = "model", threshold=0.5) {
  auc <- calcAUC(ypred, ytrue)
  dev.norm <- -2 * logLikelihood(ytrue, ypred)/length(ypred)
  cmat <- table(actual = ytrue, predicted = ypred >= threshold)
  accuracy <- sum(diag(cmat)) / sum(cmat)
  precision <- cmat[2, 2] / sum(cmat[, 2])
  recall <- cmat[2, 2] / sum(cmat[2, ])
  f1 <- 2 * precision * recall / (precision + recall)
  list(perf = data.frame(Model = model.name, AUC = auc, Precision = precision,
                         Recall = recall, F1 = f1, Dev_Norm = dev.norm),
       confusion_matrix = cmat)
}
compare_perf_table <- function(model1, x, y, knnTrain, knnCl, nK, model1_name="Model 1", threshold=0.5) {
  pred1 <- predict(model1, newdata=x)
  pred2 <- knnPredict(x, knnTrain, knnCl, nK)
  perf1 <- performanceMeasures(y, pred1, model.name="Decision Tree", threshold=threshold)
  perf2 <- performanceMeasures(y, pred2, model.name="KNN Model", threshold=threshold)
  perftable <- rbind(perf1$perf, perf2$perf)
  perftable
}
find_convex_hull <- function(proj2Ddf, groups) {
  do.call(rbind,
          lapply(unique(groups),
                 FUN = function(c) {
                   f <- subset(proj2Ddf, cluster==c);
                   f[chull(f),]
                 }
          )
  )
}
knnPredict <- function(df, knnTrain, knnCl, k) {
  knnDecision <- knn(knnTrain, df, knnCl, k=k, prob=TRUE)
  ifelse(knnDecision == TRUE,
         attributes(knnDecision)$prob,
         1 - attributes(knnDecision)$prob)
}
format_numbers <- function(num) {
  sapply(num, function(x) {
    if (x >= 1e9) {
      paste(format(round(x / 1e9, 2), nsmall = 2), "B")
    } else if (x >= 1e6) {
      paste(format(round(x / 1e6, 2), nsmall = 2), "M")
    } else if (x >= 1e3) {
      paste(format(round(x / 1e3, 2), nsmall = 2), "K")
    } else {
      as.character(x)
    }
  })
}
colors = c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 11, name = "RdBu"))
plot_theme <- function() {
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 25)),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    panel.grid.major = element_line(color = "grey", linetype = "dashed", size = 0.2),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),  
    plot.margin = margin(t = 20, r = 20, b = 15, l = 20),
    axis.ticks.length = unit(7, "pt"),  
    axis.text = element_text(size = 11), 
    axis.title = element_text(size = 12),  
    axis.title.x = element_text(margin = margin(t = 14)),  
    axis.title.y = element_text(margin = margin(r = 14))   #
  )
}

# --- UI ---
ui <- fluidPage(
  titlePanel(
    tags$div(
      tags$h3("Countries and Death Causes Data", style = "margin: 0; background-color: #666666; color: white; padding: 10px; text-align: center;")
    )
  ),
  tabsetPanel(
    tabPanel("Exploratory Data Analytics",
             sidebarPanel(
               radioButtons("eda_section", "Choose Analysis Type:", choices = c("Country Analysis", "Cause Analysis", "Correlation Analysis"), selected = "Country Analysis"),
               conditionalPanel(
                 condition = "input.eda_section == 'Country Analysis'",
                 radioButtons("plot_type", "Choose Plot Type:", choices = c("BarPlot", "BoxPlot"), selected = "BarPlot"),
                 conditionalPanel(
                   condition = "input.plot_type == 'BarPlot'",
                   selectInput("catVar", "Select y-axis Attribute:", c("Continent", "Country"), selected = "Country"),
                   uiOutput("numVar_ui")
                 ),
                 conditionalPanel(
                   condition = "input.plot_type == 'BoxPlot'",
                   selectInput("catVar_1", "Select y-axis Attribute:", c("Country", "Continent"), selected = "Continent"),
                   selectInput("numVar_1", "Select x-axis Attribute:", 
                               choices = list(
                                 "Sum" = c("Total_Deaths", "Death_Rate"), 
                                 "Category" = c("Behavioral", "Environmental", "Metabolic"), 
                                 "Specific Cause" = setdiff(names(numeric_var), c("Total_Deaths", "Death_Rate", "Behavioral", "Environmental", "Metabolic"))
                               ), 
                               selected = "Total_Deaths"),
                   radioButtons(inputId = "y_trans_1", label = "Log Transformation", choices = c("Original", "Log"), selected = "Original")
                 )
               ),
               conditionalPanel(
                 condition = "input.eda_section == 'Cause Analysis'",
                 radioButtons("cause_or_category", "Choose Analysis Level:", choices = c("Specific Cause", "General Category"), selected = "Specific Cause"),
                 selectInput("year_selection", "Filter by Year:", choices = c("All", levels(df_shiny$Year)), selected = "All")
               ),
               conditionalPanel(
                 condition = "input.eda_section == 'Correlation Analysis'",
                 selectInput("x_attr", "Select x-axis Attribute:", 
                             choices = list(
                               "Category" = c("Behavioral", "Environmental", "Metabolic"), 
                               "Specific Cause" = setdiff(names(numeric_var), c("Total_Deaths", "Death_Rate", "Behavioral", "Environmental", "Metabolic"))
                             ), 
                             selected = "Behavioral"),
                 selectInput("y_attr", "Select y-axis Attribute:", 
                             choices = list(
                               "Category" = c("Behavioral", "Environmental", "Metabolic"), 
                               "Specific Cause" = setdiff(names(numeric_var), c("Total_Deaths", "Death_Rate", "Behavioral", "Environmental", "Metabolic"))
                             ), 
                             selected = "Environmental"),
                 radioButtons(inputId = "transformation", label = "Log Transformation", choices = c("Original", "Log"), selected = "Log"),
                 checkboxInput("facet", "Facet Plot by Continent", FALSE)
               )
             ),
             mainPanel(
               conditionalPanel(
                 condition = "input.eda_section == 'Country Analysis'",
                 plotOutput("plot_by_country", height = "520px", width = "620px")
               ),
               conditionalPanel(
                 condition = "input.eda_section == 'Cause Analysis'",
                 plotOutput("plot_by_cause", height = "520px", width = "620px")
               ),
               conditionalPanel(
                 condition = "input.eda_section == 'Correlation Analysis'",
                 plotOutput("scatter_plot", height = "520px", width = "500px")
               )
             )
    ),
    tabPanel("Single Variable Models",
             sidebarLayout(
               sidebarPanel(
                 pickerInput("variable", "Select Variables for Comparison:", 
                             choices = numericVars,
                             selected = "High.systolic.blood.pressure",
                             multiple = TRUE,
                             options = list(`actions-box` = TRUE))
               ),
               mainPanel(
                 plotOutput("rocPlot", width = "500px", height = "500px"),
                 tableOutput("aucOutput")
               )
             )
    ),
    tabPanel("Classification",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("feature_set", "Choose Feature Set:",
                              choices = list("Feature1 (by Literature Review)" = "feature1", "Feature2 (by Single Model) " = "feature2", "Feature3 (by PCA)" = "feature3"),
                              selected = "feature1"),
                 tableOutput("selectedFeaturesTable"),
                 radioButtons("table_type", "Choose Evaluation Table:", 
                              choices = list("Performance Summary" = "performance", 
                                             "Confusion Matrix" = "confusion"), 
                              selected = "performance"),
               ),
               mainPanel(
                 plotOutput("classificationROCPlot", width = "500px", height = "500px"),
                 conditionalPanel(
                   condition = "input.table_type === 'confusion'",
                   fluidRow(
                     column(6, tags$div("Decision Tree", style = "font-size: 11px; font-weight: bold;"), tableOutput("treeConfusionMatrix")),
                     column(6, tags$div("KNN", style = "font-size: 11px; font-weight: bold;"), tableOutput("knnConfusionMatrix"))
                   )
                 ),
                 conditionalPanel(
                   condition = "input.table_type === 'performance'",
                   tableOutput("performanceTable")
                 )
               )
             )
    ),
    tabPanel("Clustering",
             sidebarLayout(
               sidebarPanel(
                 sliderInput("num_clusters", "Choose Number of Clusters (k):", min = 1, max = 6, value = 3),
                 tableOutput("cluster_avg_death_rate")
               ),
               mainPanel(
                 plotOutput("clusteringPlot", width = "550px", height = "500px")
               )
             )
    )
  ),
  tags$style(HTML(
    "table { font-size: 11px; }\n    th, td { padding: 1px; }"
  ))
)

# --- SERVER ---
server <- function(input, output) {
  options(warn = -1)
  
  # EDA PLOT
  ranges <- reactiveValues(x = NULL, y = NULL)
  output$numVar_ui <- renderUI({
    choices <- list(
      "Sum" = c("Total_Deaths", "Death_Rate"),
      "Category" = c("Behavioral", "Environmental", "Metabolic"),
      "Specific Cause" = setdiff(names(numeric_var), c("Total_Deaths", "Death_Rate", "Behavioral", "Environmental", "Metabolic"))
    )
    if (input$catVar %in% c("Cause Category", "Specific Cause")) {
      choices <- list("Sum" = c("Total_Deaths", "Death_Rate"))
    }
    selectInput("numVar", "Select x-axis Attribute:", choices = choices, selected = "Total_Deaths")
  })
  output$plot_by_country <- renderPlot({
    numVar <- input$numVar
    catVar <- input$catVar
    if (input$plot_type == "BarPlot") {
      df_plot <- aggregate(df_shiny[, numVar], by = list(df_shiny[, catVar]), FUN = sum)
      if (ncol(df_plot) == 2) { 
        names(df_plot) <- c(catVar, "y_attr")
      } else {
        return()
      }
      if (catVar == "Country") {
        names(df_plot) <- c("Country", "y_attr")
        if (nrow(df_plot) > 0) {
          df_plot <- df_plot[order(-df_plot$y_attr), ]
          df_plot <- head(df_plot, 20)
        } else {
          df_plot <- data.frame()
        }
        cols <- c("Africa" = "#1B9E77", "Europe" = "#7570B3", "Asia" = "#D95F02", 
                  "Americas" = "#E7298A", "Oceania" = "#66A61E")
        df_plot <- df_plot %>%
          left_join(df_shiny %>% distinct(Country, Continent), by = c("Country" = "Country"))
        ggplot(df_plot, aes(x = reorder(Country, y_attr), y = y_attr, fill = Continent)) +
          geom_col() +
          coord_flip() +
          scale_fill_manual(values = cols) +
          geom_text(aes(label = format_numbers(y_attr)), hjust = -0.1, size = 3) +
          labs(title = paste(numVar, "by", catVar), x = "", y = numVar, fill = "Continent") +
          scale_y_continuous(limits = c(0, max(df_plot$y_attr) * 1.1)) +
          plot_theme() +
          theme(legend.position = "bottom")
      } else if (catVar == "Continent") {
        names(df_plot) <- c("Continent", "y_attr")
        cols <- c("Africa" = "#1B9E77", "Europe" = "#7570B3", "Asia" = "#D95F02", 
                  "Americas" = "#E7298A", "Oceania" = "#66A61E")
        ggplot(df_plot, aes(x = reorder(Continent, y_attr), y = y_attr, fill = Continent)) +
          geom_col() +
          coord_flip() +
          scale_fill_manual(values = cols) +
          geom_text(aes(label = format_numbers(y_attr)), hjust = -0.1, size = 3) +
          labs(title = paste(numVar, "by", catVar), x = "", y = numVar) +
          scale_y_continuous(limits = c(0, max(df_plot$y_attr) * 1.1)) +
          plot_theme() +
          theme(legend.position = "bottom")
      }
    } else if (input$plot_type == "BoxPlot") {
      df_plot <- df_shiny
      catVar_1 <- input$catVar_1
      numVar_1 <- input$numVar_1
      if (catVar_1 == "Country") {
        set.seed(4009)
        sample_countires <- sample(unique(df_plot$Country), 20)
        df_plot <- subset(df_plot, Country %in% sample_countires)
      }
      if (input$y_trans_1 == "Log") {
        df_plot$y_attr <- log(df_plot[,numVar_1])
        df_plot <- subset(df_plot, !is.infinite(y_attr) & !is.nan(y_attr))
        y_attr <- "y_attr"
      } else {
        y_attr <- numVar_1
      }
      ggplot(df_plot, aes_string(x = catVar_1, y = y_attr, fill = "Continent")) +
        geom_boxplot() +
        coord_flip() +
        labs(x = catVar_1, y = numVar_1, title = paste(numVar_1, "by", catVar_1),x = "", y = numVar_1) +
        plot_theme() +
        theme(legend.position = "bottom") 
    }
  })
  output$plot_by_cause <- renderPlot({
    df_plot <- df_shiny_long
    year_selection <- input$year_selection
    cols <- c("Behavioral" = "#1B9E77", "Environmental" = "#7570B3", "Metabolic" = "#D95F02")
    if (year_selection != "All") {
      df_plot <- subset(df_plot, Year == year_selection)
    }
    if (input$cause_or_category == "Specific Cause") {
      df_plot <- aggregate(Deaths ~ Cause + Category, df_plot, sum)
      ggplot(df_plot, aes(x = reorder(Cause, Deaths), y = Deaths, fill = Category)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = cols) +
        geom_text(aes(label = format_numbers(Deaths)), hjust = -0.1, size = 3) +
        labs(title = paste("Total Deaths by Cause for Year:", year_selection), fill = "Category",x = "", y = "Total deaths") +
        scale_y_continuous(limits = c(0, max(df_plot$Deaths) * 1.1)) +
        plot_theme() +
        theme(legend.position = "bottom") 
    } else if (input$cause_or_category == "General Category") {
      df_plot <- aggregate(Deaths ~ Category, df_plot, sum)
      ggplot(df_plot, aes(x = reorder(Category, Deaths), y = Deaths, fill = Category)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = cols) +
        geom_text(aes(label = format_numbers(Deaths)), hjust = -0.1, size = 3) +
        labs(title = paste("Total Deaths by Category for Year:", year_selection),x = "", y = "Total deaths") +
        scale_y_continuous(limits = c(0, max(df_plot$Deaths) * 1.1)) +
        plot_theme() +
        theme(legend.position = "bottom") 
    }
  })
  output$scatter_plot <- renderPlot({
    df_plot <- df_shiny
    cols <- c("Africa" = "#1B9E77", "Europe" = "#7570B3", "Asia" = "#D95F02", 
              "Americas" = "#E7298A", "Oceania" = "#66A61E")
    if (input$transformation == "Log") {
      df_plot$x_attr <- log(df_plot[[input$x_attr]], base = exp(1))
      df_plot$y_attr <- log(df_plot[[input$y_attr]], base = exp(1))
    } else {
      df_plot$x_attr <- df_plot[[input$x_attr]]
      df_plot$y_attr <- df_plot[[input$y_attr]]
    }
    p <- ggplot(df_plot, aes(x = x_attr, y = y_attr, color = Continent)) +
      geom_point(alpha = 0.25) +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
      labs(title = paste(input$x_attr, "and", input$y_attr), x = input$x_attr, y = input$y_attr) +
      plot_theme() +
      scale_color_manual(values = cols) + 
      theme(legend.position = "bottom")
    if (input$facet) {
      p <- p + facet_wrap(~ Continent)
    }
    print(p)
  })

  # SINGLE VARIABLE MODELS 
  # rocPlot
  output$rocPlot <- renderPlot({
    selected_variable <- input$variable
    if (length(selected_variable) == 0) {
      plot.new()
      text(x = 0.5, y = 0.5, "Please choose at least one variable.", cex = 1.2)
      return()
    }
    n <- 1
    labels <- vector("character", length(selected_variable))
    cols <- vector("numeric", length(selected_variable))
    for (var in selected_variable) {
      pi <- paste('pred', var, sep = '')
      predictor <- dTest[, pi]
      plot_roc(predictor, dTest[, outcome], colour_id = colors[n], label = var, overlaid = TRUE)
      labels[n] <- var
      cols[n] <- colors[n]
      n <- n + 1
    }
    legend("bottomright", legend = labels, col = cols, lwd = 2)
    title("ROC Curve Comparison for Test Data")
  })
  
  # aucTable
  output$aucOutput <- renderTable({
    selected_variable <- input$variable
    auc_df <- data.frame() 
    for (var in selected_variable) {
      result_row <- calcAUD_Deviance(var, dTrain, dCal, dTest, outcome, pos, logNull)
      auc_df <- rbind(auc_df, result_row)  
    }
    auc_df <- auc_df[order(-auc_df$Deviance_Reduction), ]
    auc_df
  }, digits = 3)
  
  # CLASSIFICATION
  # selectedFeature
  output$selectedFeaturesTable <- renderTable({
    if (input$feature_set == "feature1") {
      feature_set <- vars1
    } else if (input$feature_set == "feature2") {
      feature_set <- vars2
    } else if (input$feature_set == "feature3") {
      feature_set <- vars3
    }
    data.frame(
      Features = paste(feature_set, collapse = ", "),
      Num = length(feature_set)
    )
  })
  # rocPlot
  output$classificationROCPlot <- renderPlot({
    feature_set <- input$feature_set
    if (feature_set == "feature1") {
      tmodel.pred <- predict(tmodel, newdata = dTest[vars1])
      knn.pred <- knnPredict(dTest[vars1], dTrain[vars1], dTrain[, outcome] == pos, nK)
    } else if (feature_set == "feature2") {
      tmodel.pred <- predict(tmodel2, newdata = dTest[vars2])
      knn.pred <- knnPredict(dTest[vars2], dTrain[vars2], dTrain[, outcome] == pos, nK)
    } else if (feature_set == "feature3") {
      tmodel.pred <- predict(tmodel3, newdata = dTest[vars3])
      knn.pred <- knnPredict(dTest[vars3], dTrain[vars3], dTrain[, outcome] == pos, nK)
    }
    dTest.gt <- dTest[, outcome] == pos
    cols <- brewer.pal(n = 4, name = "Dark2")
    roc_tmodel <- rocit(score = tmodel.pred, class = dTest.gt)
    roc_knn <- rocit(score = knn.pred, class = dTest.gt)
    plot(roc_tmodel, col = c(cols[1],"black"), lwd = 3, legend = FALSE, YIndex = FALSE, values = TRUE, asp = 1,
         cex.lab = 3, cex.axis = 2, cex.main = 3)
    lines(roc_knn$TPR ~ roc_knn$FPR, col = cols[2], lwd = 3)
    legend("bottomright", col = c("black", cols[1], cols[2]), legend = c("Null Model", "Decision Tree", "KNN"), 
           lwd = 2, lty = c(2, 1, 1))
    title("ROC Curve for Classification Models (Test Data)")
  })
  output$performanceTable <- renderTable({
    if (input$feature_set == "feature1") {
      perftable <- compare_perf_table(tmodel, dTest[vars1], dTest[, outcome] == pos, dTrain[, vars1], dTrain[, outcome] == pos, nK)
    } else if (input$feature_set == "feature2") {
      perftable <- compare_perf_table(tmodel2, dTest[vars2], dTest[, outcome] == pos, dTrain[, vars2], dTrain[, outcome] == pos, nK)
    } else if (input$feature_set == "feature3") {
      perftable <- compare_perf_table(tmodel3, dTest[vars3], dTest[, outcome] == pos, dTrain[, vars3], dTrain[, outcome] == pos, nK)
    }
    perftable  
  }, digits = 3)
  
  # confusion matrix
  output$treeConfusionMatrix <- renderTable({
      if (input$feature_set == "feature1") {
        pred_tree <- predict(tmodel, newdata = dTest[vars1])
      } else if (input$feature_set == "feature2") {
        pred_tree <- predict(tmodel2, newdata = dTest[vars2])
      } else if (input$feature_set == "feature3") {
        pred_tree <- predict(tmodel3, newdata = dTest[vars3])
      }
      actual <- dTest[, outcome] == pos
      cm_tree <- table(Actual = actual, Predicted = pred_tree >= 0.5)
      colnames(cm_tree) <- c("Pred Pos", "Pred Neg")
      rownames(cm_tree) <- c("Actual Pos", "Actual Neg")
      as.data.frame.matrix(cm_tree)
  }, rownames = TRUE)
  output$knnConfusionMatrix <- renderTable({
      if (input$feature_set == "feature1") {
        pred_knn <- knnPredict(dTest[vars1], dTrain[vars1], dTrain[, outcome] == pos, nK)
      } else if (input$feature_set == "feature2") {
        pred_knn <- knnPredict(dTest[vars2], dTrain[vars2], dTrain[, outcome] == pos, nK)
      } else if (input$feature_set == "feature3") {
        pred_knn <- knnPredict(dTest[vars3], dTrain[vars3], dTrain[, outcome] == pos, nK)
      }
      actual <- dTest[, outcome] == pos
      cm_knn <- table(Actual = actual, Predicted = pred_knn >= 0.5)
      colnames(cm_knn) <- c("Pred Pos", "Pred Neg")
      rownames(cm_knn) <- c("Actual Pos", "Actual Neg")
      as.data.frame.matrix(cm_knn)
  }, rownames = TRUE)
  
  # CLUSTERING
  output$clusteringPlot <- renderPlot({
    k <- input$num_clusters
    kmeans_result <- kmeans(scaled_df, k, nstart = 100, iter.max = 100)
    groups <- kmeans_result$cluster
    cluster_sizes <- table(groups)
    ordered_labels <- order(cluster_sizes, decreasing = TRUE)
    label_mapping <- setNames(seq_along(ordered_labels), ordered_labels)
    groups_mapped <- as.factor(label_mapping[as.character(groups)])
    kmclust.project2D <- cbind(project2D, cluster = groups_mapped, country = df$Country)
    kmclust.hull <- find_convex_hull(kmclust.project2D, groups_mapped)
    colors <- brewer.pal(n = k, name = "Dark2")
    centroids <- aggregate(cbind(PC1, PC2) ~ cluster, data = kmclust.project2D, mean)
    centroids$size <- cluster_sizes[as.numeric(levels(centroids$cluster))][centroids$cluster]
    shape_values <- c(16, 17, 18, 15, 3, 7)[1:k]
    ggplot(kmclust.project2D, aes(x = PC1, y = PC2, group = cluster)) +
      geom_point(aes(shape = cluster, color = cluster), alpha = 0.7) +
      geom_polygon(data = kmclust.hull, aes(fill = as.factor(cluster)), alpha = 0.4, linetype = 0) +
      scale_color_manual(values = colors, name = "Cluster") +
      scale_fill_manual(values = colors, name = "Cluster") +
      scale_shape_manual(values = shape_values, name = "Cluster") + 
      labs(title = sprintf("k-Means Clustering Visualization (k = %d)", k),
           x = "PC1",
           y = "PC2") +
      plot_theme() +
      theme(legend.position = "right") +
      geom_text(data = centroids, 
                aes(x = PC1, y = PC2, label = paste(size)), 
                vjust = -1, hjust = 0.5, size = 4, color = "black")
  })
}

shinyApp(ui = ui, server = server)