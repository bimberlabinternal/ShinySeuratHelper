
# Shiny Seurat Helper (SSH)

library(shiny)
library(rclipboard)
library(shinydashboard)
library(ggplot2)
library(data.table)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(grid)
library(gridExtra) 
library(dplyr)
library(Seurat)
col_vector <- OOSAP::ColorTheme()$col_vector

#source("./Fxs.R")

# BS <- readRDS("./data/BigSeurat.rds")
# SS <- BS[, sample(colnames(BS), 300, replace = F)]
# saveRDS(SS, "./data/SmallSeurat.rds")

SerInputObj <- readRDS("./data/SmallSeurat.rds")





ui <- dashboardPage(skin="red",
  dashboardHeader(title = "ShinySeuratHelper",
                  dropdownMenu(type = "messages",
                               messageItem(
                                 from = "A",
                                 message = "Task 1 due."
                               ),
                               messageItem(
                                 from = "B",
                                 message = "question?",
                                 icon = icon("question"),
                                 time = "13:45"
                               ),
                               messageItem(
                                 from = "C",
                                 message = "Woohoo!",
                                 icon = icon("life-ring"),
                                 time = "2014-12-01"
                               )
                  ),
                  dropdownMenu(type = "notifications",
                               notificationItem(
                                 text = "ABC",
                                 icon("users")
                               ),
                               notificationItem(
                                 text = "XYZ",
                                 icon("rep"),
                                 status = "success"
                               ),
                               notificationItem(
                                 text = "IJK",
                                 icon = icon("exclamation-triangle"),
                                 status = "warning"
                               )
                  ),
                  #red, yellow, aqua, blue, light-blue, green, navy, teal, olive, lime, orange, fuchsia, purple, maroon, black.
                  dropdownMenu(type = "tasks", badgeStatus = "success",
                               taskItem(value = 99, color = "green",
                                        "Preprocessing"
                               ),
                               taskItem(value = 99, color = "aqua",
                                        "SDA processing"
                               ),
                               taskItem(value = 99, color = "yellow",
                                        "DE analysis"
                               ),
                               taskItem(value = 70, color = "olive",
                                        "Documentation"
                               ),
                               taskItem(value = 15, color = "black",
                                        "Manuscript"
                               ),
                               taskItem(value = 70, color = "red",
                                        "Overall project"
                               )
                  )
                  
  ),
  
  
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Main Tab", tabName = "MainDash", icon = icon("dashboard"),
               badgeLabel = "underconst.", badgeColor = "yellow"),
      menuItem("Tab 2. (dev)", tabName = "SecondDash", icon = icon("affiliatetheme"),
               badgeLabel = "soon", badgeColor = "red"),
      menuItem("Tab 3. (dev)", tabName = "ThirdDash", icon = icon("allergies"),
               badgeLabel = "soon", badgeColor = "red"),
      menuItem("Tab 4. (dev)", tabName = "FourthDash", icon = icon("arrows-alt"),
               badgeLabel = "soon", badgeColor = "red"),
      menuItem("Enrichment Analysis", tabName = "Enrichment", icon = icon("dashboard"),
               badgeLabel = "underconst.", badgeColor = "yellow"),
      menuItem("Source code", icon = icon("file-code-o"), 
               href = "https://")
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
                      .content-wrapper {
                      background-color: black !important;
                      }
                      .main-sidebar {
                      background-color: black !important;
                      }
                      "))),
    tabItems(
      # First tab content
      tabItem(tabName = "MainDash",
              fluidRow(
                
                box(
                  title = "UMAP", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("UMAP_2D"), 
                  width = 5, background = "black"
                ),
                box(
                  title = "tSNE", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("tSNE_2D"), 
                  width = 5, background = "black"
                ),
                box(
                  title = "PCA_C1vC2", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("PCA_12"), 
                  width = 5, background = "black"
                ),
                box(
                  title = "Legend", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("Legend"), 
                  width = 5, background = "black"
                ),
                box(
                  title = "PCA_HM", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("PCA_HM"), 
                  width = 10, background = "black"
                )
        
              ) 
      ),
      
      # SecondDash content
      tabItem(tabName = "SecondDash",
              h2("SecondDash content")
      ),
      
      # ThirdDash content
      tabItem(tabName = "ThirdDash",
              h2("ThirdDash content")
      ),
      
      # FourthDash tab content
      tabItem(tabName = "FourthDash",
              h2("FourthDash tab content")
      ),
      
      # Enrichment tab content
      tabItem(tabName = "Enrichment",
              fluidRow(

              )
      )
      
    )
  )
  
  
)


server <- function(input, output, session) {
  
  output$UMAP_2D <- renderPlot({
    DimPlot(SerInputObj)+ theme_bw() +
      theme(legend.position = "none", aspect.ratio=1,
            legend.title = element_blank())  +
      scale_color_manual(values=rev(col_vector)) + guides(colour = guide_legend(override.aes = list(size=2, alpha=1), ncol =2))  +
      coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)
    
  })
  
  output$tSNE_2D <- renderPlot({
    DimPlot(SerInputObj, reduction = "tsne")+ theme_bw() +
      theme(legend.position = "none", aspect.ratio=1,
            legend.title = element_blank())  +
      scale_color_manual(values=rev(col_vector)) + 
      guides(colour = guide_legend(override.aes = list(size=2, alpha=1), ncol =2))  +
      coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)
    
  })
  
  output$PCA_12 <- renderPlot({
    PCAPlot(SerInputObj) + theme_bw() +
      theme(legend.position = "none", aspect.ratio=1,
            legend.title = element_blank())  +
      scale_color_manual(values=rev(col_vector)) + 
      guides(colour = guide_legend(override.aes = list(size=2, alpha=1), ncol =2))  +
      coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)
    
  })
  
  output$Legend <- renderPlot({
    
    ggFph <- DimPlot(SerInputObj)+ theme_bw() +
      theme(legend.position = "bottom", aspect.ratio=1,
            legend.title = element_blank())  +
      scale_color_manual(values=rev(col_vector)) + 
      guides(colour = guide_legend(override.aes = list(size=2, alpha=1), ncol =2))  +
      coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)
    

    grid.draw(cowplot::get_legend(ggFph))
    
  })
  
  output$PCA_HM <- renderPlot({
    PCHeatmap(SerInputObj, dims = 1, nfeatures = 30, balanced = T) 
    #Seurat::TopFeatures(SerInputObj, dims = 1, nfeatures = 30, balanced = T)

    #Loadings(SerInputObj, reduction = "pca")
  })
  
  
}


shinyApp(ui, server)

