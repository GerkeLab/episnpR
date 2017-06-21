#.libPaths(c("~/Documents/Rpackages",.libPaths()))
library(shinydashboard)
library(shiny)
library(haploR)
###################################################################################################
ui <- dashboardPage(dashboardHeader(title="Linkage Disequilibrium"),
                    dashboardSidebar(sidebarMenu(id="tabs",
                                                 menuItem("SNPs",tabName = "tab1"),
                                                 menuItem("Info",tabName = "tab2"))),
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "tab1",
                                fluidRow(
                                  box(title = "SNP List",
                                      h5(helpText("Fill in SNP")),
                                      textInput("snpList","SNPs",value=""),
                                      h5(helpText("Upload SNP List (one SNP per line)")),
                                      fileInput("file1","Choose a file",
                                                accept=c('text/csv', 
                                                         'text/comma-separated-values,text/plain', 
                                                         '.csv')),
                                      tags$hr(),
                                      h5(helpText("Select file parameters:")),
                                      checkboxInput(inputId = 'header1', label= 'Header', value= TRUE),
                                      br(),
                                      radioButtons(inputId = 'sep1', label = 'Seperator', 
                                                   choices = c(Comma=',', Semicolon=';', Tab='\t', Space= ' '), 
                                                   selected= ','),
                                      actionButton("upload1", "Upload File")),
                                  box(title="Select Output",
                                      selectInput("pop","Population",c("EUR","AFR","AMR","ASN"), selected="EUR"),
                                      sliderInput("value","LD threshold",min=0,max=1,value=0.8),
                                      selectInput("parameters","Additional Output",c("chr","pos_hg38",                         
                                                                           "D'","is_query_snp",                       
                                                                           "ref","alt","AFR",                        
                                                                           "AMR","ASN","EUR",                        
                                                                           "GERP_cons","SiPhy_cons","Chromatin_States",           
                                                                           "Chromatin_States_Imputed","Chromatin_Marks","DNAse",                      
                                                                           "Proteins","eQTL","gwas",                       
                                                                           "grasp","Motifs","GENCODE_id",                 
                                                                           "GENCODE_name","GENCODE_direction","GENCODE_distance",           
                                                                           "RefSeq_id","RefSeq_name","RefSeq_direction",           
                                                                           "RefSeq_distance","dbSNP_functional_annotation"),
                                                  multiple=TRUE),
                                      actionButton("update1","Show Table")
                                      )
                                ),
                                fluidRow(
                                  column(12,align="center",offset=2,
                                  box(title = "LD",
                                      tableOutput("LDtable1"),
                                      h5(helpText("See link below for TAD information!")),
                                      uiOutput("hic1")))
                                )),
                        tabItem(tabName = "tab2",
                                fluidRow(
                                  box(title="App Details",
                                      h5(helpText("All of the data from this site come from various online resources such as the 1000 Genome Project and ENCODE.
                                                  To perform similar queries in R please check out the haploR package!")),
                                      h5(helpText("Contact: Jordan.H.Creed@moffitt.org")),
                                      h5(helpText("Website: http://travisgerke.com"))),
                                  box(title="Other resources",
                                      h5(helpText("Aiden Lab")),
                                      a("Juicebox", href="http://www.aidenlab.org/juicebox/", target="_blank"))
                                ))
                      ))
                    )

###################################################################################################
server <- function(input, output) {
  sample<-eventReactive(input$upload1,{
    samplefile<-input$file1
    if(is.null(samplefile)){return()}
    sample1<-read.table(file=samplefile$datapath, sep= input$sep1, header= input$header1,stringsAsFactors= FALSE)
  })

  dat<-eventReactive(input$update1,{
    if(input$snpList==""){
      dat<-sample()
      snps<-dat[,1]
      x<-queryHaploreg(query = snps,ldThresh = as.numeric(input$value), ldPop = input$pop)
      #x<-x[,input$parameters]
      return(x)
    }
    if(input$snpList!=""){
      snps<-as.character(unlist(strsplit(input$snpList,",")))
      x<-queryHaploreg(query = snps,ldThresh = input$value, ldPop = input$pop)
      return(x)
    }
  })
  
  snps<-eventReactive(input$update1,{
    if(input$snpList==""){
      dat<-sample()
      snps<-dat[,1]
      return(snps)
    }
    if(input$snpList!=""){
      snps<-as.character(unlist(strsplit(input$snpList,",")))
      return(snps)
    }
  })
  
  output$hic1<-renderUI({
    x<-snps()
    y<-dat()
    if (length(x)>1){
    a("Take me to HIC", href=paste0("http://promoter.bx.psu.edu/hi-c/view.php?species=human&assembly=hg19&source=inside&tissue=GM12878&type=Lieberman-raw&c_url=&transfer=&chr=chr",y$chr[1],"&start=",min(y$pos_hg38),"&end=",max(y$pos_hg38),"&sessionID=&browser=none"), target="_blank")
    } else if (length(x)==1){
      a("Take me to HIC", href=paste0("http://promoter.bx.psu.edu/hi-c/view.php?species=human&assembly=hg19&source=inside&tissue=GM12878&type=Lieberman-raw&resolution=25&c_url=&transfer=&gene=",x,"&sessionID=&browser=none"), target="_blank")
    }
    })
  
  output$LDtable1<-renderTable({
    x<-dat()
    #x<-x[,input$parameters]
    x[,c("query_snp_rsid","rsID","r2",input$parameters)]
    #return(x[,input$parameters])
  })
  
}

################################################################################################### 
shinyApp(ui = ui, server = server)

