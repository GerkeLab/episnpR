.libPaths(c("~/Documents/Rpackages",.libPaths()))
library(shinydashboard)
library(shiny)
library(haploR)
###################################################################################################
ui <- dashboardPage(dashboardHeader(title="episnpR"),
                    dashboardSidebar(sidebarMenu(id="tabs",
                                                 menuItem("SNPs",tabName = "tab1"),
                                                 menuItem("Info",tabName = "tab2"))),
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "tab1",
                                fluidRow(
                                  box(title = "Query SNPs",
                                      textInput("snpList","Enter SNP rsIDs (comma seperated)",value=""),
                                      h5(helpText("Upload SNP List (one SNP per line)")),
                                      fileInput("file1","Choose a file",
                                                accept=c('text/csv', 
                                                         'text/comma-separated-values,text/plain', 
                                                         '.csv')),
                                      tags$hr(),
                                      actionButton("update1", "Perform query")),
                                  box(title="Select Output",
                                      selectInput("pop","Population",c("EUR","AFR","AMR","ASN"), selected="EUR"),
                                      sliderInput("value","LD threshold",min=0,max=1,value=0.8),
                                      checkboxGroupInput("parameters","Additional Output",c("Chromosome"="chr","Position"="pos_hg38",                         
                                                                                            "D'"="D'","Query SNP"="is_query_snp",                       
                                                                          "Reference allele"="ref","Alternative allele"="alt","MAF(AFR)"="AFR",                        
                                                                          "MAF(AMR)"="AMR","MAF(ASN)"="ASN","MAF(EUR)"="EUR",                        
                                                                           "GERP_cons","SiPhy_cons","Chromatin_States",           
                                                                           "Chromatin_States_Imputed","Chromatin_Marks","DNAse",                      
                                                                           "Proteins","eQTL","gwas",                       
                                                                           "grasp","Motifs","GENCODE_id",                 
                                                                           "GENCODE_name","GENCODE_direction","GENCODE_distance",           
                                                                           "RefSeq_id","RefSeq_name","RefSeq_direction",           
                                                                           "RefSeq_distance","dbSNP_functional_annotation"),inline = TRUE)
                                      )
                                ),
                                fluidRow(
                                  column(12,align="center",offset=2,
                                  box(title = "LD",
                                      tableOutput("LDtable1"),
                                      h5(helpText("See link below for TAD information!")),
                                      uiOutput("hic1"),
                                      uiOutput("clinical1"),
                                      uiOutput("ucsc1")))
                                )),
                        tabItem(tabName = "tab2",
                                fluidRow(
                                  box(title="App Details",
                                      h5(helpText("LD is calculated from 1000 Genomes Phase 1 (http://www.internationalgenome.org), and queried from HaploReg (http://archive.broadinstitute.org/mammals/haploreg/haploreg.php)
                                                  To perform similar queries in R please check out the haploR package!
                                                  TAD visualization comes from the Yue Lab (http://promoter.bx.psu.edu/hi-c/)"))),
                                  box(title="Development Team",
                                      h5(helpText("Programming: Jordan Creed, Travis Gerke")),
                                      h5(helpText("Scientific Input: Alvaro Monteiro")),
                                      h5(helpText("Website: http://travisgerke.com"))),
                                  box(title="Other resources",
                                      h5(helpText("Aiden Lab")),
                                      a("Juicebox", href="http://www.aidenlab.org/juicebox/", target="_blank"))
                                ))
                      ))
                    )

###################################################################################################
server <- function(input, output) {
  sample<-eventReactive(input$update1,{
    samplefile<-input$file1
    if(is.null(samplefile)){return()}
    sample1<-read.table(file=samplefile$datapath, sep= "\t", header= FALSE,stringsAsFactors= FALSE)
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
  
  output$clinical1<-renderUI({
    x<-snps()
    y<-dat()
    a("Take me to ClinVar", href=paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=",y$chr[1],"%5Bchr%5D+AND+",min(y$pos_hg38),"%3A",max(y$pos_hg38),"%5Bchrpos37%5D"), target="_blank")
  })
  
  output$ucsc1<-renderUI({
    x<-snps()
    y<-dat()
    a("Take me to Genome Browser", href=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr",y$chr[1],"%3A",min(y$pos_hg38),"%2D",max(y$pos_hg38),"&hgsid=596717155_di6qMTAMSs8fhJcRiuqsjlcsIxKA"), target="_blank")
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

