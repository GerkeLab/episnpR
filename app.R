.libPaths(c("~/Documents/Rpackages",.libPaths()))
library(shinydashboard)
library(shiny)
library(haploR)
library(data.table)
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
                                      textInput("snpList","Enter SNP rsIDs (comma separated)",value=""),
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
                                                                                            "Reference allele"="ref","Alternative allele"="alt","LD(AFR)"="AFR",                        
                                                                                            "LD(AMR)"="AMR","LD(ASN)"="ASN","LD(EUR)"="EUR",                        
                                                                                            "GERP scores"="GERP_cons","SiPhy scores"="SiPhy_cons","Chromatin States" ="Chromatin_States",           
                                                                                            "Imputed Chromatin States"="Chromatin_States_Imputed","Chromatin Marks"="Chromatin_Marks","DNAse"="DNAse",                      
                                                                                            "Proteins","eQTL","GWAS study name"="gwas",                       
                                                                                            "GRASP study name"="grasp","Motifs","GENCODE transcript ID"="GENCODE_id",                 
                                                                                            "GENCODE gene name"="GENCODE_name","GENCODE direction"="GENCODE_direction","GENCODE distance"="GENCODE_distance",           
                                                                                            "NCBI Reference Sequence Accession number"="RefSeq_id","NCBI Reference Sequence name"="RefSeq_name","NCBI Reference Sequence direction"="RefSeq_direction",           
                                                                                            "NCBI Reference Sequence distance"="RefSeq_distance","Annotated proteins"="dbSNP_functional_annotation"),inline = TRUE)
                                      )
                                ),
                                fluidRow(
                                  column(12,align="center",offset=2,
                                  tabBox(title = "Output",
                                    tabPanel("LD",
                                      tableOutput("LDtable1")),
                                    tabPanel("TADs",
                                      actionButton("tadButton","Look up TADs"),
                                      textOutput("tadBoundaries"),
                                      uiOutput("hic1")),
                                    tabPanel("Other",
                                      uiOutput("clinical1"),
                                      uiOutput("ucsc1"))))
                                )),
                        tabItem(tabName = "tab2",
                                fluidRow(
                                  box(title="App Details",
                                      h5(helpText("LD is calculated from 1000 Genomes Phase 1 (http://www.internationalgenome.org), and queried from HaploReg (http://archive.broadinstitute.org/mammals/haploreg/haploreg.php)
                                                  To perform similar queries in R please check out the haploR package!
                                                  TAD visualization comes from the Yue Lab (http://promoter.bx.psu.edu/hi-c/).
                                                  TAD locations are based off of those defined by Dixon et al in 'Topological domains in mammalian genomes identified by analysis of chromatin interactions'"
                                                  ))),
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
  
  in_tad<-eventReactive(input$tadButton,{
    load("./Data/tad_hg19.Rdata")
    snps<-snps()
    dat<-dat()
    dat<-dat[dat$rsID %in% snps,]
    snp_pos<-as.numeric(dat$pos_hg38)
    tad<-tad[tad$chr==max(as.numeric(dat$chr),na.rm=TRUE),]
    in_tad<-tad[between(snp_pos,tad$start_position, tad$end_position)]
    return(in_tad)
  })

  output$tadBoundaries<-renderText({
    in_tad<-in_tad()
    if(nrow(in_tad)<1){return(paste0("In a TAD!"))}
    else({return(paste0("Not in a TAD"))})
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
    if (length(x)>1){
    a("Take me to ClinVar", href=paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=",max(as.numeric(y$chr),na.rm = TRUE),"%5Bchr%5D+AND+",min(as.numeric(y$pos_hg38),na.rm=TRUE),"%3A",max(as.numeric(y$pos_hg38),na.rm=TRUE),"%5Bchrpos37%5D"), target="_blank")
    } else if (length(x)==1){
      a("Take me to ClinVar", href=paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=",max(as.numeric(y$chr),na.rm = TRUE),"%5Bchr%5D+AND+",min((as.numeric(y$pos_hg38)),na.rm=TRUE)-53500,"%3A",max((as.numeric(y$pos_hg38)),na.rm=TRUE)+53500,"%5Bchrpos37%5D"), target="_blank")
    }
  })
  
  output$ucsc1<-renderUI({
    x<-snps()
    y<-dat()
    if (length(x)>1){
    a("Take me to Genome Browser", href=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr",max(as.numeric(y$chr),na.rm = TRUE),"%3A",min(as.numeric(y$pos_hg38),na.rm=TRUE),"%2D",max(as.numeric(y$pos_hg38),na.rm=TRUE),"&hgsid=596717155_di6qMTAMSs8fhJcRiuqsjlcsIxKA"), target="_blank")
    }else if (length(x)==1){
      a("Take me to Genome Browser", href=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr",max(as.numeric(y$chr),na.rm = TRUE),"%3A",min((as.numeric(y$pos_hg38)),na.rm=TRUE)-53500,"%2D",max((as.numeric(y$pos_hg38)),na.rm=TRUE)+53500,"&hgsid=596717155_di6qMTAMSs8fhJcRiuqsjlcsIxKA"), target="_blank")
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

