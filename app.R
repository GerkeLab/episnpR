library(shinydashboard)
library(shiny)
library(haploR)
library(data.table)
library(biomaRt)
library(ggplot2)
library(shinycssloaders)
#library(ggrepel)
library(plotly)
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
                                      br(),
                                      checkboxGroupInput("parameters","HaploR",c("Chromosome"="chr","Position"="pos_hg38",                         
                                                                                            "D'"="D'","Query SNP"="is_query_snp",                       
                                                                                            "Reference allele"="ref","Alternative allele"="alt","LD(AFR)"="AFR",                        
                                                                                            "LD(AMR)"="AMR","LD(ASN)"="ASN","LD(EUR)"="EUR",                        
                                                                                            "GERP scores"="GERP_cons","SiPhy scores"="SiPhy_cons","Chromatin States" ="Chromatin_States",           
                                                                                            "Imputed Chromatin States"="Chromatin_States_Imputed","Chromatin Marks"="Chromatin_Marks","DNAse"="DNAse",                      
                                                                                            "Proteins","eQTL","GWAS study name"="gwas",                       
                                                                                            "GRASP study name"="grasp","Motifs","GENCODE transcript ID"="GENCODE_id",                 
                                                                                            "GENCODE gene name"="GENCODE_name","GENCODE direction"="GENCODE_direction","GENCODE distance"="GENCODE_distance",           
                                                                                            "NCBI Reference Sequence Accession number"="RefSeq_id","NCBI Reference Sequence name"="RefSeq_name","NCBI Reference Sequence direction"="RefSeq_direction",           
                                                                                            "NCBI Reference Sequence distance"="RefSeq_distance","Annotated proteins"="dbSNP_functional_annotation"),inline = TRUE),
                                      br(),
                                      checkboxGroupInput("parameters2","Regulome",c("Chromosome"="#chromosome",
                                                                                           "Coordinates"="coordinates",
                                                                                           "Hits"="hits", "Score"="score_anno")),
                                      h5(helpText("Tissues")),
                                      uiOutput("eTissues")
                                  )
                                ),
                                fluidRow(
                                  column(12,align="center",offset=3,
                                         tabBox(title = "Output",
                                                tabPanel("HaploReg",
                                                         tableOutput("LDtable1")),
                                                tabPanel("RegulomeDB",
                                                         tableOutput("LDtable2")),
                                                tabPanel("eQTL",
                                                         tableOutput("eTable1"),
                                                         uiOutput("eqtl1")), 
                                                tabPanel("TADs",
                                                         actionButton("tadButton","Look up TADs"),
                                                         textOutput("tadBoundaries"),
                                                         uiOutput("hic1")),
                                                tabPanel("Other",
                                                         uiOutput("clinical1"),
                                                         uiOutput("ucsc1")),
                                                tabPanel("Visuals",
                                                         withSpinner(plotlyOutput("plot1",height="450px"),color = "#00ffff", type = 6)))
                                ))),
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
  
  ensembl54=useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")
  dbsnp = useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
  
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
  
  output$eTissues<-renderUI({
    dat<-dat()
    etest<-unlist(strsplit(dat$eQTL,";"))
    etest2<-unlist(strsplit(etest,","))
    etest3<-matrix(etest2,nrow=length(etest), ncol=4, byrow = TRUE)
    etest3<-as.data.frame(etest3)
    etest3<-etest3[!duplicated(etest3$V2),]
    opt<-etest3$V2
    return(checkboxGroupInput("tissue","Tissues",choices=opt, inline = TRUE))
  })
  
  dat2<-eventReactive(input$update1,{
    if(input$snpList==""){
      dat<-sample()
      snps<-dat[,1]
      x<-queryRegulome(query = snps)
      x<-as.data.frame(x$res.table)
      x$score<-as.character(x$score)
      x$score_anno<-NA
      for (i in nrow(x)){
        if (x$score[i]=="1a"){x$score_anno[i]<-"eQTL + TF binding + matched TF motif + matched DNase Footprint + DNase peak"}
        else if (x$score[i]=="1b"){x$score_anno[i]<-"eQTL + TF binding + any motif + DNase Footprint + DNase peak"}
        else if (x$score[i]=="1c"){x$score_anno[i]<-"eQTL + TF binding + matched TF motif + DNase peak"}
        else if (x$score[i]=="1d"){x$score_anno[i]<-"eQTL + TF binding + any motif + DNase peak"}
        else if (x$score[i]=="1e"){x$score_anno[i]<-"eQTL + TF binding + matched TF motif"}
        else if (x$score[i]=="1f"){x$score_anno[i]<-"eQTL + TF binding / DNase peak"}
        else if (x$score[i]=="2a"){x$score_anno[i]<-"TF binding + matched TF motif + matched DNase Footprint + DNase peak"}
        else if (x$score[i]=="2b"){x$score_anno[i]<-"TF binding + any motif + DNase Footprint + DNase peak"}
        else if (x$score[i]=="2c"){x$score_anno[i]<-"TF binding + matched TF motif + DNase peak"}
        else if (x$score[i]=="3a"){x$score_anno[i]<-"TF binding + any motif + DNase peak"}
        else if (x$score[i]=="3b"){x$score_anno[i]<-"TF binding + matched TF motif"}
        else if (x$score[i]=="4"){x$score_anno[i]<-"TF binding + DNase peak"}
        else if (x$score[i]=="5"){x$score_anno[i]<-"TF binding or DNase peak"}
        else {x$score_anno[i]<-"Other"}
      }
      #x<-x[,input$parameters]
      return(x)
    }
    if(input$snpList!=""){
      snps<-as.character(unlist(strsplit(input$snpList,",")))
      x<-queryRegulome(query = snps)
      x<-as.data.frame(x$res.table)
      x$score<-as.character(x$score)
      x$score_anno<-NA
      for (i in 1:nrow(x)){
        if (x$score[i]=="1a"){x$score_anno[i]<-"eQTL + TF binding + matched TF motif + matched DNase Footprint + DNase peak"}
        else if (x$score[i]=="1b"){x$score_anno[i]<-"eQTL + TF binding + any motif + DNase Footprint + DNase peak"}
        else if (x$score[i]=="1c"){x$score_anno[i]<-"eQTL + TF binding + matched TF motif + DNase peak"}
        else if (x$score[i]=="1d"){x$score_anno[i]<-"eQTL + TF binding + any motif + DNase peak"}
        else if (x$score[i]=="1e"){x$score_anno[i]<-"eQTL + TF binding + matched TF motif"}
        else if (x$score[i]=="1f"){x$score_anno[i]<-"eQTL + TF binding / DNase peak"}
        else if (x$score[i]=="2a"){x$score_anno[i]<-"TF binding + matched TF motif + matched DNase Footprint + DNase peak"}
        else if (x$score[i]=="2b"){x$score_anno[i]<-"TF binding + any motif + DNase Footprint + DNase peak"}
        else if (x$score[i]=="2c"){x$score_anno[i]<-"TF binding + matched TF motif + DNase peak"}
        else if (x$score[i]=="3a"){x$score_anno[i]<-"TF binding + any motif + DNase peak"}
        else if (x$score[i]=="3b"){x$score_anno[i]<-"TF binding + matched TF motif"}
        else if (x$score[i]=="4"){x$score_anno[i]<-"TF binding + DNase peak"}
        else if (x$score[i]=="5"){x$score_anno[i]<-"TF binding or DNase peak"}
        else {x$score_anno[i]<-"Other"}
      }
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
    tad<-fread("http://compbio.med.harvard.edu/modencode/webpage/hic/IMR90_domains_hg19.bed")
    colnames(tad)<-c("chr","start_position","end_position")
    tad$chr<-gsub("chr","",tad$chr)
    tad$chr<-as.numeric(tad$chr)
    tad<-tad[!is.na(tad$chr),]
    snps<-snps()
    dat<-dat()
    dat<-dat[dat$rsID %in% snps,]
    snp_pos<-as.numeric(dat$pos_hg38)
    tad<-tad[tad$chr==max(as.numeric(dat$chr),na.rm=TRUE),]
    in_tad<-tad[between(snp_pos,tad$start_position, tad$end_position)]
    # if(nrow(in_tad>=1)){
    #   tad_region<-paste0(in_tad$chr,":",in_tad$start_position,":",in_tad$end_position)
    #   genes<-getBM(attributes = c("hgnc_symbol","start_position","end_position"),
    #                filters=c("chromosomal_region"), values=tad_region,mart = ensembl54)
    # }
    return(in_tad)
  })
  
  
  output$tadBoundaries<-renderText({
    in_tad<-in_tad()
    if(nrow(in_tad)<1){return(paste0("Not in a TAD!"))}
    else({return(paste0("In a TAD! The TAD ranges from ",in_tad$start_position," to ",in_tad$end_position))})
  })
  
  output$eTable1<-renderTable({
    dat<-dat()
    etest<-unlist(strsplit(dat$eQTL,";"))
    etest2<-unlist(strsplit(etest,","))
    etest3<-matrix(etest2,nrow=length(etest), ncol=4)
    etest3<-as.data.frame(etest3)
    colnames(etest3)<-c("Gene","Tissue","Source","p")
    return(etest3[etest3$Tissue %in% input$tissue,])
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
  
  output$eqtl1<-renderUI({
    x<-snps()
    y<-dat()
    if (length(x)>1){
      a("Take me to GTEx", href=paste0("https://www.gtexportal.org/home/browseEqtls?location=chr",max(as.numeric(y$chr),na.rm = TRUE),":",min((as.numeric(y$pos_hg38)),na.rm=TRUE),"-",max((as.numeric(y$pos_hg38)),na.rm=TRUE)), target="_blank")
    } else if (length(x)==1){
      a("Take me to GTEx", href=paste0("https://www.gtexportal.org/home/snp/",x), target="_blank")
    }
  })
  
  
  
  output$LDtable1<-renderTable({
    x<-dat()
    #x<-x[,input$parameters]
    x[,c("query_snp_rsid","rsID","r2",input$parameters)]
    #return(x[,input$parameters])
  })
  
  output$LDtable2<-renderTable({
    x<-dat2()
    x[,c("rsid",input$parameters2)]
  })
  
  output$plot1<-renderPlotly({
    # create plotting without having to search 
    ld<-dat()
    query_snps<-ld[ld$is_query_snp==1,]
    ld_snps<-ld[ld$is_query_snp==0,]
    
    tad<-in_tad()
    
    if(nrow(tad)>=1){
    x_min<-min(c(min(ld$pos_hg38,na.rm=TRUE),min(tad$start_position)))
    } else {x_min<-min(ld$pos_hg38,na.rm=TRUE)}
    
    if(nrow(tad)>=1){
    x_max<-max(c(max(ld$pos_hg38,na.rm=TRUE),max(tad$end_position)))
    } else {x_max<-max(ld$pos_hg38,na.rm=TRUE)}
    
    genes<-getBM(attributes = c("hgnc_symbol","start_position","end_position"),
                 filters=c("chromosomal_region"), values=paste0(max(ld$chr,na.rm = TRUE),":",x_min,":",x_max),mart = ensembl54)
    
    ldBlocks<-ggplot(ld)+
      geom_segment(data=ld,aes(x=min(ld$pos_hg38,na.rm = TRUE),y=1,xend=max(ld$pos_hg38,na.rm = TRUE),yend=1, color=as.factor(ld$query_snp_rsid), size=30))+ 
      geom_vline(data = query_snps, aes(xintercept=pos_hg38)) +
      geom_vline(data = ld_snps, aes(xintercept=pos_hg38, alpha=0.1), color="grey")+ # add different color from query snps to make more visible?
      geom_segment(data=genes,aes(x=start_position,y=3,xend=end_position,yend=3,color=hgnc_symbol,size=30),alpha=0.5)+
      annotate("text",x=x_min, y=1.25, label="LD", color="purple", angle=90)+
      annotate("text",x=x_min, y=3.25, label="Genes", color="purple", angle=90)+
      annotate("text",x=x_min, y=2.25, label="TADs", color="purple", angle=90)+
      theme(legend.position = "none",axis.text.y=element_blank(),axis.title.y=element_blank(),panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) + 
      xlab("BP") + coord_cartesian(ylim=c(0.75,3.5))+
      scale_y_continuous(breaks = c(1,2,3)) 
    
    if(nrow(tad)>=1){
      ldBlocks<-ldBlocks + geom_segment(data=tad, aes(x=tad$start_position, y=2, xend=tad$end_position, yend=2, size=30))
    } else {ldBlocks<-ldBlocks}

    #return(ldBlocks)
    ggplotly(ldBlocks) %>% 
      layout(autosize=TRUE)
  })
  
}

################################################################################################### 
shinyApp(ui = ui, server = server)

