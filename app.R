library(shinydashboard)
library(shiny)
library(haploR)
library(data.table)
library(biomaRt)
library(ggplot2)
library(shinycssloaders)
#library(ggrepel)
library(plotly)
library(jsonlite)
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
                                  tabBox(title="Select Output",
                                      tabPanel("Source",
                                               selectInput("pop","Population",c("EUR","AFR","AMR","ASN"), selected="EUR"),
                                               sliderInput("value","LD threshold",min=0,max=1,value=0.8)),
                                      tabPanel("HaploReg",
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
                                                                                          "NCBI Reference Sequence distance"="RefSeq_distance","Annotated proteins"="dbSNP_functional_annotation"),inline = TRUE)),
                                      tabPanel("RegulomeDB",
                                               checkboxGroupInput("parameters2","Regulome",c("Chromosome"="#chromosome",
                                                                                             "Coordinates"="coordinate",
                                                                                             "Hits"="hits", "Score"="score_anno"))),
                                      tabPanel("eQTL",
                                               uiOutput("eTissues")),
                                      tabPanel("Oncotator",
                                               checkboxGroupInput("oncoParameters1","General",c("ACHILLES_Lineage_Results_Top_Genes","CCLE_By_Gene_total_mutations_in_gene",  
                                                                                                 "COSMIC_FusionGenes_fusion_genes",           
                                                                                                 "COSMIC_Tissue_tissue_types_affected","COSMIC_Tissue_total_alterations_in_gene",   
                                                                                                 "Familial_Cancer_Genes_Reference","Familial_Cancer_Genes_Syndrome",            
                                                                                                 "Familial_Cancer_Genes_Synonym","HGNC_Accession.Numbers",                   
                                                                                                 "HumanDNARepairGenes_Role","MutSig.Published.Results_Published_Results",
                                                                                                 "TCGAScape_Amplification_Peaks","TCGAScape_Deletion_Peaks",                  
                                                                                                 "TUMORScape_Amplification_Peaks","TUMORScape_Deletion_Peaks",                 
                                                                                                 "alt_allele",                                
                                                                                                 "build","chr",                                      
                                                                                                 "class","end",                                       
                                                                                                 "gene","protein_change",                            
                                                                                                 "ref_allele","start",                                     
                                                                                                 "strand","transcripts" ), inline=TRUE),
                                               checkboxGroupInput("oncoParameters2","Cancer Gene Census",c("CGC_Cancer.Germline.Mut","CGC_Cancer.Molecular.Genetics",             
                                                                                                           "CGC_Cancer.Somatic.Mut","CGC_Cancer.Syndrome",                       
                                                                                                           "CGC_Chr","CGC_Chr.Band",                              
                                                                                                           "CGC_GeneID","CGC_Mutation.Type",                         
                                                                                                           "CGC_Name","CGC_Other.Germline.Mut",                    
                                                                                                           "CGC_Other.Syndrome.Disease","CGC_Tissue.Type",                           
                                                                                                           "CGC_Translocation.Partner","CGC_Tumour.Types...Somatic.Mutations.",     
                                                                                                           "CGC_Tumour.Types..Germline.Mutations."),inline = TRUE),
                                               checkboxGroupInput("oncoParameters3","HUGO Gene Nomenclature Committee",c("HGNC_Approved.Name","HGNC_CCDS.IDs",                             
                                                                                             "HGNC_Chromosome","HGNC_Date.Modified",                        
                                                                                             "HGNC_Date.Name.Changed","HGNC_Date.Symbol.Changed",                  
                                                                                             "HGNC_Ensembl.Gene.ID","HGNC_Ensembl.ID.supplied.by.Ensembl.",      
                                                                                             "HGNC_Entrez.Gene.ID","HGNC_Entrez.Gene.ID.supplied.by.NCBI.",     
                                                                                             "HGNC_Enzyme.IDs","HGNC_Gene.family.description",              
                                                                                             "HGNC_HGNC.ID","HGNC_Locus.Group",                          
                                                                                             "HGNC_Locus.Type","HGNC_Name.Synonyms",                        
                                                                                             "HGNC_OMIM.ID.supplied.by.NCBI.","HGNC_Previous.Names",                      
                                                                                             "HGNC_Previous.Symbols","HGNC_Primary.IDs",                          
                                                                                             "HGNC_Pubmed.IDs","HGNC_Record.Type",                          
                                                                                             "HGNC_RefSeq.IDs","HGNC_RefSeq.supplied.by.NCBI.",            
                                                                                             "HGNC_Secondary.IDs","HGNC_Status",                               
                                                                                             "HGNC_Synonyms","HGNC_UCSC.ID.supplied.by.UCSC." ,           
                                                                                             "HGNC_UniProt.ID.supplied.by.UniProt.","HGNC_VEGA.IDs"), inline=TRUE),
                                               checkboxGroupInput("oncoParameters4","UniProt",c("UniProt_AA_experimental_info","UniProt_AA_natural_variation",              
                                                                                                "UniProt_AA_region","UniProt_AA_site",                           
                                                                                                "UniProt_DrugBank","UniProt_GO_Biological_Process",             
                                                                                                "UniProt_GO_Cellular_Component","UniProt_GO_Molecular_Function",            
                                                                                                "UniProt_alt_uniprot_accessions","UniProt_uniprot_accession",                 
                                                                                                "UniProt_uniprot_entry_name"), inline=TRUE))
                                  )
                                ),
                                fluidRow(
                                  tabBox(title = "Variant Annotation",
                                                tabPanel("HaploReg",
                                                         tableOutput("LDtable1")),
                                                tabPanel("RegulomeDB",
                                                         tableOutput("LDtable2")),
                                                tabPanel("TADs",
                                                         actionButton("tadButton","Look up TADs"),
                                                         textOutput("tadBoundaries"),
                                                         uiOutput("hic1"))
                                                ),
                                  tabBox(title="Gene Annotation",
                                         tabPanel("ENSEMBL",
                                                  tableOutput("geneTable")),
                                         tabPanel("Oncotator",
                                                  tableOutput("oncoTable")),
                                         tabPanel("eQTL",
                                                  tableOutput("eTable1"),
                                                  uiOutput("eqtl1"))),
                                  column(12,align="center",offset=3,
                                         tabBox(title="Other",
                                         tabPanel("Links",
                                                  uiOutput("clinical1"),
                                                  uiOutput("ucsc1")),
                                         tabPanel("Visual",
                                                  h5(helpText("TAD boundaries need to be calculated before being plotted!")),
                                                  withSpinner(plotlyOutput("plot1",height="450px"),color = "#00ffff", type = 6))))
                                )),
                        tabItem(tabName = "tab2",
                                fluidRow(
                                  box(title="App Details",
                                      h5(helpText("LD is calculated from 1000 Genomes Phase 1 (http://www.internationalgenome.org), and queried from HaploReg (http://archive.broadinstitute.org/mammals/haploreg/haploreg.php)
                                                  To perform similar queries in R please check out the haploR package!
                                                  For TAD visualization check out the Yue Lab (http://promoter.bx.psu.edu/hi-c/).
                                                  TAD locations are based off of those defined by Dixon et al in 'Topological domains in mammalian genomes identified by analysis of chromatin interactions'."
                                      ))),
                                  box(title="Development Team",
                                      h5(helpText("Programming: Jordan Creed, Travis Gerke")),
                                      h5(helpText("Scientific Input: Alvaro Monteiro")),
                                      h5(helpText("Website: http://travisgerke.com"))),
                                  box(title="Other resources",
                                      h5(helpText("Aiden Lab")),
                                      a("Juicebox", href="http://www.aidenlab.org/juicebox/", target="_blank"))
                                      ),
                                fluidRow(
                                  box(title="Notes",
                                      h5(helpText("If no snps are in LD above the specified threshold than a range of 53500 bp is applied to either side of the snp.
                                                  This range is then applied to querying data from Oncotator, ENSEMBL, ClinVar and the Genome Browser.
                                                  If snps in LD exist than the range is based on the minimum and maximum value of all snps in LD.")))
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
    if(nrow(etest3)>1){
    return(checkboxGroupInput("tissue","Tissues",choices=opt, inline = TRUE))
    } else(print("No Tissues to Show"))
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
    etest3<-matrix(etest2,nrow=length(etest), ncol=4,byrow=TRUE)
    etest3<-as.data.frame(etest3)
    colnames(etest3)<-c("Source","Tissue","Gene","p")
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
      a("Take me to Genome Browser", href=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr",max(as.numeric(y$chr),na.rm = TRUE),"%3A",min(as.numeric(y$pos_hg38),na.rm=TRUE),"%2D",max(as.numeric(y$pos_hg38),na.rm=TRUE),"&hgsid=598506407_cis2LZUJLabCsy1N2YPEuJv8vbBZ"), target="_blank")
    }else if (length(x)==1){
      a("Take me to Genome Browser", href=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr",max(as.numeric(y$chr),na.rm = TRUE),"%3A",min((as.numeric(y$pos_hg38)),na.rm=TRUE)-53500,"%2D",max((as.numeric(y$pos_hg38)),na.rm=TRUE)+53500,"&hgsid=598506407_cis2LZUJLabCsy1N2YPEuJv8vbBZ"), target="_blank")
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
  
  output$geneTable<-renderTable({
    ld<-dat()
    chr<-min(ld$chr,na.rm=TRUE)
    if(nrow(ld)>1){
      min<-min(ld$pos_hg38,na.rm=TRUE)
      max<-max(ld$pos_hg38,na.rm=TRUE)
    } else{
      min<-min(ld$pos_hg38,na.rm=TRUE)-53500
      max<-max(ld$pos_hg38,na.rm=TRUE)+53500
    }
    
    genes<-getBM(attributes = c("hgnc_symbol","start_position","end_position"),
                 filters=c("chromosomal_region"), values=paste0(chr,":",min,":",max),mart = ensembl54)
    return(genes)
  })
  
  output$oncoTable<-renderTable({
    ld<-dat()
    chr<-min(ld$chr,na.rm=TRUE)
    if(nrow(ld)>1){
    min<-min(ld$pos_hg38,na.rm=TRUE)
    max<-max(ld$pos_hg38,na.rm=TRUE)
    } else{
      min<-min(ld$pos_hg38,na.rm=TRUE)-53500
      max<-max(ld$pos_hg38,na.rm=TRUE)+53500
    }
    
    x<-fromJSON(paste0("http://portals.broadinstitute.org/oncotator/genes/",chr,"_",min,"_",max,"/"))
    
    genes<-as.data.frame(x[[1]])
    
    for (i in 2:length(x)){
      gene_dat<-as.data.frame(x[[i]])
      genes<-rbind(genes, gene_dat)
    }
    
    genes<-genes[,c("gene",input$oncoParameters1,input$oncoParameters2,input$oncoParameters3,input$oncoParameters4)]
    return(genes)
    
  })
  
  output$plot1<-renderPlotly({
    # create plotting without having to search 
    ld<-dat()
    query_snps<-ld[ld$is_query_snp==1,]
    ld_snps<-ld[ld$is_query_snp==0,]

    x_min<-min(ld$pos_hg38,na.rm=TRUE)
    x_max<-max(ld$pos_hg38,na.rm=TRUE)
    
    ldBlocks<-ggplot(ld)+
      geom_segment(data=ld[ld$is_query_snp==1,],aes(x=as.vector(tapply(ld$pos_hg38, ld$query_snp_rsid, min)),y=1,xend=as.vector(tapply(ld$pos_hg38, ld$query_snp_rsid, max)),yend=1,color=as.factor(query_snp_rsid), size=30, text=paste0("SNP: ",query_snp_rsid)), alpha=0.5)+
      geom_vline(data = query_snps, aes(xintercept=pos_hg38, text=paste0("SNP: ",query_snp_rsid))) +
      geom_vline(data = ld_snps, aes(xintercept=pos_hg38, alpha=0.1, text=paste0("SNP: ",rsID)), color="grey")+ # add different color from query snps to make more visible?
      theme(legend.position = "none",axis.text.y=element_blank(),axis.title.y=element_blank(),panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) + 
      xlab("BP") + coord_cartesian(ylim=c(0.75,3.5))+
      scale_y_continuous(breaks = c(1,2,3)) 
    
    if(input$tadButton!=0){
      tad<-in_tad()
      
      if(nrow(tad)>=1){
        x_min<-min(c(min(ld$pos_hg38,na.rm=TRUE),min(tad$start_position)))
      } else {x_min<-min(ld$pos_hg38,na.rm=TRUE)}
      
      if(nrow(tad)>=1){
        x_max<-max(c(max(ld$pos_hg38,na.rm=TRUE),max(tad$end_position)))
      } else {x_max<-max(ld$pos_hg38,na.rm=TRUE)}
      
      genes<-getBM(attributes = c("hgnc_symbol","start_position","end_position"),
                   filters=c("chromosomal_region"), values=paste0(max(ld$chr,na.rm = TRUE),":",x_min,":",x_max),mart = ensembl54)
      colnames(genes)<-c("Symbol","Start","End")
      
      ldBlocks<-ldBlocks + 
        geom_segment(data=tad, aes(x=tad$start_position, y=2, xend=tad$end_position, yend=2, size=30, text=paste0("Start: ",start_position,"\n","End",end_position)))+
        geom_segment(data=genes,aes(x=Start,y=3,xend=End,yend=3,color=Symbol,size=30, text=paste0("Symbol: ",Symbol,"\n", "Start: ",Start,"\n", "End",End)),alpha=0.5)+
        annotate("text",x=x_min, y=2.25, label="TADs", color="purple", angle=90) +
        annotate("text",x=x_min, y=1.25, label="LD", color="purple", angle=90)+
        annotate("text",x=x_min, y=3.25, label="Genes", color="purple", angle=90)
      
    } else {
      x_min<-min(ld$pos_hg38,na.rm=TRUE)
      x_max<-max(ld$pos_hg38,na.rm=TRUE)
       
      genes<-getBM(attributes = c("hgnc_symbol","start_position","end_position"),
                  filters=c("chromosomal_region"), values=paste0(max(ld$chr,na.rm = TRUE),":",x_min,":",x_max),mart = ensembl54)
      colnames(genes)<-c("Symbol","Start","End")
      ldBlocks<-ldBlocks + geom_segment(data=genes,aes(x=Start,y=3,xend=End,yend=3,color=Symbol,size=30, text=paste0("Symbol: ",Symbol,"\n", "Start: ",Start,"\n", "End",End)),alpha=0.5)+
        annotate("text",x=x_min, y=1.25, label="LD", color="purple", angle=90)+
        annotate("text",x=x_min, y=3.25, label="Genes", color="purple", angle=90)
      }

    #return(ldBlocks)
    ggplotly(ldBlocks, tooltip="text") %>% 
      layout(autosize=TRUE)
  })
  
}

################################################################################################### 
shinyApp(ui = ui, server = server)

