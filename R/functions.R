
#'GO enrichment analysis
#'
#'This function does GO enrichment on a user defined GO universe
#'
#' @usage GOanalysis(genes, universe = pt210_GO_universe,
#'                organism = "Populus trichocarpa", pv = 0.01,
#'                ontology=c("BP","MF","CC"), testDirection = "over",
#'                conditional=FALSE)
#' @param genes Provide a vector of gene names to be tested
#' @param universe Provide a user define GO universe
#' @param organism Specify the name of the organism
#' @param pv Specify a p-value threshhold for significance (default = 0.01)
#' @param ontology Specify the ontologies to test (Default = c("BP","MF","CC"))
#' @param testDirection A string which can be either "over" or "under". This determines whether the test performed detects over or under represented GO terms. See \link[Category]{GSEAGOHyperGParams} (Defalut = "over")
#' @param conditional A logical indicating whether the calculation should condition on the GO structure. (GO only)
#' @import annotate
#' @import GO.db
#' @import GSEABase
#' @import GOstats
#' @importMethodsFrom AnnotationDbi GOFrame
#' @importMethodsFrom AnnotationDbi GOAllFrame
#' @importMethodsFrom GSEABase GeneSetCollection
#' @importMethodsFrom GSEABase GOCollection
#' @importMethodsFrom Category GSEAGOHyperGParams
#' @importMethodsFrom GOstats hyperGTest
#' @return Returns a list object that contains the ontology results for Biological Processes (BP), Molecular Functions (MF) and Cellular Components (CC)
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
GOanalysis<-function (genes, universe = pt210_GO_Universe, organism = "Populus trichocarpa", pv=0.01, ontology=c("BP","MF","CC"), testDirection = "over", conditional=FALSE)
{

  goFrame<-GOFrame(universe,organism=organism)
  goAllFrame<-GOAllFrame(goFrame)


  gsc <- GSEABase::GeneSetCollection(goAllFrame, setType = GSEABase::GOCollection())
  uni<-as.character(unique(universe$frame.gene_id))
  universe$frame.gene_id <- as.character(universe$frame.gene_id)
  g_id<-as.character(intersect(universe$frame.gene_id,genes))

  output<-list()
  if("MF" %in% ontology)
  {
  params_MF <- GSEAGOHyperGParams(name="Annotation Params MF",
                                  geneSetCollection=gsc,
                                  geneIds =g_id,
                                  universeGeneIds = uni,
                                  ontology = "MF",
                                  pvalueCutoff = pv,
                                  conditional = conditional,
                                  testDirection = testDirection)

  Over_MF <- hyperGTest(params_MF)
  output$MF = Over_MF
  }

  if("BP" %in% ontology)
  {
  params_BP <- GSEAGOHyperGParams(name="Annotation Params BP",
                                  geneSetCollection=gsc,
                                  geneIds =g_id,
                                  universeGeneIds = uni,
                                  ontology = "BP",
                                  pvalueCutoff = pv,
                                  conditional = conditional,
                                  testDirection = testDirection)

  Over_BP <- hyperGTest(params_BP)
  output$BP = Over_BP
  }

  if("CC" %in% ontology)
  {
  params_CC <- GSEAGOHyperGParams(name="Annotation Params CC",
                                  geneSetCollection=gsc,
                                  geneIds =g_id,
                                  universeGeneIds = uni,
                                  ontology = "CC",
                                  pvalueCutoff = pv,
                                  conditional = conditional,
                                  testDirection = testDirection)

  Over_CC <- hyperGTest(params_CC)
  output$CC = Over_CC
  }

  return(output)
}



#'GO enrichment analysis using Arabidopsis TAIR10
#'
#'This function does GO enrichment on TAIR10 GO universe
#'
#' @usage atGOanalysis(genes, universe = TAIR10_GO_universe,
#'                organism = "Arabidopsis thaliana", pv = 0.01,
#'                ontology=c("BP","MF","CC"), testDirection = "over")
#' @param genes Provide a vector of gene names to be tested
#' @param universe Arabidopsis TAIR10 GO universe
#' @param organism Specify the name of the organism
#' @param pv Specify a p-value threshhold for significance (default = 0.01)
#' @param ontology Specify the ontologies to test (Default = c("BP","MF","CC"))
#' @param testDirection A string which can be either "over" or "under". This determines whether the test performed detects over or under represented GO terms. See \link[Category]{GSEAGOHyperGParams} (Defalut = "over")
#' @import annotate
#' @import GO.db
#' @import GSEABase
#' @import GOstats
#' @importMethodsFrom AnnotationDbi GOFrame
#' @importMethodsFrom AnnotationDbi GOAllFrame
#' @importMethodsFrom GSEABase GeneSetCollection
#' @importMethodsFrom GSEABase GOCollection
#' @importMethodsFrom Category GSEAGOHyperGParams
#' @importMethodsFrom GOstats hyperGTest
#' @return Returns a list object that contains the ontology results for Biological Processes (BP), Molecular Functions (MF) and Cellular Components (CC)
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
atGOanalysis<-function (genes, universe = TAIR10_GO_Universe, organism = "Arabidopsis thaliana", pv=0.01, ontology=c("BP","MF","CC"), testDirection = "over")
{

  #remove transcript numbers from AT names
  genes1<-sub("(AT\\d+G\\d+)\\.\\d+","\\1",genes,perl=TRUE)


  goFrame<-GOFrame(universe,organism=organism)
  goAllFrame<-GOAllFrame(goFrame)

  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

  uni<-as.character(unique(universe$frame.gene_id))
  universe$frame.gene_id <- as.character(universe$frame.gene_id)
  g_id<-as.character(intersect(universe$frame.gene_id,genes1))

  output<-list()
  if("MF" %in% ontology)
  {
    params_MF <- GSEAGOHyperGParams(name="Annotation Params MF",
                                    geneSetCollection=gsc,
                                    geneIds =g_id,
                                    universeGeneIds = uni,
                                    ontology = "MF",
                                    pvalueCutoff = pv,
                                    conditional = FALSE,
                                    testDirection = testDirection)

    Over_MF <- hyperGTest(params_MF)
    output$MF = Over_MF
  }

  if("BP" %in% ontology)
  {
    params_BP <- GSEAGOHyperGParams(name="Annotation Params BP",
                                    geneSetCollection=gsc,
                                    geneIds =g_id,
                                    universeGeneIds = uni,
                                    ontology = "BP",
                                    pvalueCutoff = pv,
                                    conditional = FALSE,
                                    testDirection = testDirection)

    Over_BP <- hyperGTest(params_BP)
    output$BP = Over_BP
  }

  if("CC" %in% ontology)
  {
    params_CC <- GSEAGOHyperGParams(name="Annotation Params CC",
                                    geneSetCollection=gsc,
                                    geneIds =g_id,
                                    universeGeneIds = uni,
                                    ontology = "CC",
                                    pvalueCutoff = pv,
                                    conditional = FALSE,
                                    testDirection = testDirection)

    Over_CC <- hyperGTest(params_CC)
    output$CC = Over_CC
  }

  return(output)

}

#'Trace GO terms to gene names
#'
#'This function maps GO terms back to gene names
#'
#' @usage atGO2Potri(go_id, mapping = MAPPING_tair10_GO_to_pt210)
#' @param go_id GO term in format GO:0000000
#' @param mapping Specify a object to map GO to genes
#' @return Returns a vector of genes that map to a GO term
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
atGO2Potri<-function(go_id, mapping = MAPPING_tair10_GO_to_pt210)
{

  if(length(go_id)==0)
  {
    return(print("NO arabidopsis GO term was recieved"))
  } else if (length(go_id)>1) {
    return(print("More than one GO term was recieved"))
  } else if (length(grep("GO\\:\\d+",go_id,perl=TRUE))!=1) {
    return(print("GO term format is incorrect and shoulf be GO:0000000"))
  } else {
    return(as.character(unique(mapping$V1[which(mapping$frame.go_id==go_id)])))
  }

}

#'Trace GO terms to gene names and annotations
#'
#'This function maps GO terms back to gene names and includes annotations
#'
#' @usage atGO2Potri(go_id, mapping = MAPPING_tair10_GO_to_pt210
#' annotation = Ptrichocarpa_210_annotation_primary, col_index=2)
#' @param go_id GO term in format GO:0000000
#' @param mapping Specify a object to map GO to genes
#' @param annotation Specify a object to annotate genes
#' @param col_index Specify which column contains the gene names
#' @return Returns a data frame of annotated genes
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
atGO2PotriAnno<-function(go_id, mapping = MAPPING_tair10_GO_to_pt210, annotation = Ptrichocarpa_210_annotation_primary, col_index=2)
{

  if(length(go_id)==0)
  {
    return(print("NO arabidopsis GO term was recieved"))
  } else if (length(grep("GO\\:\\d+",go_id,perl=TRUE))!=length(go_id)) {
    return(print("GO term format is incorrect and shoulf be GO:0000000"))
  } else if (length(go_id)>1) {
    potri<-as.character(unique(mapping$V1[which(mapping$frame.go_id%in%go_id)]))
    return(annotation[which(annotation[,col_index] %in% potri),])
  } else {
    potri<-as.character(unique(mapping$V1[which(mapping$frame.go_id==go_id)]))
    return(annotation[which(annotation[,col_index] %in% potri),])
  }

}


#'Remove GO term from GO_object
#'
#'This function removes a GO term from GO_object
#'
#' @usage atGO2Potri(go_id, mapping = MAPPING_tair10_GO_to_pt210)
#' @param go_id GO term in format GO:0000000
#' @param GO_object A Go object from GOanalysis
#' @return Returns update GO_object
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
removeTERM<-function(GOid, GO_object)
{
  tmp<-GO_object; #GOresults[[1]]$BP

  for(g in GOid)
  {
    if(g %in% tmp@goDag@nodes)
    {
      index<-which(tmp@goDag@nodes==g)
      #remove term

      tmp@pvalue.order<-tmp@pvalue.order[which(tmp@pvalue.order!=index)]
      #update numbers
      tmp@pvalue.order[which(tmp@pvalue.order>index)]<-as.integer(tmp@pvalue.order[tmp@pvalue.order>index]-1)

      tmp@goDag<-removeNode(g,tmp@goDag)
    }
  }
  return(tmp)
}


#'No in
#'
#'This function does the opposite of in
#'
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
`%ni%`=Negate(`%in%`)

#
#'Parse GO results using grep and plot
#'
#'This function creates a summary plot of GO$BP results by parsing GO results using grep and a list of search terms
#'
#' @usage ParseGOBPnPlot(grepList,GOresults, minT=0, Rorder=NULL,
#'              my_palette=colorRampPalette(c("white", "red"))(n = 8),
#'              main = "GO Enrichment",
#'              zlim = c(0,15) )
#' @param grepList A list object where each value contains a group of search terms sperated by "|". The names and order values will be taken into account when plotting.
#' @param GOresults A list containing multiple GO$BP objects from GOanalysis
#' @param minT An integer specifying the minimum number of time a term must occur to be considered important. Default = 0
#' @param Rorder Row order for how the results shoul dbe plotted
#' @param mypalette Color palette used for heatmap
#' @param main Title of the plot
#' @param zlim Range specifying the p-values to be plotted
#' @import annotate
#' @import AnnotationDbi
#' @import GO.db
#' @import GSEABase
#' @import GOstats
#' @import WGCNA
#' @importMethodsFrom GOstats summary
#' @importMethodsFrom AnnotationDbi Term
#' @return Returns summary figure and summary table of results used to generate figure
#' @examples
#' #generate some data
#'      set1 <- sample(Ptrichocarpa_210_annotation_primary$ATG)[1:2000]
#'      set2 <- sample(Ptrichocarpa_210_annotation_primary$ATG)[1:2000]
#'
#' #calculate GO enrichment for BP
#'      GOresults<-list
#'      GOresults$set1<-atGOanalysis(set1,ontology="BP)
#'      GOresults$set2<-atGOanalysis(set2,ontology="BP)
#'
#' #create search list
#'      grepList<-list()
#'      grepList$hormone<-"hormone|gibberelli|brassino|auxin"
#'      grepList$peroxisome <- "peroxi"
#'      grepList$localization <-"protein localization"
#'
#' results<-ParseGOBPnPlot(grepList, GOresults)
#'
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
ParseGOBPnPlot<-function(grepList,GOresults, minT=0, Rorder=NULL,
                         my_palette=colorRampPalette(c("white", "red"))(n = 8),
                         main = "GO Enrichment",
                         zlim = c(0,15) )
{

  #check to make sure each GOresults has BP
  for(u in 1:length(GOresults))
  {
    if("BP" %ni% names(GOresults[[u]])) {stop(paste("BP missing from GOresults number", u,sep=" "))}
  }
    #check to make sure each GOresults has annotation
  for(y in 1:length(GOresults))
  {
  if("Term" %ni% names(GOstats::summary(GOresults[[y]]$BP))) {stop("Annotation Terms missing from GOresults: reload libraries annotate, GO.db, GSEABase and GOstats")}
  }

  key<- paste(unlist(grepList),collapse ="|")

  GOtable<-GOstats::summary(GOresults[[1]]$BP)[grep(key,GOstats::summary(GOresults[[1]]$BP)[,7],perl=TRUE),1:2]
  names(GOtable)[2]<-names(GOresults[1])

  for(l in 2:length(GOresults))
  {
    GOtmp<-GOstats::summary(GOresults[[l]]$BP)[grep(key,GOstats::summary(GOresults[[l]]$BP)[,7],perl=TRUE),1:2]
    names(GOtmp)[2]<-names(GOresults[l])
    GOtable<-merge(GOtable,GOtmp,by.x="GOBPID",by.y="GOBPID",all=T)
  }

  #remove terms that occur less than minimum times
  ind<-apply(GOtable,1,function(x) { y<-x[-1]; length(y[!is.na(y)]) })
  GOtable<-GOtable[which(ind>minT),]

  GOtable<-cbind(GOtable,AnnotationDbi::Term(GOtable$GOBPID))

  #order the GO terms based on the order of grepList
  o<-NULL
  v<-NULL
  for(G in 1:length(grepList))
  {
    oT<-grep(grepList[[G]],GOtable[,ncol(GOtable)])
    o<-c(o,oT)
    v<-c(v,length(oT))
  }

  #convert to -log10(pvalue)
  GOtable[is.na(GOtable)] <- 1
  end<-ncol(GOtable)-1
  GOtable[,2:end]<- -log10(GOtable[,2:end])


  results<-data.frame(t(GOtable[o,2:end]))
  results<-results[order(as.numeric(row.names(results))),]
  names(results)<-GOtable[o,1]

  #output results order
  if(!is.null(Rorder))
  {
    results<-results[Rorder,]
  }


  vlines<-cumsum(v)
  #make list of group names
  GOgroups<-rep(NA,ncol(results))
  for(e in 1:length(grepList))
  {
    if(e==1)
    {
      GOgroups[round(v[e]/2)]<-names(grepList)[e]
    } else {
      GOgroups[round(v[e]/2)+vlines[e-1]]<-names(grepList)[e]
    }
  }

  par(mar = c(9, 8, 2, 2));
  WGCNA::labeledHeatmap(Matrix = results,
                 xLabels = GOgroups,
                 yLabels = row.names(results),
                 yColorLabels = TRUE,
                 yColorWidth = 0.05,
                 ySymbols = row.names(results),
                 colors = my_palette,
                 #textMatrix = txt,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 cex.lab.y = 0.8,
                 cex.lab.x = 1,
                 xLabelsAngle = 90,
                 zlim = zlim,
                 main = main,
                 verticalSeparator.x = vlines
  )
  #return(GOtable)
  return(results)

}

#
#'Parse GO results using grep and plot
#'
#'This function creates a summary plot of GO$BP results by parsing GO results using grep and a list of search terms
#'
#' @usage ParseGOMFnPlot(grepList,GOresults, minT=0, Rorder=NULL,
#'              my_palette=colorRampPalette(c("white", "red"))(n = 8),
#'              main = "GO Enrichment",
#'              zlim = c(0,15) )
#' @param grepList A list object where each value contains a group of search terms sperated by "|". The names and order values will be taken into account when plotting.
#' @param GOresults A list containing multiple GO$MF objects from GOanalysis
#' @param minT An integer specifying the minimum number of time a term must occur to be considered important. Default = 0
#' @param Rorder Row order for how the results shoul dbe plotted
#' @param mypalette Color palette used for heatmap
#' @param main Title of the plot
#' @param zlim Range specifying the p-values to be plotted
#' @import annotate
#' @import AnnotationDbi
#' @import GO.db
#' @import GSEABase
#' @import GOstats
#' @import WGCNA
#' @importMethodsFrom GOstats summary
#' @importMethodsFrom AnnotationDbi Term
#' @return Returns summary figure and summary table of results used to generate figure
#' @examples
#' #generate some data
#'      set1 <- sample(Ptrichocarpa_210_annotation_primary$ATG)[1:2000]
#'      set2 <- sample(Ptrichocarpa_210_annotation_primary$ATG)[1:2000]
#'
#' #calculate GO enrichment for MF
#'      GOresults<-list
#'      GOresults$set1<-atGOanalysis(set1,ontology="MF)
#'      GOresults$set2<-atGOanalysis(set2,ontology="MF)
#'
#' #create search list
#'      grepList<-list()
#'      grepList$hormone<-"hormone|gibberelli|brassino|auxin"
#'      grepList$peroxisome <- "peroxi"
#'      grepList$localization <-"protein localization"
#'
#' results<-ParseGOMFnPlot(grepList, GOresults)
#'
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
ParseGOMFnPlot<-function(grepList,GOresults, minT=0, Rorder=NULL,
                         my_palette=colorRampPalette(c("white", "red"))(n = 8),
                         main = "GO Enrichment",
                         zlim = c(0,15) )
{

  #check to make sure each GOresults has MF
  for(u in 1:length(GOresults))
  {
    if("MF" %ni% names(GOresults[[u]])) {stop(paste("MF missing from GOresults number", u,sep=" "))}
  }
  #check to make sure each GOresults has annotation
  for(y in 1:length(GOresults))
  {
    if("Term" %ni% names(GOstats::summary(GOresults[[y]]$MF))) {stop("Annotation Terms missing from GOresults: reload libraries annotate, GO.db, GSEABase and GOstats")}
  }

  key<- paste(unlist(grepList),collapse ="|")

  GOtable<-GOstats::summary(GOresults[[1]]$MF)[grep(key,GOstats::summary(GOresults[[1]]$MF)[,7],perl=TRUE),1:2]
  names(GOtable)[2]<-names(GOresults[1])

  for(l in 2:length(GOresults))
  {
    GOtmp<-GOstats::summary(GOresults[[l]]$MF)[grep(key,GOstats::summary(GOresults[[l]]$MF)[,7],perl=TRUE),1:2]
    names(GOtmp)[2]<-names(GOresults[l])
    GOtable<-merge(GOtable,GOtmp,by.x="GOMFID",by.y="GOMFID",all=T)
  }

  #remove terms that occur less than minimum times
  ind<-apply(GOtable,1,function(x) { y<-x[-1]; length(y[!is.na(y)]) })
  GOtable<-GOtable[which(ind>minT),]

  GOtable<-cbind(GOtable,AnnotationDbi::Term(GOtable$GOMFID))

  #order the GO terms based on the order of grepList
  o<-NULL
  v<-NULL
  for(G in 1:length(grepList))
  {
    oT<-grep(grepList[[G]],GOtable[,ncol(GOtable)])
    o<-c(o,oT)
    v<-c(v,length(oT))
  }

  #convert to -log10(pvalue)
  GOtable[is.na(GOtable)] <- 1
  end<-ncol(GOtable)-1
  GOtable[,2:end]<- -log10(GOtable[,2:end])


  results<-data.frame(t(GOtable[o,2:end]))
  results<-results[order(as.numeric(row.names(results))),]
  names(results)<-GOtable[o,1]

  #output results order
  if(!is.null(Rorder))
  {
    results<-results[Rorder,]
  }


  vlines<-cumsum(v)
  #make list of group names
  GOgroups<-rep(NA,ncol(results))
  for(e in 1:length(grepList))
  {
    if(e==1)
    {
      GOgroups[round(v[e]/2)]<-names(grepList)[e]
    } else {
      GOgroups[round(v[e]/2)+vlines[e-1]]<-names(grepList)[e]
    }
  }

  par(mar = c(9, 8, 2, 2));
  WGCNA::labeledHeatmap(Matrix = results,
                        xLabels = GOgroups,
                        yLabels = row.names(results),
                        yColorLabels = TRUE,
                        yColorWidth = 0.05,
                        ySymbols = row.names(results),
                        colors = my_palette,
                        #textMatrix = txt,
                        setStdMargins = FALSE,
                        cex.text = 0.5,
                        cex.lab.y = 0.8,
                        cex.lab.x = 1,
                        xLabelsAngle = 90,
                        zlim = zlim,
                        main = main,
                        verticalSeparator.x = vlines
  )
  #return(GOtable)
  return(results)

}

#
#'Parse GO results using grep and plot
#'
#'This function creates a summary plot of GO$CC results by parsing GO results using grep and a list of search terms
#'
#' @usage ParseGOCCnPlot(grepList,GOresults, minT=0, Rorder=NULL,
#'              my_palette=colorRampPalette(c("white", "red"))(n = 8),
#'              main = "GO Enrichment",
#'              zlim = c(0,15) )
#' @param grepList A list object where each value contains a group of search terms sperated by "|". The names and order values will be taken into account when plotting.
#' @param GOresults A list containing multiple GO$CC objects from GOanalysis
#' @param minT An integer specifying the minimum number of time a term must occur to be considered important. Default = 0
#' @param Rorder Row order for how the results shoul dbe plotted
#' @param mypalette Color palette used for heatmap
#' @param main Title of the plot
#' @param zlim Range specifying the p-values to be plotted
#' @import annotate
#' @import AnnotationDbi
#' @import GO.db
#' @import GSEABase
#' @import GOstats
#' @import WGCNA
#' @importMethodsFrom GOstats summary
#' @importMethodsFrom AnnotationDbi Term
#' @return Returns summary figure and summary table of results used to generate figure
#' @examples
#' #generate some data
#'      set1 <- sample(Ptrichocarpa_210_annotation_primary$ATG)[1:2000]
#'      set2 <- sample(Ptrichocarpa_210_annotation_primary$ATG)[1:2000]
#'
#' #calculate GO enrichment for CC
#'      GOresults<-list
#'      GOresults$set1<-atGOanalysis(set1,ontology="CC)
#'      GOresults$set2<-atGOanalysis(set2,ontology="CC)
#'
#' #create search list
#'      grepList<-list()
#'      grepList$hormone<-"hormone|gibberelli|brassino|auxin"
#'      grepList$peroxisome <- "peroxi"
#'      grepList$localization <-"protein localization"
#'
#' results<-ParseGOCCnPlot(grepList, GOresults)
#'
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
ParseGOCCnPlot<-function(grepList,GOresults, minT=0, Rorder=NULL,
                         my_palette=colorRampPalette(c("white", "red"))(n = 8),
                         main = "GO Enrichment",
                         zlim = c(0,15) )
{

  #check to make sure each GOresults has CC
  for(u in 1:length(GOresults))
  {
    if("CC" %ni% names(GOresults[[u]])) {stop(paste("CC missing from GOresults number", u,sep=" "))}
  }
  #check to make sure each GOresults has annotation
  for(y in 1:length(GOresults))
  {
    if("Term" %ni% names(GOstats::summary(GOresults[[y]]$CC))) {stop("Annotation Terms missing from GOresults: reload libraries annotate, GO.db, GSEABase and GOstats")}
  }

  key<- paste(unlist(grepList),collapse ="|")

  GOtable<-GOstats::summary(GOresults[[1]]$CC)[grep(key,GOstats::summary(GOresults[[1]]$CC)[,7],perl=TRUE),1:2]
  names(GOtable)[2]<-names(GOresults[1])

  for(l in 2:length(GOresults))
  {
    GOtmp<-GOstats::summary(GOresults[[l]]$CC)[grep(key,GOstats::summary(GOresults[[l]]$CC)[,7],perl=TRUE),1:2]
    names(GOtmp)[2]<-names(GOresults[l])
    GOtable<-merge(GOtable,GOtmp,by.x="GOCCID",by.y="GOCCID",all=T)
  }

  #remove terms that occur less than minimum times
  ind<-apply(GOtable,1,function(x) { y<-x[-1]; length(y[!is.na(y)]) })
  GOtable<-GOtable[which(ind>minT),]

  GOtable<-cbind(GOtable,AnnotationDbi::Term(GOtable$GOCCID))

  #order the GO terms based on the order of grepList
  o<-NULL
  v<-NULL
  for(G in 1:length(grepList))
  {
    oT<-grep(grepList[[G]],GOtable[,ncol(GOtable)])
    o<-c(o,oT)
    v<-c(v,length(oT))
  }

  #convert to -log10(pvalue)
  GOtable[is.na(GOtable)] <- 1
  end<-ncol(GOtable)-1
  GOtable[,2:end]<- -log10(GOtable[,2:end])


  results<-data.frame(t(GOtable[o,2:end]))
  results<-results[order(as.numeric(row.names(results))),]
  names(results)<-GOtable[o,1]

  #output results order
  if(!is.null(Rorder))
  {
    results<-results[Rorder,]
  }


  vlines<-cumsum(v)
  #make list of group names
  GOgroups<-rep(NA,ncol(results))
  for(e in 1:length(grepList))
  {
    if(e==1)
    {
      GOgroups[round(v[e]/2)]<-names(grepList)[e]
    } else {
      GOgroups[round(v[e]/2)+vlines[e-1]]<-names(grepList)[e]
    }
  }

  par(mar = c(9, 8, 2, 2));
  WGCNA::labeledHeatmap(Matrix = results,
                        xLabels = GOgroups,
                        yLabels = row.names(results),
                        yColorLabels = TRUE,
                        yColorWidth = 0.05,
                        ySymbols = row.names(results),
                        colors = my_palette,
                        #textMatrix = txt,
                        setStdMargins = FALSE,
                        cex.text = 0.5,
                        cex.lab.y = 0.8,
                        cex.lab.x = 1,
                        xLabelsAngle = 90,
                        zlim = zlim,
                        main = main,
                        verticalSeparator.x = vlines
  )
  #return(GOtable)
  return(results)

}
