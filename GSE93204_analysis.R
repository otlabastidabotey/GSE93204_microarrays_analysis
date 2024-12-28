# Definim el directori de treball
setwd("C:/Users/Ot/Desktop/UOC/tfm/GSE93204_analisis")
# llibreries usades
library(GEOquery)
library(dplyr)
library(Biobase)
library(arrayQualityMetrics)
library(ggplot2)
library(ggrepel)
library(sva)
library(limma)
library(VennDiagram)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(UpSetR)




########## 1. CÀRREGA DE DADES
# Obtenim les dades que ens interessen
getGEO("GSE93204")
gse <- getGEO("GSE93204")[[1]]
# Simplifiquem el nom
phenoData <- pData(gse)
exprData <- exprs(gse)
# tenim informació sobre moltíssimes coses, però el que ens interessa realment
phenoData_summary <- phenoData %>%
  select(description.1,`treatment:ch2`, `patient id:ch2`) %>% 
  head(10)
phenoData_summary





########## 2. FILTRATGE  
# alguns pacients no tenen totes les condicions
# Vull només aquells que tinguin les 4
pacients_complets <- phenoData %>%
  filter(`treatment:ch2` %in% c('Baseline', 'C1D1', 'C1D15', 'Surgery')) %>%
  group_by(`patient id:ch2`) %>%
  filter(n() == 4) %>%
  ungroup() %>%
  distinct(`patient id:ch2`)

# Obtenim la llista dels complets
pacient_ids <- pacients_complets$`patient id:ch2`
# fem el filtratge de la phenodata
filtered_phenoData <- phenoData %>%
  filter(`patient id:ch2` %in% pacient_ids)
# filtratge de les dades d'expressió per només aquests pacients
filtered_exprData <- exprData[, colnames(exprData) %in% filtered_phenoData$geo_accession]
# no volem NAs
sum(is.na(filtered_exprData))
filtered_exprData<- apply(filtered_exprData, 2, function(x) { ifelse(is.na(x), mean(x, na.rm = TRUE), x) })
sum(is.na(filtered_exprData))

# Creem un nou objecte ExpressionSet amb les dades filtrades
# Necessitem una nova matriu de dades d'expressió 
# (es necessita com una matriu, no com un dataframe)
filtered_exprData <- as.matrix(filtered_exprData)
# Fem de la phenoData un AnnotatedDataFrame
filtered_phenoData <- AnnotatedDataFrame(filtered_phenoData)
# Creem l'ExpresionSet
gse <- ExpressionSet(assayData=filtered_exprData, phenoData=filtered_phenoData)
# Comprovem el nou ExpressionSet
print(gse)





########## 3. ANÀLISI DESCRIPTIU
# Amb la llibreria arrayQualityMetrics comprovem les nostres dades
arrayQualityMetrics(gse)
# També podem fer un boxplot d'intensitat, per veure l'expressió de les dades
boxplot(filtered_exprData, cex.axis=0.5, las=2,  which="all", 
        col = rep(c("chartreuse3", "salmon", "darkturquoise", "magenta4"), 8),
        main="Distribució de les intensitats")

# Fem també un PCA per veure la distribució de les nostres dades
# Definim la funció plotPCA3 (obtinguda de Gonzalo, Ricardo and Sanchez-Pla, Alex (2019))
plotPCA3 <- function (datos, labels, factor, title, scale, colores, size = 1.5, glineas = 0.25) {
  data <- prcomp(t(datos), scale = scale)
  dataDf <- data.frame(data$x)
  Grup <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100, 1)
  p1 <- ggplot(dataDf, aes(x = PC1, y = PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Grup), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[, 1]) - 5, max(data$x[, 1]) + 5)) +
    scale_fill_discrete(name = "Grup")
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels), segment.size = 0.25, size = size) + 
    labs(x = c(paste("PC1", loads[1], "%")), y = c(paste("PC2", loads[2], "%"))) +  
    ggtitle(paste("PCA:", title, sep = " ")) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colores)
}
# Executem la funció amb les dades i paràmetres definits
plotPCA3(filtered_exprData, filtered_phenoData$`patient id:ch2`, filtered_phenoData$`treatment:ch2`, "Mostres identificades per grup i pacient", TRUE, 
         colores = c("Baseline" = "chartreuse3", "C1D1" = "salmon", "C1D15" = "darkturquoise", "Surgery" = "magenta4"))

# No sembla un patró clar, tenim un clar efecte batch per pacient
# Corregim per dades aparellades
phenoData_df <- pData(filtered_phenoData)
# Definim el model de disseny per a la correcció 
mod <- model.matrix(~`treatment:ch2`, data=phenoData_df)
# Apliquem la correcció batch amb ComBat
corrected_exprData <- ComBat(dat=filtered_exprData, batch=phenoData_df$`patient id:ch2`, mod=mod)

# Visualització de PCA amb les dades corregides
plotPCA3(corrected_exprData, filtered_phenoData$`patient id:ch2`, filtered_phenoData$`treatment:ch2`, "Mostres corregides", TRUE, 
         colores = c("Baseline" = "chartreuse3", "C1D1" = "salmon", "C1D15" = "darkturquoise", "Surgery" = "magenta4"))


# Podem veure la variabilitat global de tots els gens
# Fem la funcó sd a tots els gens per saber-ne la desv. est.
sds <- apply (corrected_exprData, 1, sd)
sdsO<- sort(sds)
# ara ja fem el plot per veure la desviació dels gens
plot(1:length(sdsO), sdsO, main="Distribució de la variabilitat per a tots els gens",
     sub="Les línies verticals representen percentils del 90% i 95%",
     xlab="Índex del gen (de menys a més variable)", ylab="Desviació estàndard") + abline(v=length(sds)*c(0.9,0.95))


# Ara el que farem és refer el nostre ExpressionSet amb les dades corregides
# necessitem 'corrected_exprData' com a matriu
corrected_exprData <- as.matrix(corrected_exprData)
# Creem un AnnotatedDataFrame amb les dades fenotípiques
corrected_phenoData <- AnnotatedDataFrame(phenoData_df)
# Refem l'objecte ExpressionSet amb les dades corregides
gse <- ExpressionSet(assayData=corrected_exprData, phenoData=corrected_phenoData)
# Comprovem el nou ExpressionSet
print(gse)





########## 4. ANÀLISI D'EXPRESSIÓ
# Definim el disseny experimental amb les dades corregides
design <- model.matrix(~ 0 + phenoData_df$`treatment:ch2`)
colnames(design) <- c("Baseline", "C1D1", "C1D15", "Surgery")
# Ajustem el model lineal i anàlisi diferencial amb les dades corregides
fit <- lmFit(corrected_exprData, design)

# Definim els contrastos per comparar les condicions d'interès
contrast.matrix <- makeContrasts(
  Baseline_vs_C1D1 = C1D1 - Baseline,
  Baseline_vs_C1D15 = C1D15 - Baseline,
  Baseline_vs_Surgery = Surgery-Baseline,
  C1D1_vs_C1D15 = C1D15-C1D1,
  C1D1_vs_Surgery = Surgery - C1D1,
  C1D15_vs_Surgery = Surgery - C1D15,
  levels = design
)
# Ara ja podem ajustar el model segons les comparacions
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Resultats per a cada contrast
resultats_Baseline_vs_C1D1 <- topTable(fit2, coef = "Baseline_vs_C1D1", adjust = "fdr", number = Inf)
resultats_Baseline_vs_C1D15 <- topTable(fit2, coef = "Baseline_vs_C1D15", adjust = "fdr", number = Inf)
resultats_Baseline_vs_Surgery <- topTable(fit2, coef = "Baseline_vs_Surgery", adjust = "fdr", number = Inf)
resultats_C1D1_vs_C1D15 <- topTable(fit2, coef = "C1D1_vs_C1D15", adjust = "fdr", number = Inf)
resultats_C1D1_vs_Surgery <- topTable(fit2, coef = "C1D1_vs_Surgery", adjust = "fdr", number = Inf)
resultats_C1D15_vs_Surgery <- topTable(fit2, coef = "C1D15_vs_Surgery", adjust = "fdr", number = Inf)

# Podem veure el head d'un resultat
head(resultats_Baseline_vs_C1D15)
write.csv(resultats_Baseline_vs_C1D15, file = "resultats_Baseline_vs_C1D15.csv", row.names = TRUE)







########## 5. COMPARACIONS MÚLTIPLES
# per entendre com s'expressen els gens en les diferents condicions
res<-decideTests(fit2, method="separate", adjust.method="fdr", p.val=0.05, lfc=1)
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))

# De les sondes que han estat seleccionades per les comparacions múltiples
# Podem fer un HEATMAP
# Obtenim les sondes seleccionades per al heatmap
sondesHeatmap <- rownames(res.selected)
# Filtrem les dades d'expressió per les sondes seleccionades
dadesHM <- corrected_exprData[rownames(corrected_exprData) %in% sondesHeatmap,]
# Afegim informació dels grups
informacioGrups <- corrected_phenoData@data$`treatment:ch2`
names(informacioGrups) <- corrected_phenoData@data$geo_accession

# Ordre de les columnes de dadesHM segons timeline
columnes_ordenades <- order(informacioGrups)
dadesHM <- dadesHM[, columnes_ordenades]
informacioGrups <- informacioGrups[columnes_ordenades]
# Colors 
colorsPersonalitzats <- colorRampPalette(c("navyblue", "white", "red3"))(100)
# Ara ja podrem crear el heatmap 
pheatmap(dadesHM, 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         cluster_cols = FALSE, 
         annotation_col = data.frame(Grup = informacioGrups),
         color = colorsPersonalitzats,
         show_rownames = FALSE)


# Per saber en detall quines són les sondes destacades entre grups
# Fem un diagrama de Venn
venn.plot <- venn.diagram(
  x = list(
    "Baseline vs C1D15" = which(res.selected[, 2] != 0),
    "Baseline vs Surgery" = which(res.selected[, 3] != 0),
    "C1D1 vs C1D15" = which(res.selected[, 4] != 0),
    "C1D15 vs Surgery" = which(res.selected[, 6] != 0)
  ),
  category.names = c("Baseline vs C1D15", "Baseline vs Surgery", "C1D1 vs C1D15", "C1D15 vs Surgery"),
  filename = NULL,
  fill = c("cornflowerblue", "forestgreen", "red2", "orange"), # Colors de les àrees
  cex = 1,
  cat.cex = 0.7,
  cat.pos = 0,
  main = "Gens en comú en les 4 comparacions",
  fontfamily = "sans", # Tipus de lletra per als nombres
  cat.fontfamily = "sans",
  main.fontfamily = "sans"
)
grid.draw(venn.plot)

# Podem saber les sondes significatives per cada comparació
Baseline_C1D1 = rownames(res.selected)[res.selected[, 1] != 0]
Baseline_C1D15 = rownames(res.selected)[res.selected[, 2] != 0]
Baseline_Surgery = rownames(res.selected)[res.selected[, 3] != 0]
C1D1_C1D15 = rownames(res.selected)[res.selected[, 4] != 0]
C1D1_Surgery = rownames(res.selected)[res.selected[, 5] != 0]
C1D15_Surgery = rownames(res.selected)[res.selected[, 6] != 0]

# Volem saber els gens comuns entre les comparacions
sondes_comunes <- Reduce(intersect, list(Baseline_C1D15, Baseline_Surgery, C1D1_C1D15, C1D15_Surgery))
# Anotem els noms dels gens comuns
print(sondes_comunes)
# Tenim un fitxer que conté la conversió entre els aligent id i els gene symbol
taula_conversio <- read.table("C:/Users/Ot/Desktop/UOC/tfm/GSE93204_analisis/GPL6480-9577.txt", 
                              header=TRUE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
# Les columnes que necessitem: ID, GENE, GENE_SYMBOL, ENSEMBL_ID
taula_conversio <- taula_conversio %>% select(ID, GENE, GENE_SYMBOL, ENSEMBL_ID)
# Podem saber a què fan referencia les sondes compartides per les 4 comparacions
genes_comuns_annotats <- taula_conversio[taula_conversio$ID %in% sondes_comunes, ]
genes_comuns_annotats










########## 6. ANOTACIÓ DE LES SONDES
# Podem veure que hi ha sondes que no fan referència a cap gen
# I dues sondes poden fer referència al mateix gen
# Per tant, és important fer l'anotació dels resultats
# Anotem els resultats topTable amb gene symbols, gene i ensemblID
annotate_topTable <- function(topTable, taula_conversio) {
  topTable$GENE_SYMBOL <- taula_conversio$GENE_SYMBOL[match(rownames(topTable), taula_conversio$ID)]
  topTable$GENE <- taula_conversio$GENE[match(rownames(topTable), taula_conversio$ID)]
  topTable$ENSEMBL_ID <- taula_conversio$ENSEMBL_ID[match(rownames(topTable), taula_conversio$ID)]
  return(topTable)
}
# així tenim cada topTable amb les noves columnes
anot_Baseline_vs_C1D1 <- annotate_topTable(resultats_Baseline_vs_C1D1, taula_conversio)
anot_Baseline_vs_C1D15 <- annotate_topTable(resultats_Baseline_vs_C1D15, taula_conversio)
anot_Baseline_vs_Surgery <- annotate_topTable(resultats_Baseline_vs_Surgery, taula_conversio)
anot_C1D1_vs_C1D15 <- annotate_topTable(resultats_C1D1_vs_C1D15, taula_conversio)
anot_C1D1_vs_Surgery <- annotate_topTable(resultats_C1D1_vs_Surgery, taula_conversio)
anot_C1D15_vs_Surgery <- annotate_topTable(resultats_C1D15_vs_Surgery, taula_conversio)

# Podem veure el head d'un resultat (alguns NA)
head(anot_Baseline_vs_C1D15)
# Per tant, haurem de netejar la taula
# Eliminar les files amb GENE_SYMBOL NA o en blanc, duplicats i assignar rownames
resultats_nets <- function(resultats_annot) {
  resultats_nets <- resultats_annot[!is.na(resultats_annot$GENE_SYMBOL) & resultats_annot$GENE_SYMBOL != "", ]
  resultats_nets <- resultats_nets[!duplicated(resultats_nets$GENE_SYMBOL), ]
  rownames(resultats_nets) <- resultats_nets$GENE_SYMBOL
  return(resultats_nets)
}

# Així, podem netejar les taules d'anotacions
Baseline_vs_C1D1_net <- resultats_nets(anot_Baseline_vs_C1D1)
Baseline_vs_C1D15_net <- resultats_nets(anot_Baseline_vs_C1D15)
Baseline_vs_Surgery_net <- resultats_nets(anot_Baseline_vs_Surgery)
C1D1_vs_C1D15_net <- resultats_nets(anot_C1D1_vs_C1D15)
C1D1_vs_Surgery_net <- resultats_nets(anot_C1D1_vs_Surgery)
C1D15_vs_Surgery_net <- resultats_nets(anot_C1D15_vs_Surgery)

# Podem observar el resultat de la neteja de la taula
head(Baseline_vs_C1D15_net)
write.csv(Baseline_vs_C1D15_net, file = "Baseline_vs_C1D15_net.csv", row.names = TRUE)






########## 7. COMPARACIONS INDIVIDUALS ENTRE GRUPS
# Creem el volcano plot per a Baseline vs C1D1
volcano_Baseline_vs_C1D1 <- ggplot(Baseline_vs_C1D1_net, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4) +
  theme_classic() +
  labs(title = "Volcano Plot: Baseline vs C1D1", x = "Log Fold Change", y = "-Log10 P-value") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed")
# Veiem els gens diferencialment expressats
print(volcano_Baseline_vs_C1D1)

# Creem el volcano plot per a Baseline vs C1D15
volcano_Baseline_vs_C1D15 <- ggplot(Baseline_vs_C1D15_net, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4) +
  theme_classic() +
  labs(title = "Volcano Plot: Baseline vs C1D15", x = "Log Fold Change", y = "-Log10 P-value") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed")
# Veiem els gens diferencialment expressats
print(volcano_Baseline_vs_C1D15)

# Creem el volcano plot per a Baseline vs Surgery
volcano_Baseline_vs_Surgery <- ggplot(Baseline_vs_Surgery_net, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4) +
  theme_classic() +
  labs(title = "Volcano Plot: Baseline vs Surgery", x = "Log Fold Change", y = "-Log10 P-value") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed")
# Veiem els gens diferencialment expressats
print(volcano_Baseline_vs_Surgery)

# Creem el volcano plot per a C1D1 vs C1D15
volcano_C1D1_vs_C1D15 <- ggplot(C1D1_vs_C1D15_net, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4) +
  theme_classic() +
  labs(title = "Volcano Plot: C1D1 vs C1D15", x = "Log Fold Change", y = "-Log10 P-value") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed")
# Veiem els gens diferencialment expressats
print(volcano_C1D1_vs_C1D15)

# Creem el volcano plot per a C1D1 vs Surgery
volcano_C1D1_vs_Surgery <- ggplot(C1D1_vs_Surgery_net, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4) +
  theme_classic() +
  labs(title = "Volcano Plot: C1D1 vs Surgery", x = "Log Fold Change", y = "-Log10 P-value") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed")
# Veiem els gens diferencialment expressats
print(volcano_C1D1_vs_Surgery)

# Creem el volcano plot per a C1D15 vs Surgery
volcano_C1D15_vs_Surgery <- ggplot(C1D15_vs_Surgery_net, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4) +
  theme_classic() +
  labs(title = "Volcano Plot: C1D15 vs Surgery", x = "Log Fold Change", y = "-Log10 P-value") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed")
# Veiem els gens diferencialment expressats
print(volcano_C1D15_vs_Surgery)





########## 8. SIGNIFICACIÓ BIOLÒGICA
# Podem fer un GSEA
# El primer que necessitem és definir les llistes de gens per cada comparació
llista_baseline_vs_c1d1 <- Baseline_vs_C1D1_net$logFC
names(llista_baseline_vs_c1d1) <- rownames(Baseline_vs_C1D1_net)
llista_baseline_vs_c1d1 <- sort(llista_baseline_vs_c1d1, decreasing = TRUE)

llista_baseline_vs_c1d15 <- Baseline_vs_C1D15_net$logFC
names(llista_baseline_vs_c1d15) <- rownames(Baseline_vs_C1D15_net)
llista_baseline_vs_c1d15 <- sort(llista_baseline_vs_c1d15, decreasing = TRUE)

llista_baseline_vs_surgery <- Baseline_vs_Surgery_net$logFC
names(llista_baseline_vs_surgery) <- rownames(Baseline_vs_Surgery_net)
llista_baseline_vs_surgery <- sort(llista_baseline_vs_surgery, decreasing = TRUE)

llista_c1d1_vs_c1d15 <- C1D1_vs_C1D15_net$logFC
names(llista_c1d1_vs_c1d15) <- rownames(C1D1_vs_C1D15_net)
llista_c1d1_vs_c1d15 <- sort(llista_c1d1_vs_c1d15, decreasing = TRUE)

llista_c1d1_vs_surgery <- C1D1_vs_Surgery_net$logFC
names(llista_c1d1_vs_surgery) <- rownames(C1D1_vs_Surgery_net)
llista_c1d1_vs_surgery <- sort(llista_c1d1_vs_surgery, decreasing = TRUE)

llista_c1d15_vs_surgery <- C1D15_vs_Surgery_net$logFC
names(llista_c1d15_vs_surgery) <- rownames(C1D15_vs_Surgery_net)
llista_c1d15_vs_surgery <- sort(llista_c1d15_vs_surgery, decreasing = TRUE)

# Ara ja podem realitzar cada GSEA amb els resultats net
gsea_baseline_vs_c1d1 <- gseGO(geneList = llista_baseline_vs_c1d1, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                               ont = "ALL", minGSSize = 15, maxGSSize = 500, pAdjustMethod = "fdr")
gsea_baseline_vs_c1d15 <- gseGO(geneList = llista_baseline_vs_c1d15, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                                ont = "ALL", minGSSize = 15, maxGSSize = 500, pAdjustMethod = "fdr")
gsea_baseline_vs_surgery <- gseGO(geneList = llista_baseline_vs_surgery, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                                  ont = "ALL", minGSSize = 15, maxGSSize = 500, pAdjustMethod = "fdr")
gsea_c1d1_vs_c1d15 <- gseGO(geneList = llista_c1d1_vs_c1d15, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                            ont = "ALL", minGSSize = 15, maxGSSize = 500, pAdjustMethod = "fdr")
gsea_c1d1_vs_surgery <- gseGO(geneList = llista_c1d1_vs_surgery, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                              ont = "ALL", minGSSize = 15, maxGSSize = 500, pAdjustMethod = "fdr")
gsea_c1d15_vs_surgery <- gseGO(geneList = llista_c1d15_vs_surgery, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                               ont = "ALL", minGSSize = 15, maxGSSize = 500, pAdjustMethod = "fdr")
# Podem veure els termes més enriquits per cada comparació
dotplot(gsea_baseline_vs_c1d1) + ggtitle("GSEA: Baseline vs C1D1")
dotplot(gsea_baseline_vs_c1d15) + ggtitle("GSEA: Baseline vs C1D15")
dotplot(gsea_baseline_vs_surgery) + ggtitle("GSEA: Baseline vs Surgery")
dotplot(gsea_c1d1_vs_c1d15) + ggtitle("GSEA: C1D1 vs C1D15")
dotplot(gsea_c1d1_vs_surgery) + ggtitle("GSEA: C1D1 vs Surgery")
dotplot(gsea_c1d15_vs_surgery) + ggtitle("GSEA: C1D15 vs Surgery")


# Potser per alguna comparació seria interessant fer
cnetplot(gsea_baseline_vs_c1d15) + ggtitle("Cnetplot: Baseline vs C1D15")




########## 9. GENS CANDIDATS DE RESISTÈNCIA
# El que sí que volem saber són quins gens estan sobreexpressats
# Filtrem els gens sobreexpressats per a cada comparació utilitzant adjusted p-value
gens_Baseline_vs_C1D1_UP <- rownames(subset(Baseline_vs_C1D1_net, (logFC > 1) & (adj.P.Val < 0.05)))
gens_Baseline_vs_C1D15_UP <- rownames(subset(Baseline_vs_C1D15_net, (logFC > 1) & (adj.P.Val < 0.05)))
gens_Baseline_vs_Surgery_UP <- rownames(subset(Baseline_vs_Surgery_net, (logFC > 1) & (adj.P.Val < 0.05)))
gens_C1D1_vs_C1D15_UP <- rownames(subset(C1D1_vs_C1D15_net, (logFC > 1) & (adj.P.Val < 0.05)))
gens_C1D1_vs_Surgery_UP <- rownames(subset(C1D1_vs_Surgery_net, (logFC > 1) & (adj.P.Val < 0.05)))
gens_C1D15_vs_Surgery_UP <- rownames(subset(C1D15_vs_Surgery_net, (logFC > 1) & (adj.P.Val < 0.05)))
# Mostrem els resultats
gens_Baseline_vs_C1D1_UP
gens_Baseline_vs_C1D15_UP
gens_Baseline_vs_Surgery_UP
gens_C1D1_vs_C1D15_UP
gens_C1D1_vs_Surgery_UP
gens_C1D15_vs_Surgery_UP

# Gens comuns entre Baseline_vs_C1D1 i Baseline_vs_C1D15
comuns_baseline_c1d1_c1d15 <- intersect(gens_Baseline_vs_C1D1_UP, gens_Baseline_vs_C1D15_UP)
print(comuns_baseline_c1d1_c1d15)


# Gens comuns entre Baseline_vs_Surgery i gens_C1D1_vs_Surgery_UP
comuns_c1d1_c1d15_surgery <- intersect(gens_C1D1_vs_Surgery_UP, gens_C1D15_vs_Surgery_UP)
print(comuns_c1d1_c1d15_surgery)


# Per veure les interseccions entre els gens sobreexpressats
# Primer fem una llista amb tots els gens que tenim
gens_total <- unique(c(
  gens_Baseline_vs_C1D1_UP,
  gens_Baseline_vs_C1D15_UP,
  gens_Baseline_vs_Surgery_UP,
  gens_C1D1_vs_C1D15_UP,
  gens_C1D1_vs_Surgery_UP,
  gens_C1D15_vs_Surgery_UP
))

# Crear una taula booleana TRUE/FALSE 
# així sabem si  el gen està, o no, en cada comparació
dades_upset <- data.frame(
  Gen = gens_total,
  `Baseline vs C1D1` = gens_total %in% gens_Baseline_vs_C1D1_UP,
  `Baseline vs C1D15` = gens_total %in% gens_Baseline_vs_C1D15_UP,
  `Baseline vs Surgery` = gens_total %in% gens_Baseline_vs_Surgery_UP,
  `C1D1 vs C1D15` = gens_total %in% gens_C1D1_vs_C1D15_UP,
  `C1D1 vs Surgery` = gens_total %in% gens_C1D1_vs_Surgery_UP,
  `C1D15 vs Surgery` = gens_total %in% gens_C1D15_vs_Surgery_UP
)

# Convertir TRUE/FALSE a 1/0 per totes les columnes menys la columna 'Gen'
# Ho necessitem per la funció upset
dades_upset <- dades_upset %>% mutate(across(-Gen, as.integer))
# a més, em canvia en nom de les columnes, no volem punts, vull espais
colnames(dades_upset) <- gsub("\\.", " ", colnames(dades_upset))

# Fem el plot
upset(dades_upset,
  sets = c("Baseline vs C1D1", "Baseline vs C1D15", "Baseline vs Surgery",
           "C1D1 vs C1D15", "C1D1 vs Surgery", "C1D15 vs Surgery"),
  order.by = "freq")


#### Gens d'interès
gens_interes <- unique(c(
  gens_C1D1_vs_Surgery_UP,
  gens_C1D15_vs_Surgery_UP
))



###############  COMPARACIÓ AMB INVITRO

# Dades de gens BT-474 i ZR75 CRISPR/CAS9 SCREEN CANDIDATES
gens_bt474 <- c("ACVRL1", "ID1", "PLCB1", "ITGB6", "TP53", "ADAM12", "PAX5", "THBS4", "TIMP4", "PDCD1", 
                "COL27A1", "CD8A", "SNAI1", "NSD1", "DUSP4", "SCARA5", "B3GNT3", "PYCARD", "TLR4", "LEF1",
                "CDKN1A", "PIK3CA", "AGT", "TTYH1", "CACNG1", "LTB", "LEPR", "TSPAN7", "MAPK1", "SMC1B", 
                "HSPA2", "MSR1", "LFNG", "TCF4", "PIK3R2", "HOXA9", "PSMB10", "GJB2", "BLM", "CNTFR", "LAMA3",
                "PRF1", "COLEC12", "EFNA5", "SOCS1", "EPAS1", "PAX8", "ITPR1", "DTX3", "DKK2", "RASGRP1", 
                "CCNE1", "LTBP1", "GRIN1", "NEIL2", "PARP1", "GNLY", "BNIP3", "HLA-DQB1", "AREG", "HLA-DMB", 
                "BBOX1", "PDE9A", "E2F1", "TUBA4A", "DLL3", "HLA-B", "TPSAB1", "EGLN3", "IL22RA2", "SLC2A11", 
                "PTEN", "FSTL1", "BTG2", "NETO2", "CCL5", "IL1B", "RORA", "SFRP1", "CDKN3", "GSK3B", "MT1G", 
                "DDR2", "CD27", "LAMB3", "MMP3", "STAT1", "CPA3", "TSPAN1", "NOTCH3", "HLA-E", "RPS6KA5", 
                "CCND2", "NFATC1", "FZD9", "IFT140", "CD19", "CCR5", "BMP8A", "WNT11", "HIF1A", "IDO1", 
                "IL2RA")

gens_zr75 <- c("BMP6", "TCF4", "GADD45A", "DSC2", "HDC", "BCAS1", "PIK3R3", "SOCS2", "ATP10B", "PIK3R2",
               "PIK3CG", "SMAD3", "PDCD1", "STC1", "CAV1", "CSF3R", "CDH1", "ID4", "CDC25A", "CCL21",
               "PDK4", "LEF1", "CD19", "CXCL5", "LTB", "ECM2", "POLQ", "GJB2", "PSMB9", "BMP4", "TNKS2",
               "NOD2", "ERBB4", "DUSP4", "LIF", "AGTR1", "ROCK2", "ADAM12", "ATOSA", "PTTG1", "PBX3",
               "CDH2", "HDAC11", "TNFAIP6", "CLDN3", "MUS81", "PRKACA", "SMURF2", "HBB", "DHRS2", "FGL2",
               "CKS1B", "RASGRF2", "TYK2", "DLGAP5", "CCR5", "FHL1", "GZMB")

gens_invitro <- union(gens_zr75, gens_bt474)
gens_invitro


gens_interseccio_microarrays<-intersect(gens_interes,gens_invitro)
gens_interseccio_microarrays





########## què passa si ara ho faig amb p.val 0.1
# Filtrem els gens sobreexpressats per a cada comparació utilitzant adjusted p-value
gens_Baseline_vs_C1D1_UP_0.1 <- rownames(subset(Baseline_vs_C1D1_net, (logFC > 1) & (adj.P.Val < 0.1)))
gens_Baseline_vs_C1D15_UP_0.1 <- rownames(subset(Baseline_vs_C1D15_net, (logFC > 1) & (adj.P.Val < 0.1)))
gens_Baseline_vs_Surgery_UP_0.1 <- rownames(subset(Baseline_vs_Surgery_net, (logFC > 1) & (adj.P.Val < 0.1)))
gens_C1D1_vs_C1D15_UP_0.1 <- rownames(subset(C1D1_vs_C1D15_net, (logFC > 1) & (adj.P.Val < 0.1)))
gens_C1D1_vs_Surgery_UP_0.1 <- rownames(subset(C1D1_vs_Surgery_net, (logFC > 1) & (adj.P.Val < 0.1)))
gens_C1D15_vs_Surgery_UP_0.1 <- rownames(subset(C1D15_vs_Surgery_net, (logFC > 1) & (adj.P.Val < 0.1)))
# Mostrem els resultats
gens_Baseline_vs_C1D1_UP_0.1
gens_Baseline_vs_C1D15_UP_0.1
gens_Baseline_vs_Surgery_UP_0.1
gens_C1D1_vs_C1D15_UP_0.1
gens_C1D1_vs_Surgery_UP_0.1
gens_C1D15_vs_Surgery_UP_0.1

# Per veure les interseccions entre els gens sobreexpressats
# Primer fem una llista amb tots els gens que tenim
gens_total_0.1 <- unique(c(
  gens_Baseline_vs_C1D1_UP_0.1,
  gens_Baseline_vs_C1D15_UP_0.1,
  gens_Baseline_vs_Surgery_UP_0.1,
  gens_C1D1_vs_C1D15_UP_0.1,
  gens_C1D1_vs_Surgery_UP_0.1,
  gens_C1D15_vs_Surgery_UP_0.1
))

# Crear una taula booleana TRUE/FALSE 
# així sabem si  el gen està, o no, en cada comparació
dades_upset_0.1 <- data.frame(
  Gen = gens_total_0.1,
  `Baseline vs C1D1` = gens_total_0.1 %in% gens_Baseline_vs_C1D1_UP_0.1,
  `Baseline vs C1D15` = gens_total_0.1 %in% gens_Baseline_vs_C1D15_UP_0.1,
  `Baseline vs Surgery` = gens_total_0.1 %in% gens_Baseline_vs_Surgery_UP_0.1,
  `C1D1 vs C1D15` = gens_total_0.1 %in% gens_C1D1_vs_C1D15_UP_0.1,
  `C1D1 vs Surgery` = gens_total_0.1 %in% gens_C1D1_vs_Surgery_UP_0.1,
  `C1D15 vs Surgery` = gens_total_0.1 %in% gens_C1D15_vs_Surgery_UP_0.1
)

# Convertir TRUE/FALSE a 1/0 per totes les columnes menys la columna 'Gen'
# Ho necessitem per la funció upset
dades_upset_0.1 <- dades_upset_0.1 %>% mutate(across(-Gen, as.integer))
# a més, em canvia en nom de les columnes, no volem punts, vull espais
colnames(dades_upset_0.1) <- gsub("\\.", " ", colnames(dades_upset_0.1))

# Fem el plot
upset(dades_upset_0.1,
      sets = c("Baseline vs C1D1", "Baseline vs C1D15", "Baseline vs Surgery",
               "C1D1 vs C1D15", "C1D1 vs Surgery", "C1D15 vs Surgery"),
      order.by = "freq")


#### Gens d'interès
gens_interes_0.1 <- unique(c(
  gens_C1D1_vs_Surgery_UP_0.1,
  gens_C1D15_vs_Surgery_UP_0.1
))














