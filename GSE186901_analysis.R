# Definim el directori de treball
setwd("C:/Users/Ot/Desktop/UOC/tfm/GSE186901_analisis")

# Definim les llibreries que necessitarem per a l'anàlisi
library(GO.db)
library(GEOquery)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)
library(randomcoloR)
library(limma)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)


##### 1. CÀRREGA DE LES DADES
# Obtenim les dades fenotípiques de GEO
gse <- getGEO("GSE186901", GSEMatrix = TRUE)[[1]]
# N'extraiem la informació fenotípica
pheno_data <- pData(gse)
print(head(pheno_data))
dim(pheno_data)
# Carreguem la matriu de comptatges RNA-Seq
# No la podem obtenir directament del GSE
seqdata <- read.delim("./GSE186901_palbo_geo_tpm_mat.txt", stringsAsFactors = FALSE)
# La primera columna només indica el número d'entrada (1,2,3...)
countdata <- seqdata[,-1]  
dim(countdata)
# Anotació a Entrez ID per poder fer l'anàlisi
# Obtenim els Entrez IDs 
gens_ids <- countdata$gene_name
# BiomaRt molts problemes, ho fem amb AnnotationDbi
# important determinar que select és d'aquest paquet
anotacions <- AnnotationDbi::select(org.Hs.eg.db, keys = gens_ids, columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")
anotacions
# Ajuntem els Entrez IDs amb la matriu de RNA-Seq
countdata <- merge(countdata, anotacions, by.x = "gene_name", by.y = "SYMBOL")
dim(countdata)
# si ho fem amb annotationnDBI
colnames(countdata)[colnames(countdata) == "ENTREZID"] <- "entrezgene_id"




########## 2. FILTRATGE DE LES DADES
# Eliminem els les files amb Entrez IDs duplicats
# No volem 2 sondes pel mateix gen
duplicats <- countdata$entrezgene_id[duplicated(countdata$entrezgene_id)]
countdata <- countdata[!countdata$entrezgene_id %in% duplicats, ]
# Assignem els entrez IDs com a noms de les files
rownames(countdata) <- countdata$entrezgene_id
# Eliminem les columnes que no tenen info numèrica
# l'anomenem counts.TPM per recordar que és en TPM
counts.TPM <- countdata[, -which(colnames(countdata) %in% c("gene_name", "accession_number","entrezgene_id"))]
dim(counts.TPM)


##### FILTRATGE PER NOMÉS DADES APARELLADES
# Obtenim els noms dels pacients que té el GSE
noms_pacients <- colnames(counts.TPM)
# Extraiem l'identificador del pacient de les mostres
id_pacients <- sub("_.*", "", noms_pacients)
# Identifiquem NOMÉS els pacients amb mostres tant en Baseline com en PD
baseline_samples <- grep("Baseline", noms_pacients, value = TRUE)
pd_samples <- grep("PD", noms_pacients, value = TRUE)
pacients_seleccionats <- intersect(sub("_.*", "", baseline_samples), sub("_.*", "", pd_samples))
# del pacient 66 no tenim tota la informació, l'eliminem
pacients_seleccionats <- pacients_seleccionats[pacients_seleccionats != "BRO7F.066"]
# Seleccionem només els comptatges dels pacients amb dades tant de Baseline com de PD
samples_seleccionats <- noms_pacients[id_pacients %in% pacients_seleccionats]
counts.TPM <- counts.TPM[, samples_seleccionats]
# Les mostres seleccionades (44, de 22 pacients)
print(samples_seleccionats)
length(samples_seleccionats)
# També hem de filtrar la phenodata
# Com que el format no és el mateix, hem de canviar punts per guions en els noms de les mostres 
samples_seleccionats <- gsub("\\.", "-", samples_seleccionats)
# Filtrar pheno_data per les mostres seleccionades utilitzant %in%
pheno_data <- pheno_data[pheno_data$title %in% samples_seleccionats, ]
# Comrpovem que tenim les 44 mostres (amb 42 dades clíniques)
dim(pheno_data)



##### 2. OBERTURA DE LES DADES (sense filtratge per dades aparellades)

# Obtenim els noms dels pacients que té el GSE
###########noms_pacients <- colnames(counts.TPM)

# Ja no filtrem les mostres per pacients aparellats
# Comprovem que tenim totes les mostres de les dades
###########3print(noms_pacients)

# A continuació, no fem cap filtre sobre les mostres.
# Mantenim totes les mostres que tenim
# Comprovem que tenim totes les mostres
########dim(pheno_data)






##### 3. PREPROCESSAT DE LES DADES (I)
# Per facilitar l'anàlisi, ens quedem només amb els gens més expressats
thresh <- counts.TPM > 0.5
# Obtenim TRUE o FALSE
head(thresh)
# Ens quedem els gens que com a mínim tinguin alta expressió en 2 mostres
keep <- rowSums(thresh) >= 2
# Filtrem per aquest processament i ens quedem la matriu 'counts.keep'
counts.keep <- counts.TPM[keep,]
# Veiem que tenim les 44 mostres però ara comptem amb menys gens
dim(counts.keep)




########## 4. ANÀLISI DE QUALITAT
# Per l'anàlisi serà més fàcil amb conversió log2
logcounts <- as.matrix(log2(counts.TPM+0.1))

# Ajustem els noms (són massa llargs)
ajustar_noms_mostres <- function(noms_mostres) {
  noms_mostres_modificats <- sub("Baseline_", "B", noms_mostres)
  noms_mostres_modificats <- sub("PD_", "PD", noms_mostres_modificats)
  noms_mostres_modificats <- sub("WTS", "", noms_mostres_modificats)
  return(noms_mostres_modificats)
}
# cannviem el noms a la matriu de comptatges
colnames(logcounts) <- ajustar_noms_mostres(colnames(logcounts))
# afegim colors per grup
colors <- rep(c("darkorchid1", "darkgreen"), length.out = ncol(logcounts))

# Distribució de tots els nostres comptatges
par(mar = c(7, 4, 2, 2))
boxplot(logcounts, ylab="Log2-TPM", las=2, xlab="", col=colors, cex.axis=0.8, main="Boxplots de logTPMs")
abline(h=median(logcounts), col="blue")


# Per la notació, ho deixem com ho teníem
par(mar = c(5, 4, 4, 2))
logcounts <- as.matrix(log2(counts.TPM+0.1))



########## 5. ANÀLISI DESCRIPTIVA

##  Podem fer un PCA per presentar les nostres dades
pca_result <- prcomp(t(logcounts), scale. = FALSE)
# Convertim els resultats a un data frame
pca_data <- data.frame(pca_result$x)
# Incloem la informació fenotípica
pca_data$grup <- ifelse(grepl("Baseline", rownames(pca_data)), "Baseline", "PD")
pca_data$pacient <- sub("_.*", "", rownames(pca_data))

# PCA amb colors diferents per a 'Baseline' i 'PD'
pca_baseline_pd <- ggplot(pca_data, aes(x = PC1, y = PC2, color = grup)) +
  geom_point(size = 3, alpha = 0.6) +
  theme_classic() +
  labs(title = "PCA de les mostres en les dues condicions", x = paste0("PC1: ", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 1), "% variància"),
       y = paste0("PC2: ", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 1), "% variància")) +
  scale_color_manual(values = c("Baseline" = "darkorchid1", "PD" = "darkgreen"))

print(pca_baseline_pd)

# No hi ha un patró clar, veiem com s'agruppen entre pacients
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = grup)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_line(aes(group = pacient), color = "gray", alpha = 0.5) +
  theme_classic() +
  labs(title = "PCA connectant mostres per pacients",
       x = paste0("PC1: ", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 1), "% variància"),
       y = paste0("PC2: ", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 1), "% variància")) +
  scale_color_manual(values = c("Baseline" = "darkorchid1", "PD" = "darkgreen"))

print(pca_plot)


# Per expemplificar-ho, fem PCA per pacients (necessitem més colors)
# Per 44 mostres, necessitem 2 colors distintius
colors <- distinctColorPalette(23)
# PCA amb colors diferents per a cada pacient
pca_pacients <- ggplot(pca_data, aes(x = PC1, y = PC2, color = pacient)) +
  geom_point(size = 3, alpha = 0.6) +
  theme_classic() +
  labs(title = "PCA per Pacients", x = paste0("PC1: ", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 1), "% variança"),
       y = paste0("PC2: ", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 1), "% variança")) +
  scale_color_manual(values = colors)

print(pca_pacients)


# Podem veure-ho també veient-ne les diferències
ggplot(pca_data, aes(x = grup, y = PC1, group = pacient)) +
  geom_point(aes(color = pacient), size = 3, alpha = 0.6) +
  geom_line(color = "grey", alpha = 0.5) +
  theme_classic() +
  labs(title = "Canvi en PC1 entre Baseline i PD",
       x = "Condició", y = "PC1") +
  scale_color_manual(values = colors)

# I per representar les diferències, podem veure-ho en colors, també
pca_data <- pca_data %>%
  arrange(pacient, grup) %>%
  group_by(pacient) %>%
  mutate(direccio = ifelse(lead(PC1) > PC1, "asc", "desc")) %>%
  ungroup()

# Podem veure com els canvis són majoritàriament decreixents
ggplot(pca_data, aes(x = grup, y = PC1, group = pacient)) +
  geom_point(aes(color = pacient), size = 3, alpha = 0.6) +
  geom_line(aes(color = direccio), size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("asc" = "green", "desc" = "red")) +
  theme_classic() +
  labs(title = "Canvi en PC1 entre Baseline i PD",
       x = "Grup", y = "PC1")





############# 6. DISSENY DEL MODEL 
# Creem una variable que identifiqui cada pacient
pheno_data$pacient <- sub("_.*", "", pheno_data$title)
# i una que identifiqui cada grup
pheno_data$grup <- ifelse(grepl("Baseline", pheno_data$title), "Baseline", "PD")
# definim el model que analitzarem, segons el grup
design <- model.matrix(~ 0 + grup, data = pheno_data)
colnames(design) <- c("Baseline", "PD")
design


## això és el que teniem abans de corregir per mostres aparellades
# Definim el grup basant-se en time:ch1 (Baseline vs PD)
###### group <- pheno_data$`time:ch1`
# matriu de disseny
###### design <- model.matrix(~0 + group)
###### colnames(design) <- c("Baseline","PD")
###### rownames(design) <- pheno_data$title
###### print(design)
# voom per transformar les dades
###### dge <- DGEList(counts = counts.keep)
###### v <- voom(dge, design)
# Ajustem el model lineal
###### fit <- lmFit(v, design)




########## 7. ESTIMACIÓ DEL MODEL
# No podem treballar dades log-transformades 
# amb logcounts tenim error per dades negatives
# La funció voom ajusta els comptatges
voom_data <- voom(counts.keep, design)
# Ajustem el model de limma per a dades aparellades (block per pacient)
block <- pheno_data$pacient 
dupcor <- duplicateCorrelation(voom_data, design, block=block)
# Ajustem el model tenint en compte la correlació 
fit <- lmFit(voom_data, design, block=block, correlation=dupcor$consensus)
# Definim els contrastos d'interès
contrast <- makeContrasts(Baseline_vs_PD = PD - Baseline, levels = design)
# Apliquem els contrastos a l'objecte 'fit'
fit2 <- contrasts.fit(fit, contrast)
# Estimació Bayesiana empírica per ajustar per variància
fit2 <- eBayes(fit2)
# Resultats per significació per p valor
TopTab <- topTable(fit2, adjust = "fdr", sort.by= "p",number = Inf)
print(head(TopTab))


########## 8. ANOTACIÓ DELS GENS
# Obtenim els Gene Symbols per als Entrez IDs seleccionats 
anotacions_gens <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(TopTab), columns = c("SYMBOL", "ENTREZID"), keytype = "ENTREZID")

########## 9. VISUALITZACIÓ ENTRE LES CONDICIONS
# Fem un volcano plot
volcano <- ggplot(TopTab, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4) +
  theme_classic() +
  labs(title = "Volcano Plot Baseline vs PD", x = "Log Fold Change", y = "-Log10 P-value") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed")

print(volcano)

# Heatmap expressió pels gens més significatius
# No som tan restrictius (P.value < 0.1 per obtenir més gens)
top_gens <- rownames(subset(TopTab, (abs(logFC)> 1.0) & (P.Value < 0.05)))
length(top_gens)
# Només ens interessen les dades d'expressió per als gens seleccionats
expression_data <- logcounts[top_gens, ]
# Centrem les dades per fer el heatmap
mat <- expression_data - rowMeans(expression_data)

# Ordenem les columnes per grups (baseline i PD)
noms_pacients <- colnames(mat)
baseline_samples <- noms_pacients[grepl("Baseline", noms_pacients)]
pd_samples <- noms_pacients[grepl("PD", noms_pacients)]
samples_ordenats <- c(baseline_samples, pd_samples)
mat <- mat[, samples_ordenats]

# Afegim informació fenotípica (condicions) com a anotacions de columna
annotation_col <- data.frame(Condition = ifelse(grepl("Baseline", colnames(mat)), "Baseline", "PD"))
rownames(annotation_col) <- colnames(mat)
# Per últim, afegim el nom en format SYMBOL
rownames(mat) <- anotacions_gens$SYMBOL[match(rownames(mat), anotacions_gens$ENTREZID)]
# Definim colors per a les condicions
colors_condicions_heatmap <- list(Condition = c(Baseline = "darkorchid1", PD = "darkgreen"))
# Colors del heatmap per incloure valors negatius
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(14)
# Creem el heatmap
pheatmap(mat, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         color = heatmap_colors, 
         breaks = seq(min(mat), max(mat), length.out = 19), 
         annotation_col = annotation_col, 
         annotation_colors = colors_condicions_heatmap, 
         main = "Gens amb abs(logFC)> 1.0 i P.Value < 0.05")




####### 10. Anàlisi de significació biològica
# Fem un GSEA per veure en quines categories trobem diferències
geneList <- TopTab$logFC
names(geneList) <- rownames(TopTab)
# Ordenem perlogFC
geneList <- sort(geneList, decreasing = TRUE)
# GSEA amb GO per processos biològics
gsea_results <- gseGO(geneList = geneList, 
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP",  
                      pvalueCutoff = 0.05, 
                      verbose = TRUE)
class(gsea_results)
# Per veure els resultats farem un dotplot
dotplot(gsea_results, showCategory = 10)
cnetplot(gsea_results,
         showCategory = 9, 
         node_label='category', 
         cex_gene = 0.1,
         cex_label_category = 2,
         circular = FALSE,
         colorEdge = TRUE)




######### 11.	Definició dels gens candidats 
gens_UP <- rownames(subset(TopTab, (logFC > 1) & (P.Value < 0.05)))
length(gens_UP)
# Obtenim els Gene Symbols per als Entrez IDs seleccionats 
anotacions_gens_UP<- AnnotationDbi::select(org.Hs.eg.db, keys = gens_UP, columns = c("SYMBOL", "ENTREZID"), keytype = "ENTREZID")
# Símbols dels gens seleccionats 
gens_UP_symbol <- anotacions_gens_UP$SYMBOL
gens_UP_symbol



#si som menys restrictius
gens_UP_0.1 <- rownames(subset(TopTab, (logFC > 0.7) & (P.Value < 0.1)))
length(gens_UP_0.1)
# Obtenim els Gene Symbols per als Entrez IDs seleccionats 
anotacions_gens_UP_0.1<- AnnotationDbi::select(org.Hs.eg.db, keys = gens_UP_0.1, columns = c("SYMBOL", "ENTREZID"), keytype = "ENTREZID")
# Símbols dels gens seleccionats 
gens_UP_symbol_0.1 <- anotacions_gens_UP_0.1$SYMBOL
gens_UP_symbol_0.1


# Dades de gens BT-474 i ZR75 CRISPR/CAS9 SCREEN CANDIDATES
# confidencial
gens_bt474 <- c("x","x","x","...","x")
# confidencial
gens_zr75 <- c("x","x","x","...","x")
  
gens_invitro <- union(gens_zr75, gens_bt474)
gens_invitro


gens_interseccio_rnaseq<-intersect(gens_UP_symbol_0.1,gens_invitro)
gens_interseccio_rnaseq

