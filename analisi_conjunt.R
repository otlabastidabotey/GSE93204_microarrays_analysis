setwd("C:/Users/Ot/Desktop/UOC/tfm")


####################################################
###################################### GENS IN VITRO

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

####################################################
######################################## GENS RNASEQ



# Llindar a 0,05
gens_interes_rnaseq_0.05 <- c("TNFRSF11B", "FGB", "FGG", "SERPINA1", "SPP1", 
                    "CP", "HP", "CRP", "FGA", "ITIH2", "FGL1")

# Llindar a 0,1

gens_interes_rnaseq_0.1 <- c("TNFRSF11B", "CDKN3", "CTSL", "UBE2T", "HMGB3", 
                        "CKS2", "CDC6", "TK1", "FGB", "FGG", "SERPINA1", 
                        "SPP1", "CP", "MB", "EFEMP1", "HP", "CRP", "FGA", 
                        "SPTSSB", "AADAC", "ITIH2", "NQO1", "EGLN3", "FGL1", 
                        "ORM2", "ORM1", "APOA2", "APCS", "SSX2", "GC", 
                        "FABP1", "APOB", "ALB", "VTN", "C19orf33", "UGT1A1", 
                        "ITIH3", "SLC2A2", "TM4SF4", "C4BPA", "HSD11B1", 
                        "LBP", "AMBP", "APOC3")


# comparacions
interseccio_rnaseq_0.05 <-intersect(gens_interes_rnaseq_0.05,gens_invitro)
interseccio_rnaseq_0.05

interseccio_rnaseq_0.1 <-intersect(gens_interes_rnaseq_0.1,gens_invitro)
interseccio_rnaseq_0.1




####################################################
################################### GENS MICROARRAYS

# Llindar a 0,05

gens_interes_micro_0.05 <- c("FOSB", "DUSP1", "NR4A3", "ZFP36", 
                  "EGR1", "IL6", "RASD1", "ATF3", 
                  "SERPINB9", "RGS1", "CCL2", "PDE4D", 
                  "PTGS2", "KLRC2", "POU2AF1", "KLRC1", 
                  "KLHDC7B", "LOC642838", "PHLDA1", "FOS", 
                  "TRAT1", "MALAT1", "LOC100653245", "LTF", 
                  "HPR", "SAA2", "SHISA7", "LINC00239", 
                  "PDK4", "SLC9B2", "FAM169A", "HIST1H1B", 
                  "ASPM", "DUSP2", "PTTG2", "CDCA2", 
                  "CASC5", "KIAA0101", "BUB1B", "DEPDC1", 
                  "PTTG1", "BIRC5", "CENPW", "STAR", 
                  "CCNB2", "LOC100507307", "FBXL5", "CDC45", 
                  "TYMS", "NELL2", "FANCI", "MND1", 
                  "JUN", "FAM72D", "CXCR4", "TOP2A", 
                  "CDC25C", "PBK", "MS4A1", "CFP", 
                  "TMEM38B", "ZWINT", "RRM2", "FOXM1", 
                  "E2F2", "BUB1", "HIST2H2AB", "PTPRZ1", 
                  "KIF23", "CD8B", "MCM10", "CXorf65", 
                  "SGOL1", "CD38", "TPX2", "POLQ", 
                  "SOCS3", "HIST1H1D", "TTK", "FAM54A", 
                  "ESCO2", "NCR3", "IRF4", "CCR2", 
                  "KIF18A", "CIT", "ZNF90", "SPIB", 
                  "HIST1H3H", "NDC80", "TNFRSF17", "CDCA7", 
                  "HIST1H3D", "CDC6", "GPR183", "PNOC", 
                  "DIAPH3", "LAMP3", "CD2", "CR2")


# Llindar a 0,1
gens_interes_micro_0.1 <- c("FOSB", "DUSP1", "NR4A3", "ZFP36", 
                      "EGR1", "IL6", "RASD1", "ATF3", 
                      "SERPINB9", "RGS1", "CCL2", "PDE4D", 
                      "PTGS2", "KLRC2", "POU2AF1", "KLRC1", 
                      "KLHDC7B", "LOC642838", "PHLDA1", "FOS", 
                      "TRAT1", "MALAT1", "LOC100653245", "LTF", 
                      "HPR", "SAA2", "SHISA7", "LINC00239", 
                      "PDK4", "SLC9B2", "FAM169A", "MS4A1", 
                      "CXorf65", "GPR183", "AIM2", "CYP19A1", 
                      "C4orf19", "SPIB", "DGKG", "LOC100507043", 
                      "TSPAN8", "CD38", "HSD11B1", "MYBL1", 
                      "ZNF483", "HIST1H1B", "ASPM", "DUSP2", 
                      "PTTG2", "CDCA2", "CASC5", "KIAA0101", 
                      "BUB1B", "DEPDC1", "PTTG1", "BIRC5", 
                      "CENPW", "STAR", "CCNB2", "LOC100507307", 
                      "FBXL5", "CDC45", "TYMS", "NELL2", 
                      "FANCI", "MND1", "JUN", "FAM72D", 
                      "CXCR4", "TOP2A", "CDC25C", "PBK", 
                      "CFP", "TMEM38B", "ZWINT", "RRM2", 
                      "FOXM1", "E2F2", "BUB1", "HIST2H2AB", 
                      "PTPRZ1", "KIF23", "CD8B", "MCM10", 
                      "SGOL1", "TPX2", "POLQ", "SOCS3", 
                      "HIST1H1D", "TTK", "FAM54A", "ESCO2", 
                      "NCR3", "IRF4", "CCR2", "KIF18A", 
                      "CIT", "ZNF90", "HIST1H3H", "NDC80", 
                      "TNFRSF17", "CDCA7", "HIST1H3D", "CDC6", 
                      "PNOC", "DIAPH3", "LAMP3", "CD2", 
                      "CR2", "CD1A", "ANLN", "SPTSSB", 
                      "SCGB1D2", "C8orf4", "HS6ST3", "C15orf48", 
                      "E2F8", "CXCL1", "FAM83D", "EFNA5")


# comparacions
interseccio_microarrays_0.05 <-intersect(gens_interes_micro_0.05,gens_invitro)
interseccio_microarrays_0.05

interseccio_microarrays_0.1 <-intersect(gens_interes_micro_0.1,gens_invitro)
interseccio_microarrays_0.1




