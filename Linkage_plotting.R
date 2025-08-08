library(openxlsx)
library(dplyr)
library(LinkageMapView)


folder<- "/Users/marina/Library/CloudStorage/OneDrive-Menntaský/QTL_mapping/Salpinus_align_LinkageMaps/H19-111/DeNovo/separatechromosome/genotypes_2"

fichiers <- list.files(path = folder, pattern = "^genotypes_\\d+\\.txt$", full.names = TRUE)

for (fichier in fichiers) {
  nom_fichier <- tools::file_path_sans_ext(basename(fichier))  # ex: genotypes_42
  assign(nom_fichier, read.delim(fichier, header = FALSE))
}

objets<- ls(pattern = "(^genotypes)")

for (gen in objets) {
  # Extraire le numéro depuis le nom (ex : 42 dans genotypes_42)
  num <- sub("genotypes_", "", gen)
  
  # Récupère l'objet
  df <- get(gen)
  
  # Vérifie qu’il a bien une colonne V2
  if (is.data.frame(df) && "V2" %in% names(df)) {
    # Remplace les "0" par le numéro
    df$V2[df$V2 == "0"] <- num
    
    # Réassigner le data.frame modifié
    assign(gen, df)
  }
}

for (gen in objets) {
  df<- get(gen)
  df<- df[-c(1:6),]
  assign(gen, df)
}

genotypes_list <- lapply(objets, function(obj) get(obj))

genotypes <- bind_rows(genotypes_list)

genotypes<- genotypes[ , -c(5:88)]

genotypes_male111<-genotypes[ , c(1:3)]
genotypes_female111<- genotypes[ , c(1,2,4)]

colnames(genotypes_male111)<- c("Marker", "Group", "Position")
colnames(genotypes_female111)<- c("Marker", "Group", "Position")

genotypes_male111<- genotypes_male111[ ,c(2,3,1)]
genotypes_female111<- genotypes_female111[ ,c(2,3,1)]

genotypes_female111$Position <- as.numeric(genotypes_female111$Position)
genotypes_male111$Position <- as.numeric(genotypes_male111$Position)

lmv.linkage.plot(mapthis = genotypes_male111, outfile = "Map_males_111.pdf")
lmv.linkage.plot(mapthis = genotypes_female111, outfile = "Map_females_111.pdf")

# Move male by female ####

data<- genotypes_female111

genotypes_female111$Position<- as.numeric(genotypes_female111$Position)

maxpos <- floor(max(genotypes_female111$Position))
at.axis <- seq(0, maxpos)

axlab <- vector()
for (lab in 0:maxpos) {
  if (!lab %% 10) {
    axlab <- c(axlab, lab)
  }
  else {
    axlab <- c(axlab, NA)
  }
}


outfile = file.path("denmap_female.pdf")
lmv.linkage.plot(genotypes_female111,outfile,denmap=TRUE, cex.axis = 1, at.axis = at.axis, labels.axis = axlab, main = "Male", pdf.height = 16)                         
