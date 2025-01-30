
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
BiocManager::install("GO.db")
if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
BiocManager::install("impute")
if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
BiocManager::install("preprocessCore")
if (!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA")
if (!requireNamespace("gplots", quietly = TRUE)) install.packages("gplots")
if (!requireNamespace("cluster", quietly = TRUE)) install.packages("cluster")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("biomaRt")

library(biomaRt)
library(clusterProfiler)  # Analiza wzbogacania
library(org.Hs.eg.db)     # Baza danych dla genów ludzkich
library(enrichplot)       # Wizualizacje dla wyników wzbogacenia
library(pheatmap)
library(cluster)
library(ggplot2)
library(readxl)
library(dplyr)
library(WGCNA)
library(gplots)
library(WGCNA)
library(ggplot2)
library(readxl)
library(dplyr)


options(stringsAsFactors = FALSE)

# Ładowanie danych
file_path <- "C:/Users/mczarkow/Downloads/summary_htseq_norm.xlsx"
summary_htseq <- read_excel(file_path, sheet = "summary_htseq")
gene_column <- "Name"
available_samples <- names(summary_htseq)[-(1:5)]  # Pomiń pierwsze 5 kolumn

num_lncRNA_genes <- sum(summary_htseq$lncRNA == 1, na.rm = TRUE)
cat("Liczba genów z lncRNA = 1:", num_lncRNA_genes, "\n")

# Filtrowanie danych
filtered_genes <- summary_htseq[summary_htseq$var > 10, ]
gene_expression_matrix <- as.data.frame(filtered_genes[, available_samples])
rownames(gene_expression_matrix) <- filtered_genes[[gene_column]]
if (!(gene_column %in% colnames(filtered_genes))) {
  stop("Kolumna z nazwami genów nie istnieje! Sprawdź nazwę kolumny.")
}

# Przypisanie nazw genów do wierszy
gene_names <- filtered_genes[[gene_column]]  # Pobierz nazwy genów z kolumny
if (anyDuplicated(gene_names)) {
  gene_names <- make.unique(gene_names)  # Sprawdź, czy nazwy genów są unikalne
}
datExpr <- as.data.frame(t(gene_expression_matrix))

# Sprawdzenie jakości genów i próbek
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Dobór optymalnej wartości soft-thresholding power
powers <- c(1:30)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Ustawienie nazwy pliku wyjściowego
output_file <- "soft_thresholding_analysis.png"

# Tworzenie pliku PNG
png(filename = output_file, width = 1800, height = 900, res = 150)

# Ustawienie układu graficznego na dwa wykresy obok siebie
par(mfrow = c(1,2))

# Sprawdzenie, czy dane istnieją
if (!is.null(sft$fitIndices)) {
  
  # Wykres zależności R^2 od power (dostosowany)
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], 
       type = "b", 
       col = "blue",
       pch = 19, 
       xlab = "Soft Threshold (power)", 
       ylab = "Scale-Free Topology Model Fit, signed R²", 
       main = "Dobór optymalnej wartości soft-thresholding",
       cex.main = 1, cex.lab = 0.9, cex.axis = 0.8)
  
  # Dodanie poziomej linii referencyjnej dla R² = 0.9
  abline(h=0.9, col="red", lty=2)
  
  # Wykres średniego stopnia łączności w zależności od power
  plot(sft$fitIndices[,1], sft$fitIndices[,5], 
       type = "b", 
       col = "blue",
       pch = 19, 
       xlab = "Soft Threshold (power)", 
       ylab = "Mean Connectivity", 
       main = "Średnia łączność w zależności od power",
       cex.main = 1, cex.lab = 0.9, cex.axis = 0.8)
} else {
  plot.new()
  text(0.5, 0.5, "Brak danych do wygenerowania wykresu", cex = 1.5)
}

# Zamknięcie pliku PNG
dev.off()

# Sprawdzenie, czy plik został zapisany poprawnie
if (file.exists(output_file)) {
  cat("Wykres został zapisany jako:", output_file, "\n")
} else {
  cat("Błąd: Plik nie został zapisany. Sprawdź uprawnienia i katalog roboczy.\n")
}


# Wybór optymalnej wartości softPower
softPower <- sft$powerEstimate
cat("Optymalna wartość soft-thresholding power:", softPower, "\n")

##### Analiza sieci w blokach (wieloblokowa analiza) - WGCNA ######
maxBlockSize <- 556  # Maksymalny rozmiar bloku
bwnet <- blockwiseModules(
  datExpr,
  maxBlockSize = maxBlockSize,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 10,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "TOM",
  verbose = 3
)

# Wyodrębnienie kolorów modułów i eigengenów
moduleColors_blockwise <- bwnet$colors
MEs_blockwise <- bwnet$MEs
save(bwnet, moduleColors_blockwise, MEs_blockwise, file = "blockwiseModules_results.RData")

# Wykresy dendrogramów dla każdego bloku
plotDendrogramForModules <- function(datExpr, moduleColors) {
  unique_modules <- unique(moduleColors)
  
  for (module in unique_modules) {
    cat("Generowanie dendrogramu dla modułu:", module, "\n")
    
    genes_in_module <- names(moduleColors)[moduleColors == module]
    
    if (length(genes_in_module) < 2) {
      cat("Za mało genów do dendrogramu dla modułu", module, "\n")
      next
    }
    
    module_data <- datExpr[, genes_in_module, drop = FALSE]
    
    if (nrow(module_data) < 2 || ncol(module_data) < 2) {
      cat("Za mało danych w module", module, "do utworzenia dendrogramu\n")
      next
    }
    
    tryCatch({
      module_data <- t(module_data)
      dist_matrix <- as.dist(1 - cor(t(module_data), use = "pairwise.complete.obs"))
      clustering <- hclust(dist_matrix, method = "average")
      
      output_filename <- paste0("dendrogram_module_", sprintf("%02d", module), ".png")
      png(filename = output_filename, width = 800, height = 600)
      
      plot(
        clustering,
        main = paste("Dendrogram dla modułu", module),
        xlab = "Geny",
        sub = "",
        hang = 0.03
      )
      
      # Dynamically set `k`, ensuring it's at least 2
      k_value <- min(25, length(genes_in_module))  # Minimum 2, max 5 (możesz dostosować)
      
      rect.hclust(
        clustering,
        k = k_value,  
        border = moduleColors[genes_in_module][1]  
      )
      
      dev.off()
      
      cat("Dendrogram dla modułu", module, "zapisany do pliku:", output_filename, "\n")
      
    }, error = function(e) {
      dev.off()
      cat("Błąd podczas generowania dendrogramu dla modułu", module, ":", e$message, "\n")
    })
  }
  
  cat("Dendrogramy dla modułów zostały wygenerowane.\n")
}


# Wywołanie funkcji
plotDendrogramForModules(datExpr, moduleColors_blockwise)

# Extract gene names from the original data (column names of datExpr)
geneNames <- colnames(datExpr)  # Genes are in columns of datExpr

# Ensure the number of genes matches the length of moduleColors
if (length(moduleColors_blockwise) != length(geneNames)) {
  stop("Mismatch between moduleColors and the number of genes.")
}

# Create a data frame with genes and their assigned modules
modules_and_genes <- data.frame(
  Gene = geneNames,
  Module = moduleColors_blockwise
)

# Save the combined data to a CSV file
output_file <- "Modules_and_Genes.csv"
write.csv(modules_and_genes, file = output_file, row.names = FALSE)

cat("Klastry i geny zapisane do pliku:", output_file, "\n")

# Konwersja do ramki danych
module_sizes_df <- as.data.frame(module_sizes)

# Nazwanie kolumn (opcjonalne, dla przejrzystości)
colnames(module_sizes_df) <- c("Moduł", "Liczba_genów")

# Ustawienie nazwy pliku
output_file <- "module_sizes_plot.png"

# Zapis wykresu do pliku PNG
png(filename = output_file, width = 1000, height = 800, res = 150)

ggplot(data = module_sizes_df, aes(x = Moduł, y = Liczba_genów)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(
    title = "Liczba genów w modułach",
    x = "Moduł",
    y = "Liczba genów"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

cat("Wykres został zapisany jako:", output_file, "\n")

#Heatmapy
##HEATMAPY
generateHeatmapForModules <- function(datExpr, moduleColors, module) {
  # Wybierz geny przypisane do danego modułu
  genes_in_module <- names(moduleColors)[moduleColors == module]
  
  # Sprawdź, czy są geny w module
  if (length(genes_in_module) < 2) {
    cat("Za mało genów do heatmapy dla modułu", module, "\n")
    return(NULL)
  }
  
  # Sprawdź, czy wszystkie geny z modułu istnieją w danych
  missing_genes <- setdiff(genes_in_module, colnames(datExpr))
  if (length(missing_genes) > 0) {
    cat("Niektóre geny w module", module, "nie istnieją w danych ekspresji:\n")
    print(missing_genes)
    return(NULL)
  }
  
  # Ekstrakcja danych ekspresji dla genów w module
  module_data <- datExpr[, genes_in_module, drop = FALSE]
  
  # Sprawdź, czy dane mają wystarczający rozmiar
  if (nrow(module_data) < 2 || ncol(module_data) < 2) {
    cat("Za mało danych (geny lub próbki) do heatmapy dla modułu", module, "\n")
    return(NULL)
  }
  
  # Transpozycja danych, aby geny były w wierszach, a próbki w kolumnach
  module_data <- t(module_data)
  
  # Macierz odległości
  row_dist <- as.dist(1 - cor(t(module_data))) # Dla genów
  col_dist <- as.dist(1 - cor(module_data))   # Dla próbek
  
  # Klasteryzacja
  row_clustering <- hclust(row_dist, method = "average")  # Geny
  col_clustering <- hclust(col_dist, method = "average") # Próbki
  
  # Posortowanie danych na podstawie klastrowania
  module_data <- module_data[row_clustering$order, col_clustering$order]
  
  # Tworzenie heatmapy i zapisywanie do pliku
  output_filename <- paste0("heatmap_module_", module, ".png")
  png(filename = output_filename, width = 800, height = 600)
  pheatmap(
    module_data,
    cluster_rows = FALSE,  # Geny są już posortowane
    cluster_cols = FALSE,  # Próby są już posortowane
    main = paste("Heatmapa modułu", module, "\n(Geny: Wiersze, Próbki: Kolumny)"),
    color = colorRampPalette(c("blue", "white", "red"))(50),
    fontsize = 8
  )
  dev.off()
  
  # Wyświetlenie heatmapy w konsoli
  pheatmap(
    module_data,
    cluster_rows = FALSE,  # Geny są już posortowane
    cluster_cols = FALSE,  # Próby są już posortowane
    main = paste("Heatmapa modułu", module, "\n(Geny: Wiersze, Próbki: Kolumny)"),
    color = colorRampPalette(c("blue", "white", "red"))(50),
    fontsize = 8
  )
  
  cat("Heatmapa modułu", module, "zapisana do pliku:", output_filename, "\n")
}

# Generowanie heatmapy dla wszystkich modułów
unique_modules <- unique(moduleColors_blockwise)
for (module in unique_modules) {
  cat("\nGenerowanie heatmapy dla modułu:", module, "\n")
  generateHeatmapForModules(datExpr, moduleColors_blockwise, module)
}


###### Silhouette dla każdego modułu ######
# Transpozycja macierzy datExpr
# Transpozycja macierzy datExpr
datExpr_genes <- t(datExpr)


# Dopasowanie genów między datExpr_genes a modules_and_genes
common_genes <- intersect(rownames(datExpr_genes), modules_and_genes$Gene)

# Znajdź geny, które nie pasują
missing_in_datExpr <- setdiff(modules_and_genes$Gene, rownames(datExpr_genes))
missing_in_modules <- setdiff(rownames(datExpr_genes), modules_and_genes$Gene)

# Wyświetl problematyczne geny
cat("Geny w modules_and_genes, ale nie w datExpr_genes:\n")
print(missing_in_datExpr)

cat("Geny w datExpr_genes, ale nie w modules_and_genes:\n")
print(missing_in_modules)


# Filtrowanie i dopasowanie kolejności
datExpr_genes <- datExpr_genes[common_genes, ]
modules_and_genes <- modules_and_genes[modules_and_genes$Gene %in% common_genes, ]
modules_and_genes <- modules_and_genes[match(rownames(datExpr_genes), modules_and_genes$Gene), ]

# Konwersja modułów na liczby całkowite
gene_clusters <- as.integer(modules_and_genes$Module)

# Sprawdzenie długości obu struktur
if (length(gene_clusters) != nrow(datExpr_genes)) {
  stop("Liczba genów w gene_clusters i datExpr_genes nie jest zgodna!")
}

# Macierz odległości na podstawie korelacji między genami
cor_matrix <- cor(t(datExpr_genes))  # Transpozycja, aby geny były w kolumnach
dist_matrix <- as.dist(1 - cor_matrix)


# Obliczanie silhouette score dla genów
silhouette_scores <- silhouette(gene_clusters, dist_matrix)

# Wizualizacja wyników silhouette
plot(silhouette_scores, col = as.numeric(modules_and_genes$Module), border = NA)

# Wyciąganie średnich wartości silhouette dla każdego modułu
silhouette_avg <- aggregate(silhouette_scores[, "sil_width"], 
                            by = list(Module = modules_and_genes$Module), 
                            FUN = mean)
colnames(silhouette_avg) <- c("Module", "Average_Silhouette_Score")

# Wyświetlenie średnich wartości silhouette
print(silhouette_avg)

# Zapis wyników do pliku CSV
output_silhouette_file <- "Silhouette_Scores_Genes.csv"
write.csv(silhouette_avg, file = output_silhouette_file, row.names = FALSE)

cat("Średnie wartości silhouette zapisane do pliku:", output_silhouette_file, "\\n")

print(dim(datExpr_genes))
print(length(gene_clusters))


#Sprawdzanie klastrów
# Utworzenie ramki danych z genami i przypisaniem do modułów
modules_and_genes <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors_blockwise
)

# Zapisanie ramki danych do pliku CSV
output_file <- "modules_and_genes.csv"
write.csv(modules_and_genes, file = output_file, row.names = FALSE)

cat("Plik z genami i modułami został zapisany jako:", output_file, "\n")

# Wyświetlenie zawartości każdego modułu
for (module in unique(modules_and_genes$Module)) {
  cat("\nModuł:", module, "\n")
  print(modules_and_genes[modules_and_genes$Module == module, ])
}

# Zlicz geny w każdym module
module_sizes <- table(moduleColors_blockwise)


# Znajdź małe moduły
small_modules <- module_sizes[module_sizes < 30]
cat("Moduły z mniej niż 30 genami:\n")
print(small_modules)

#############################################WZBOGACENIE
# Create a data frame with genes and their assigned modules
modules_and_genes <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors_blockwise
)

write.csv(modules_and_genes, "Modules_And_Genes.csv", row.names = FALSE)



module_gene_list <- split(modules_and_genes$Gene, modules_and_genes$Module)

# Lista na zbiorcze wyniki
all_results <- list()
go_enrichment_results <- list()  # Zdefiniowanie listy na wyniki analizy wzbogacenia

total_genes_mapped <- 0  # Liczba wszystkich zmapowanych genów

# Przeprowadzenie analizy wzbogacenia dla każdego modułu
for (module in names(module_gene_list)) {
  cat("\nAnaliza wzbogacenia dla modułu:", module, "\n")
  
  genes <- module_gene_list[[module]]
  
  # Wyświetlenie liczby genów przed wzbogaceniem
  cat("Liczba genów w module:", length(genes), "\n")
  
  # Pomijanie modułów z mniej niż 2 genami
  if (length(genes) < 2) {
    cat("Za mało genów do analizy dla modułu:", module, "\n")
    next
  }
  
  # Konwersja identyfikatorów genów na ENTREZID
  gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Logowanie genów, które się nie mapują
  unmapped_genes <- setdiff(genes, gene_entrez$SYMBOL)
  if (length(unmapped_genes) > 0) {
    cat("Nie udało się zmapować genów dla modułu", module, ":\n")
    print(unmapped_genes)
  }
  
  # Wyświetlenie liczby zmapowanych genów
  if (!is.null(gene_entrez)) {
    num_genes_mapped <- nrow(gene_entrez)
    total_genes_mapped <- total_genes_mapped + num_genes_mapped
    cat("Liczba genów zmapowanych na ENTREZID:", num_genes_mapped, "\n")
  } else {
    cat("Brak genów zmapowanych na ENTREZID dla modułu:", module, "\n")
    next
  }
  
  # Sprawdzenie, czy `bitr` zwrócił jakiekolwiek wyniki
  if (nrow(gene_entrez) == 0) {
    cat("Brak genów do analizy dla modułu:", module, "\n")
    next
  }
  
  # Przeprowadzenie analizy wzbogacenia
  go_results <- enrichGO(
    gene          = gene_entrez$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",  # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  # Sprawdzanie wyników i zapisywanie
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    result_df <- as.data.frame(go_results)
    result_df$Module <- module
    
    # Zapisanie wyników do pliku CSV
    write.csv(result_df, paste0("GO_Enrichment_Module_", module, ".csv"), row.names = FALSE)
    all_results[[module]] <- result_df
  } else {
    cat("Brak istotnych wyników dla modułu:", module, "\n")
  }
}

# Scal wszystkie wyniki w jedną ramkę danych
if (length(all_results) > 0) {
  combined_results <- do.call(rbind, all_results)
  
  # Zapisz zbiorcze wyniki do pliku CSV
  write.csv(combined_results, "GO_Enrichment_Combined_Results_WGCNA.csv", row.names = FALSE)
  cat("Zbiorcze wyniki zapisane w pliku: GO_Enrichment_Combined_Results.csv\n")
} else {
  cat("Brak wyników do scalenia.\n")
}

# Wyświetlenie zbiorczej liczby genów zmapowanych
cat("\nZbiorcza liczba genów zmapowanych na ENTREZID we wszystkich modułach:", total_genes_mapped, "\n")


################################## WZBOGCENIE PYTHON
  ################################## WZBOGCENIE PYTHON
  # Wczytanie danych z pliku CSV
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  # Wczytanie danych z pliku CSV wygenerowanego w Pythonie
  file_path <- "C:/Users/mczarkow/Downloads/kmeans_40_clusters-20(1).csv"
  modules_and_genes_python <- read.csv(file_path)
  
  # Przygotowanie listy genów dla każdego modułu
  module_gene_list_python <- split(modules_and_genes_python$Name, modules_and_genes_python$Module)
  
  # Lista na zbiorcze wyniki
  all_results <- list()
  go_enrichment_results <- list()  # Zdefiniowanie listy na wyniki analizy wzbogacenia
  
  # Zmienna do śledzenia liczby genów zmapowanych do ENTREZID
  total_genes_mapped <- 0
  
  # Przeprowadzenie analizy wzbogacenia dla każdego modułu
  for (module in names(module_gene_list_python)) {
    cat("\nAnaliza wzbogacenia dla modułu:", module, "\n")
    
    genes <- module_gene_list_python[[module]]
    
    # Wyświetlenie liczby genów przed wzbogaceniem
    cat("Liczba genów w module:", length(genes), "\n")
    
    # Pomijanie modułów z mniej niż 2 genami
    if (length(genes) < 2) {
      cat("Za mało genów do analizy dla modułu", module, "\n")
      next
    }
    
    # Konwersja identyfikatorów genów na ENTREZID
    gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    # Logowanie genów, które się nie mapują
    unmapped_genes <- setdiff(genes, gene_entrez$SYMBOL)
    if (length(unmapped_genes) > 0) {
      cat("Nie udało się zmapować genów dla modułu", module, ":\n")
      print(unmapped_genes)
    }
    
    # Wyświetlenie liczby zmapowanych genów
    if (!is.null(gene_entrez)) {
      num_genes_mapped <- nrow(gene_entrez)
      total_genes_mapped <- total_genes_mapped + num_genes_mapped
      cat("Liczba genów zmapowanych na ENTREZID:", num_genes_mapped, "\n")
    } else {
      cat("Brak genów zmapowanych na ENTREZID dla modułu:", module, "\n")
      next
    }
    
    # Sprawdzenie, czy `bitr` zwrócił jakiekolwiek wyniki
    if (nrow(gene_entrez) == 0) {
      cat("Brak genów do analizy dla modułu", module, "\n")
      next
    }
    
    # Przeprowadzenie analizy wzbogacenia
    go_results <- enrichGO(
      gene          = gene_entrez$ENTREZID,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    
    # Sprawdzanie wyników i zapisywanie
    if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
      result_df <- as.data.frame(go_results)
      result_df$Module <- module
      write.csv(result_df, paste0("GO_Enrichment_Module_", module, ".csv"), row.names = FALSE)
      all_results[[module]] <- result_df
    } else {
      cat("Brak istotnych wyników dla modułu", module, "\n")
    }
  }
  
  # Scal wszystkie wyniki w jedną ramkę danych
  if (length(all_results) > 0) {
    combined_results <- do.call(rbind, all_results)
    
    # Zapisz zbiorcze wyniki do pliku CSV
    write.csv(combined_results, "GO_Enrichment_Combined_Results_Python.csv", row.names = FALSE)
    cat("Zbiorcze wyniki zapisane w pliku: GO_Enrichment_Combined_Results.csv\n")
  } else {
    cat("Brak wyników do scalenia.\n")
  }
  
  # Wyświetlenie zbiorczej liczby genów zmapowanych
  cat("\nZbiorcza liczba genów zmapowanych na ENTREZID we wszystkich modułach:", total_genes_mapped, "\n")
  
  #################################### Zebranie wyników wzbogacenia i silhoutte między Python i WGCNA
  
  # Wczytanie wymaganych bibliotek
  if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
  library(openxlsx)
  library(dplyr)
  
  # Wczytanie wyników z plików CSV
  results_r <- read.csv("GO_Enrichment_Combined_Results_WGCNA.csv")
  results_python <- read.csv("GO_Enrichment_Combined_Results_Python.csv")
  
  # Wczytanie wyników silhouette z R i Pythona
  file_path <- "C:/Users/mczarkow/Downloads/kmeans_python.csv"
  silhouette_r <- read.csv("Silhouette_Scores_Genes.csv") %>%
    rename(Silhouette_R = Average_Silhouette_Score)
  
  silhouette_python <- read.csv(file_path) %>%
    rename(Silhouette_Python = Value)
  
  # Wczytanie danych o genach i modułach
  modules_and_genes <- read.csv("Modules_And_Genes.csv")  # Plik musi zawierać kolumny: Gene, Module
  
  # Zapisanie modules_and_genes do pliku CSV
  write.csv(modules_and_genes, "Modules_And_Genes_Output.csv", row.names = FALSE)
  
  # Połączenie wyników silhouette z modułami
  silhouette_r <- merge(modules_and_genes, silhouette_r, by = "Module", all.x = TRUE)
  silhouette_python <- merge(modules_and_genes_python, silhouette_python, by = "Module", all.x = TRUE)
  
  # Funkcja do tworzenia podsumowania w hierarchicznej strukturze
  create_hierarchical_summary_excel <- function(modules_and_genes, go_results, silhouette_scores, output_file) {
    summary_list <- list()
    summary_list[[1]] <- c("Klaster", "L_genów", "Miara_sil", "Term_GO", "pval", "GeneRatio", "BgRatio")
    
    sorted_modules <- sort(unique(modules_and_genes$Module))
    row_index <- 2  
    
    for (module in sorted_modules) {
      genes_in_module <- modules_and_genes[modules_and_genes$Module == module, "Gene"]
      n_genes <- length(genes_in_module)
      
      silhouette_value <- silhouette_scores[silhouette_scores$Module == module, ]$Silhouette_R
      if (length(silhouette_value) == 0) silhouette_value <- NA
      if (length(silhouette_value) > 1) silhouette_value <- mean(silhouette_value, na.rm = TRUE)
      
      go_terms <- go_results[go_results$Module == module, ]
      
      summary_list[[row_index]] <- c(paste("Klaster", module), n_genes, silhouette_value, "", "", "", "")
      row_index <- row_index + 1
      
      if (nrow(go_terms) > 0) {
        for (i in 1:nrow(go_terms)) {
          summary_list[[row_index]] <- c(
            "", "", "",
            go_terms$Description[i],
            go_terms$p.adjust[i],
            go_terms$GeneRatio[i],
            go_terms$BgRatio[i]
          )
          row_index <- row_index + 1
        }
      } else {
        summary_list[[row_index]] <- c("", "", "", "Brak wyników GO", "", "", "")
        row_index <- row_index + 1
      }
    }
    
    summary_df <- do.call(rbind, lapply(summary_list, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
    
    write.xlsx(summary_df, output_file, colNames = FALSE, rowNames = FALSE)
    cat("Hierarchiczne podsumowanie zapisano w pliku:", output_file, "\n")
  }
  
  # Tworzenie plików dla wyników WGCNA i Pythona
  create_hierarchical_summary_excel(
    modules_and_genes = modules_and_genes,
    go_results = results_r,
    silhouette_scores = silhouette_r,
    output_file = "Hierarchiczne_Podsumowanie_WGCNA.xlsx"
  )
  
  create_hierarchical_summary_excel(
    modules_and_genes = modules_and_genes_python,
    go_results = results_python,
    silhouette_scores = silhouette_python,
    output_file = "Hierarchiczne_Podsumowanie_Python.xlsx"
  )
  
  
  
  results_r <- read.csv("GO_Enrichment_Combined_Results_WGCNA.csv")
  
  # Podliczenie częstości występowania poszczególnych opisów (Description)
  top_descriptions <- results_r %>%
    group_by(Description) %>%
    summarise(Frequency = n()) %>%
    arrange(desc(Frequency)) %>%
    top_n(20, wt = Frequency)  # Pobranie 20 najczęstszych terminów GO
  
  # Wykres najczęstszych terminów GO
  ggplot(top_descriptions, aes(x = reorder(Description, Frequency), y = Frequency)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    coord_flip() +
    labs(title = "Top 20 najczęstszych terminów GO", x = "Opis GO", y = "Liczba wystąpień") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),  # Białe tło
      panel.background = element_rect(fill = "white", color = NA), # Białe tło panelu
      legend.background = element_rect(fill = "white", color = NA)
    )
  
  # Zapis wykresu do pliku
  ggsave("Top_20_GO_Descriptions_R.png", width = 10, height = 6, dpi = 300)
  
  
  
  
  # Wczytanie danych wyników wzbogacenia wygenerowanych w Pythonie
  results_python <- read.csv("GO_Enrichment_Combined_Results_Python.csv")
  
  # Sprawdzenie struktury danych
  str(results_python)
  
  # Podliczenie częstości występowania poszczególnych opisów (Description)
  library(dplyr)
  top_descriptions_python <- results_python %>%
    group_by(Description) %>%
    summarise(Frequency = n()) %>%
    arrange(desc(Frequency)) %>%
    top_n(20, wt = Frequency)  # Pobranie 20 najczęstszych terminów GO
  
  # Wykres najczęstszych terminów GO dla danych z Pythona
  library(ggplot2)
  ggplot(top_descriptions_python, aes(x = reorder(Description, Frequency), y = Frequency)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    coord_flip() +
    labs(title = "Top 20 najczęstszych terminów GO (Dane z Pythona)", 
         x = "Opis GO", 
         y = "Liczba wystąpień") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),  # Białe tło
      panel.background = element_rect(fill = "white", color = NA), # Białe tło panelu
      legend.background = element_rect(fill = "white", color = NA)
    )
  
  # Zapis wykresu do pliku PNG
  ggsave("Top_20_GO_Descriptions_From_Python.png", width = 10, height = 6, dpi = 300)
  