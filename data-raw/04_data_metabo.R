library(data.table)
library(openxlsx)
library(dplyr)
library(tidyr)
library(tibble)
library(janitor)

transpose_df <- function(df) {
    t_df <- data.table::transpose(df)
    colnames(t_df) <- rownames(df)
    rownames(t_df) <- colnames(df)
    t_df <- t_df %>%
        tibble::rownames_to_column(.data = .) %>%
        tibble::as_tibble(.)
    return(t_df)
}

# QC data ====
## FIA ====
FIA_QC <-  openxlsx::read.xlsx(
    xlsxFile = "data-raw/ADMC_DUKE_Q500_FIA_QC_20230522.xlsx",
    skipEmptyRows = TRUE, sheet = 'ADMC_DUKE_Q500_FIA_QC') |>
    tibble::tibble()

FIA_LODs <- FIA_QC |>
    dplyr::slice(1:8) |>
    dplyr::select_if( ~ !any(is.na(.))) |>
    transpose_df() |>
    janitor::row_to_names(row_number = 1) |>
    dplyr::rename("Short.Name/Injection.Number" = "Injection.Number")

FIA_QC_samples <-  FIA_QC |>
    dplyr::slice(9:n()) |>
    dplyr::select(-c(Customer.Sample.Identification,
                     Sample.Bar.Code, Sample.Type,
                     Sample.Identification, Species,
                     Material, Well.Position,
                     Sample.Volume, Run.Number,
                     Injection.Number))

K <- FIA_QC_samples |>
    process_data_with_report(LODs = FIA_LODs,
                             inj_type = "FIA",
                             qc_samples = TRUE)

FIA_QC_samples <- K$processed_data
K$report$removed_columns
K$report$imputed_columns
rm(K)

## UPLC ====
UPLC_QC <-  openxlsx::read.xlsx(
    xlsxFile = "data-raw/ADMC_DUKE_Q500_UPLC_QC_20230522.xlsx",
    skipEmptyRows = TRUE, sheet = 'ADMC_DUKE_Q500_UPLC_QC') |>
    tibble::tibble()

UPLC_LODs <- UPLC_QC |>
    dplyr::slice(1:5) |>
    dplyr::select_if(~ !any(is.na(.))) |>
    transpose_df() |>
    janitor::row_to_names(row_number = 1) |>
    dplyr::rename("Short.Name/Injection.Number" = "Injection.Number")

UPLC_QC_samples <- UPLC_QC |>
    dplyr::slice(6:n()) |>
    dplyr::select(-c(Customer.Sample.Identification,
                     Sample.Bar.Code, Sample.Type,
                     Sample.Identification, Species,
                     Material, Well.Position,
                     Sample.Volume, Run.Number,
                     Injection.Number))

K <- UPLC_QC_samples |>
    process_data_with_report(LODs = UPLC_LODs,
                             inj_type = "UPLC",
                             qc_samples = TRUE)

UPLC_QC_samples <- K$processed_data
K$report$removed_columns
K$report$imputed_columns
rm(K)
# Metabo data ====
## FIA ====
FIA_Metabo <- data.table::fread(
    input = "data-raw/ADMC_DUKE_Q500_FIA_11Nov2024.csv",
    sep = ",") |> tibble::tibble()

FIA_Metabo <- FIA_Metabo |>
    dplyr::select(-c(`Customer Sample Identification`,
                     `Sample Bar Code`, `Sample Type`,
                     `Sample Identification`, Species,
                     Material, `Well Position`,
                     `Sample Volume`, `Run Number`,
                     `Injection Number`,
                     EXAMDATE, `Plate Bar Code`, update_stamp))

FIA_processed <- FIA_Metabo |>
    dplyr::filter(VISCODE2 == "bl") |>
    ## add_metadata(other_info = adnimerge) |>
    process_data_with_report(LODs = FIA_LODs, inj_type = "FIA")

## UPLC ====
UPLC_Metabo <- data.table::fread(
    input = "data-raw/ADMC_DUKE_Q500_UPLC_11Nov2024.csv",
    sep = ",") |> tibble::tibble()

UPLC_Metabo <- UPLC_Metabo |>
    dplyr::select(-c(`Customer Sample Identification`,
                     `Sample Bar Code`, `Sample Type`,
                     `Sample Identification`, Species,
                     Material, `Well Position`,
                     `Sample Volume`, `Run Number`,
                     `Injection Number`,
                     EXAMDATE, `Plate Bar Code`, update_stamp))

UPLC_processed <- UPLC_Metabo |>
    dplyr::filter(VISCODE2 == "bl") |>
    ## add_metadata(other_info = adnimerge) |>
    process_data_with_report(LODs = UPLC_LODs, inj_type = "UPLC")


# Merge FIA and UPLC with their corresponding QCs ====

All_Metabo <- UPLC_processed$processed_data |> dplyr::full_join(
    FIA_processed$processed_data, by = c("RID" = "RID", "VISCODE2" = "VISCODE2"))


All_Metabo |> dplyr::filter(VISCODE2 == "bl") |>
    dplyr::select(-c(VISCODE2)) |>
    add_metadata(other_info = adnimerge)


# PCA for FIA raw ====

# replace space by point in metabolite names
colnames(FIA_Metabo) <-
    gsub(" ", ".", FIA_Metabo |> colnames())

raw_data <- FIA_processed$processed_data |>
    cbind(UPLC_processed$processed_data[-c(1:9)]) |>
    data.frame()

raw_data1 <- raw_data[-c(1:9)] |>
    dplyr::mutate_at(colnames(raw_data[-c(1:9)]), as.numeric)
raw_data1 <- raw_data1[sapply(raw_data1, is.numeric)]

dplyr::bind_rows(raw_data1, cbind(FIA_QC_samples |> dplyr::select(-1),
                         UPLC_QC_samples |> dplyr::select(-1)))

mypca <-  prcomp(raw_data1)
pca <- cbind(raw_data[1:9], mypca$x[], 1:2)
ggplot(pca, aes(PC1,PC2, col = DX_bl, fill = DX_bl)) +
    stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
    geom_point(shape = 21, col = "black")


biomarkers <- raw_data[1:9]
acp <- PCA(raw_data1, graph = FALSE)
pca_coords <- as.data.frame(acp$ind$coord[, 1:5])
colnames(pca_coords) <- paste0("Dim", 1:5)

pca_coords$Cohort <- as.factor(raw_data$DX_bl)
pca_coords$sampleID <- raw_data$RID

coul <- rainbow(length(unique(pca_coords$Cohort)))
color_map <- setNames(coul, levels(pca_coords$Cohort))
pca_coords$Color <- color_map[pca_coords$Cohort]

axis = list(showline = FALSE, zeroline = FALSE, gridcolor = '#ffff',
            ticklen = 4, titlefont = list(size = 13))


library(plotly)
splom_plot <- plot_ly(
    data = pca_coords,
    type = 'splom',
    dimensions = list(
        list(label = paste("PC1: (", round(acp$eig[1, 2], 2), "%)", sep = ""), values = ~Dim1),
        list(label = paste("PC2: (", round(acp$eig[2, 2], 2), "%)", sep = ""), values = ~Dim2),
        list(label = paste("PC3: (", round(acp$eig[3, 2], 2), "%)", sep = ""), values = ~Dim3),
        list(label = paste("PC4: (", round(acp$eig[4, 2], 2), "%)", sep = ""), values = ~Dim4)
    ),
    color = ~pca_coords$Cohort, colors = c('#636EFA','#EF553B'),
    text = ~paste("sample: ", sampleID)
) %>%
    layout(
        title = "Scatterplot Matrix of PCA 4 Dimensions",
        dragmode = "select",
        legend = list(title = list(text = 'Group')),
        hovermode = 'closest',
        dragmode = 'select',
        plot_bgcolor = 'rgba(240,240,240, 0.95)',
        xaxis = list(domain = NULL, showline = F, zeroline = F, gridcolor = '#ffff', ticklen = 4),
        yaxis = list(domain = NULL, showline = F, zeroline = F, gridcolor = '#ffff', ticklen = 4),
        xaxis2 = axis,
        xaxis3 = axis,
        xaxis4 = axis,
        yaxis2 = axis,
        yaxis3 = axis,
        yaxis4 = axis
    )  %>% style(diagonal = list(visible = FALSE))

splom_plot

rm("FIA_QC", "UPLC_QC")

# Save the data sets ====
usethis::use_data(FIA_LODs, overwrite = TRUE)
usethis::use_data(FIA_QC_samples, overwrite = TRUE)
usethis::use_data(FIA_Metabo, overwrite = TRUE)
usethis::use_data(UPLC_LODs, overwrite = TRUE)
usethis::use_data(UPLC_QC_samples, overwrite = TRUE)
usethis::use_data(UPLC_Metabo, overwrite = TRUE)





