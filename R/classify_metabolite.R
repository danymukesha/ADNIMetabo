#' Classify Metabolites into Biocrates MxP® Quant 500 Classes
#'
#' This function classifies metabolite names into their corresponding
#' Biocrates MxP® Quant 500 metabolite classes based on regex pattern matching.
#' It covers up to 26 biochemical classes including Acylcarnitines, Ceramides,
#' Cholesteryl Esters, Diglycerides, Triglycerides, and many others.
#'
#' @param met A character string representing the metabolite name.
#'
#' @return A character string specifying the metabolite class. Returns
#' "Unknown" if no class matches.
#'
#' @examples
#' classify_metabolite("C5-M-DC") # "Acylcarnitines"
#' classify_metabolite("Cer(d18:1/24:1)") # "Ceramides"
#' classify_metabolite("TG(18:1_34:2)") # "Triglycerides"
#' classify_metabolite("UnknownMet") # "Unknown"
#'
#' \dontrun{
#' library(dplyr)
#' library(stringr)
#'
#' all_features <- ADNIMetabo::reformatted_data_all |>
#'     dplyr::select(-c("RID", "VISCODE", "Plate")) |>
#'     colnames()
#'
#' metabolite_table <- tibble(
#'     Metabolite = all_features,
#'     Class = sapply(all_features, classify_metabolite)
#' )
#'
#' # View unknown metabolites
#' metabolite_table |> filter(Class == "Unknown")
#'
#' # Save to CSV
#' write.csv(metabolite_table, "metabolite_class_table.csv", row.names = FALSE)
#' }
#'
#' @import dplyr
#' @import stringr
#' @export
classify_metabolite <- function(met) {
    if (str_detect(met, "^C\\d+(:\\d+)?(-DC|(-M)?-DC-M|-OH|-M-DC)?(\\s*\\(.*\\))?$")) {
        return("Acylcarnitines")
    } else if (str_detect(met, "^Cer\\(.*\\)$")) {
        return("Ceramides")
    } else if (str_detect(met, "^CE\\(\\d+:\\d+\\)$")) {
        return("Cholesteryl Esters")
    } else if (str_detect(met, "^DG\\([^)]+\\)$") || str_detect(met, "^DG-O?\\([^)]+\\)$")) {
        return("Diglycerides")
    } else if (str_detect(met, "^TG\\([^)]+\\)$")) {
        return("Triglycerides")
    } else if (str_detect(met, "^PC (aa|ae) C\\d+:\\d+$")) {
        return("Phosphatidylcholines")
    } else if (str_detect(met, "^lysoPC a C\\d+:\\d+$")) {
        return("Lysophosphatidylcholines")
    } else if (str_detect(met, "^SM( \\(OH\\))? C\\d+:\\d+$")) {
        return("Sphingomyelins")
    } else if (str_detect(met, "^HexCer\\(.*\\)$")) {
        return("Hexosylceramides")
    } else if (str_detect(met, "^Hex2Cer\\(.*\\)$")) {
        return("Dihexosylceramides")
    } else if (str_detect(met, "^Hex3Cer\\(.*\\)$")) {
        return("Trihexosylceramides")
    } else if (str_detect(met, "^FA\\(\\d+:\\d+\\)$") || met %in% c("DHA", "EPA")) {
        return("Fatty Acids")
    } else if (str_detect(met, "^(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|AA)$")) {
        return("Amino Acids")
    } else if (str_detect(met, "^(1-Met-His|3-Met-His|5-AVA|AABA|Ac-Orn|ADMA|alpha-AAA|Anserine|BABA|Betaine|Carnosine|Cit|Creatinine|Cystine|DOPA|HArg|HCys|Kynurenine|Met-SO|Nitro-Tyr|Orn|PheAlaBetaine|ProBetaine|Sarcosine|SDMA|t4-OH-Pro|c4-OH-Pro|Taurine|TrpBetaine|beta-Ala)$")) {
        return("Amino Acid-related")
    } else if (str_detect(met, "^(Dopamine|GABA|Histamine|PEA|Putrescine|Serotonin|Spermidine|Spermine)$")) {
        return("Biogenic Amines")
    } else if (str_detect(met, "^(CA|CDCA|DCA|GCA|GCDCA|GDCA|GLCA|GLCAS|GUDCA|TCA|TCDCA|TDCA|TLCA|TMCA)$")) {
        return("Bile Acids")
    } else if (str_detect(met, "^(AconAcid|DiCA\\(\\d+:\\d+\\)|HipAcid|Lac|OH-GlutAcid|Suc|AbsAcid)$")) {
        return("Carboxylic Acids")
    } else if (str_detect(met, "^(3-IAA|3-IPA|Ind-SO4|Indole)$")) {
        return("Indoles and Derivatives")
    } else if (str_detect(met, "^(Cortisol|Cortisone|DHEAS)$")) {
        return("Hormones")
    } else if (str_detect(met, "^(Hypoxanthine|Xanthine)$")) {
        return("Nucleobases and Related")
    } else if (str_detect(met, "^TMAO$")) {
        return("Amine Oxides")
    } else if (str_detect(met, "^Trigonelline$")) {
        return("Alkaloids")
    } else if (str_detect(met, "^p-Cresol-SO4$")) {
        return("Cresols")
    } else if (str_detect(met, "^H1$")) {
        return("Carbohydrates and Related")
    } else if (str_detect(met, "^Choline$")) {
        return("Vitamins and Cofactors")
    } else {
        return("Unknown")
    }
}
