#' Transform column names of the FIA dataset
#'
#' This function applies a series of transformations to the column names
#' of a data frame. It performs specific replacements for common column name
#' patterns and handles special cases for various lipid and compound names
#' (e.g., Ceramides, HexCer, SM, PC, etc.).
#'
#' @param old_name string parameter representing the old column name.
#'
#' @return string parameter with the transformed column name.
#'
#' @details The function handles the following transformations:
#' \itemize{
#'   \item Special cases for specific column names like "PTID" to "RID",
#'   "Customer_Sample_Identification" to "Customer Sample Identification", etc.
#'   \item Replacement of underscores with colons for compound names
#'   like "C3_1" to "C3:1".
#'   \item Transformation of OH/DC notation compounds
#'   (e.g., "_OH__" to "-OH (", "_DC__" to "-DC (").
#'   \item Special formatting for Ceramide, HexCer, Hex2Cer, and Hex3Cer names.
#'   \item SM, lysoPC, and PC aa/ae transformations.
#'   \item Support for handling DG and TG names, including modifications
#'   to handle isotopic forms.
#'   \item Reformatting of CE names.
#' }
#'
#' @examples
#' # e.g. of applying the transformation function to a set of column names:
#' old_names <- c("PTID", "Customer_Sample_Identification", "C3_1",
#'     "Cer_d18:1_18:1")
#' # old_names <- colnames(FIA_Metabo)
#' new_names <- sapply(old_names, transform_colnames, USE.NAMES = FALSE)
#' print(new_names)
#'
#' @export
transform_colnames <- function(old_name) {
    if (old_name == "PTID") return("RID")
    if (old_name == "Customer_Sample_Identification")
        return("Customer Sample Identification")
    if (old_name == "Plate_Bar_Code") return("Plate Bar Code")
    if (old_name == "Sample_Bar_Code") return("Sample Bar Code")
    if (old_name == "Well_Position") return("Well Position")
    if (old_name == "Sample_Volume") return("Sample Volume")
    if (old_name == "Run_Number") return("Run Number")
    if (old_name == "Injection_Number") return("Injection Number")

    # the pattern of underscores to colons for simple compounds (e.g., C3_1 to C3:1)
    if (grepl("^C[0-9]+_[0-9]+$", old_name)) {
        return(gsub("_", ":", old_name))
    }

    # compounds with OH/DC notations
    if (grepl("_OH__", old_name)) {
        return(gsub("_OH__", "-OH (", old_name) %>% gsub("_", ":", .) %>% paste0(")"))
    }
    if (grepl("_OH_", old_name)) {
        return(gsub("_OH", "-OH", old_name) %>% gsub("_", ":", .))
    }
    if (grepl("_DC__", old_name)) {
        return(gsub("_DC__", "-DC (", old_name) %>% gsub("_", ":", .) %>% paste0(")"))
    }
    if (grepl("_DC_", old_name)) {
        return(gsub("_DC", "-DC", old_name) %>% gsub("_", ":", .))
    }

    # Ceramide names
    if (grepl("^Cer_d", old_name)) {
        return(gsub("^Cer_d", "Cer(d", old_name) %>%
                   gsub("_", ":", .) %>%
                   gsub(":([0-9]+):([0-9]+)$", "/\\1:\\2)", .) %>%
                   gsub(":([0-9]+):([0-9]+):", "/\\1:\\2(", .))
    }

    # HexCer names
    if (grepl("^HexCer_d", old_name)) {
        return(gsub("^HexCer_d", "HexCer(d", old_name) %>%
                   gsub("_", ":", .) %>%
                   gsub(":([0-9]+):([0-9]+)$", "/\\1:\\2)", .) %>%
                   gsub(":([0-9]+):([0-9]+):", "/\\1:\\2(", .))
    }

    # Hex2Cer names
    if (grepl("^Hex2Cer_d", old_name)) {
        return(gsub("^Hex2Cer_d", "Hex2Cer(d", old_name) %>%
                   gsub("_", ":", .) %>%
                   gsub(":([0-9]+):([0-9]+)$", "/\\1:\\2)", .) %>%
                   gsub(":([0-9]+):([0-9]+):", "/\\1:\\2(", .))
    }

    # Hex3Cer names
    if (grepl("^Hex3Cer_d", old_name)) {
        return(gsub("^Hex3Cer_d", "Hex3Cer(d", old_name) %>%
                   gsub("_", ":", .) %>%
                   gsub(":([0-9]+):([0-9]+)$", "/\\1:\\2)", .) %>%
                   gsub(":([0-9]+):([0-9]+):", "/\\1:\\2(", .))
    }

    # SM names
    if (grepl("^SM__OH__", old_name)) {
        return(gsub("^SM__OH__", "SM (OH) ", old_name) %>% gsub("_", ":", .))
    }
    if (grepl("^SM_", old_name)) {
        return(gsub("^SM_", "SM ", old_name) %>% gsub("_", ":", .))
    }

    # lysoPC names
    if (grepl("^lysoPC_a_", old_name)) {
        return(gsub("^lysoPC_a_", "lysoPC a ", old_name) %>% gsub("_", ":", .))
    }

    # PC aa/ae names
    if (grepl("^PC_aa_", old_name)) {
        return(gsub("^PC_aa_", "PC aa ", old_name) %>% gsub("_", ":", .))
    }
    if (grepl("^PC_ae_", old_name)) {
        return(gsub("^PC_ae_", "PC ae ", old_name) %>% gsub("_", ":", .))
    }

    # DG names
    if (grepl("^DG_", old_name)) {
        new_name <- gsub("^DG_", "DG(", old_name)
        # DG-O separately
        if (grepl("_O_", new_name)) {
            new_name <- gsub("_O_", "-O(", new_name)
        }
        new_name <- gsub("_", ":", new_name)
        new_name <- gsub(":([0-9]+):([0-9]+)$", "_\\1:\\2)", new_name)
        new_name <- gsub(":([0-9]+):([0-9]+):", "_\\1:\\2_", new_name)
        return(new_name)
    }

    # TG names
    if (grepl("^TG_", old_name)) {
        new_name <- gsub("^TG_", "TG(", old_name)
        new_name <- gsub("_", ":", new_name)
        new_name <- gsub(":([0-9]+):([0-9]+)$", "_\\1:\\2)", new_name)
        new_name <- gsub(":([0-9]+):([0-9]+):", "_\\1:\\2_", new_name)
        return(new_name)
    }

    # CE names
    if (grepl("^CE_", old_name)) {
        return(gsub("^CE_", "CE(", old_name) %>%
                   gsub("_", ":", .) %>%
                   paste0(")"))
    }

    # default case - return the original name
    return(old_name)
}
