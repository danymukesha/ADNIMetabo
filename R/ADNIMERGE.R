#### ADNIMERGE ####

#' The ADNIMERGE class
#'
#' The ADNIMERGE object is a representation of the ADNIMERGE dataset,
#' containing demographic, clinical, biomarker, imaging, cognitive, and
#' other relevant data.
#' @slot demographic_data a data.frame or tbl_df of
#' \link{demographic_data}.
#' @slot baseline_clinical_data a data.frame or tbl_df of
#' \link{baseline_clinical_data}.
#' @slot follow_up_clinical_data a data.frame or tbl_df of
#' \link{follow_up_clinical_data}.
#' @slot biomarker_data a data.frame or tbl_df of
#' \link{biomarker_data}.
#' @slot imaging_data a data.frame or tbl_df of \link{imaging_data}.
#' @slot cognitive_composite_data a data.frame or tbl_df of
#' \link{cognitive_composite_data}.
#' @slot ecog_scores_data a data.frame or tbl_df of
#' \link{ecog_scores_data}.
#' @slot RID_to_PTID a data.frame or tbl_df linking \link{RID_to_PTID}.
#' @slot logs Log information of data processing steps.
#' @slot misc_data ...
#' @importFrom methods setClass
#' @name ADNIMERGE-class
#' @rdname ADNIMERGE-class
#' @export
#' @return An ADNIMERGE class object.
#' @seealso \link{ADNIMERGE}
#'
ADNIMERGE <- setClass(
    Class = 'ADNIMERGE',
    slots = list(
        'demographic_data' = 'tbl_df',
        'baseline_clinical_data' = 'tbl_df',
        'follow_up_clinical_data' = 'tbl_df',
        'biomarker_data' = 'tbl_df',
        'imaging_data' = 'tbl_df',
        'cognitive_composite_data' = 'tbl_df',
        'ecog_scores_data' = 'tbl_df',
        'RID_to_PTID' = 'tbl_df',
        'logs' = 'character',
        'misc_data' = 'list'
    )
)
