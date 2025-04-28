library(tidyverse)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(scales)
library(missMDA) # for handling missing values
library(visdat)
library(dplyr, warn.conflicts = F)

options(width = 500)

adni_df <- read.csv('data-raw/ADNIMERGE_20Sep2024.csv') |>
    dplyr::select(-c(update_stamp))

# description of the dataset (missing values) ====
description <- "
The ADNI study is made up of four distinct phases – ADNI1, ADNIGO, ADNI2,
and ADNI3 – in each of which new participants were recruited while existing
participants from earlier phases continued to be monitored.
Our main dataset is the ADNI Merge dataset, from the Alzheimer's Disease
Neuroimaging Initiative. This is essentially a dataset combining key predictors
from all four phases, assembled using various sources of data within the ADNI
repository.

Categories include:
1. Cognitive Tests: some of the predictors, such as:
    - MMSE (Mini-Mental State Examination),
    - CDRSB (Clinical Dementia Rating Sum of Boxes),
    - ADAS (Alzheimer’s Disease Assessment Scale–Cognitive subscale) and
    - RAVLT (Rey Auditory Verbal Learning Test) come from cognitive tests
    that clinicians used to base the patient diagnoses on.
2. Demographic Variables: such as:
    - Age, Gender, Marriage status and Education levels.
3. Brain-Related Variables: these variables, such as
    - Hippocampus, Ventricles, WholeBrain, Entorhinal, Fusiform and MidTemp)
    measure various aspects of the brain.
4. Important Biomarkers: biomarkers, such as
    - A-beta (Amyloid Beta), Tau, APOE4 and FDG, are important proteins
    or biomarkers that are associated with Alzheimer’s disease
    or Mild Cognitive Impairment in a lot of medical literature
    about the disease.
"

"The table below provides a more detailed overview of the key variables
present in our dataset, the range of values that these variables take on
and the percentage of missing values, which is a key component to take
into consideration with this dataset."

adni_df <- adni_df |>
    dplyr::mutate_at(c("DX", "DX_bl"), ~na_if(.,"")) |>
    dplyr::mutate_at(c("DX", "DX_bl"), ~replace_na(.,"Unknown"))

adni_df_missingness <- adni_df |>
    dplyr::select(!ends_with("_bl"))

keep_cols <- c(
    colnames(adni_df_missingness)[colSums(is.na(adni_df_missingness)) > 100],
    "VISCODE", "COLPROT"
)

## missing Values
"Given the strict data-collection protocol followed and strict criteria
for selecting patients so as to prevent drop-out, we expect that missing
values are typically Missing At Random (MAR); missingness is due solely
to observable factors such as the follow-up time (i.e. not all variables
are re-collected at every follow-up), the patient’s initial diagnosis, and
the particular ADNI phase (given that the data are currently in longform).
Especially where the missingness is due to differing procedures carried out
for patients with different baseline diagnoses, we must ensure that the way
we handle missing values does not introduce bias in our model."

# baseline observations and sort by protocol (ADNI stage)
adni_anymissing <- adni_df_missingness |>
    dplyr::filter(VISCODE == "bl") |>
    dplyr::select(all_of(keep_cols)) |>
    dplyr::arrange(VISCODE, COLPROT)

visdat::vis_miss(adni_anymissing |>
                dplyr::select(-c(VISCODE)),
                cluster = TRUE,
                sort_miss = FALSE,
                facet = COLPROT)

"The plot above shows that certain categories, such as
the Everyday Cognition (Ecog) are missing in general for the ADNI 1 phase,
whereas the brain related predictors are missing for ADNI 3. One important
observation is that the important biolarkers (such as Amyloid Beta, Tau
and PTau) are missing partially throughout the data. The reason for this
is that these are taken from a patient's Cerebraospinal Fluid (CSF),
which requires an invasive procedure to obtain. Therefore, not all patients
undergo this procedure during a baseline test and we see a lot of missing
values here."

a <- 2

# demographic Information ====
description_demographics <- "
The patient demographics helped shape some of our goals and research questions.
First, we noticed that all participants were at or above the age of 55.
This means that our ability to make an 'early' diagnosis is limited, since
many of the participants got a screening done since they exhibited symptoms
of cognitive impairment of some form. Second, we notice that a majority
of our population is white and married We understand that this means that
our findings do not extend to a larger population.
"

age_data <- na.omit(adni_df$AGE)
p1 <- ggplot(data.frame(AGE = age_data), aes(x = AGE)) +
    geom_histogram(bins = 10, fill = "skyblue", alpha = 0.7) +
    labs(title = "Age", y = "Frequency") +
    theme_minimal(base_size = 10)

# EDUCATION
education_data <- na.omit(adni_df$PTEDUCAT)
p2 <- ggplot(data.frame(PTEDUCAT = education_data), aes(x = PTEDUCAT)) +
    geom_histogram(bins = 10, fill = "lightgreen", alpha = 0.7) +
    labs(title = "Education (years)", y = "Frequency") +
    theme_minimal()

# GENDER
gender_data <- na.omit(adni_df$PTGENDER)
p3 <- ggplot(data.frame(PTGENDER = gender_data), aes(x = PTGENDER)) +
    geom_bar(fill = "lightcoral", alpha = 0.7) +
    labs(title = "Gender", x = "Gender", y = "Count") +
    theme_minimal()

# RACE
race_data <- na.omit(adni_df$PTRACCAT)
p4 <- ggplot(data.frame(PTRACCAT = race_data), aes(x = PTRACCAT)) +
    geom_bar(fill = "gold", alpha = 0.7) +
    labs(title = "Race", x = "Race", y = "Count") +
    theme_minimal() +
    coord_flip()

# MARITAL STATUS
marital_data <- na.omit(adni_df$PTMARRY)
p5 <- ggplot(data.frame(PTMARRY = marital_data), aes(x = PTMARRY)) +
    geom_bar(fill = "orchid", alpha = 0.7) +
    labs(title = "Marital Status", x = "Status", y = "Count") +
    theme_minimal() +
    coord_flip()

# Combine plots in a grid
grid.arrange(p1, p2, p3, p4, p5, nrow = 1)


options(repr.plot.width = 5, repr.plot.height = 2)

library(patchwork)

combined_plot <- p1 + p2 + p3 + p4 + p5 +
    plot_layout(guides = "collect") &
    theme(legend.position = "right",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 15, face = "bold"),
          plot.title = element_text(size = 15, face = "bold"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "grey85"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 15,face = "bold"))

combined_plot <- combined_plot + plot_annotation(tag_levels = c('A', '1'))

combined_plot

# description diagnosis ====
description_diagnoses <- "
Initial Diagnoses
The baseline diagnoses are encoded in two ways in the ADNI Merge Dataset.
The first encoding is as follows:

- CN: Cognitively Normal
- MCI: Mild Cognitive Impairment
- Dementia: Alzheimer's Disease or other Dementia.

The second encoding is as follows:
- CN: Cognitively Normal
- EMCI: Early Mild Cognitive Impairment
- LMCI: Late Mild Cognitive Impairment
- SMC: Significant Memory Concerns
- AD: Alzheimer's Disease

I choose the first encoding, changing 'Dementia' to 'AD' since there is
an equivalency in the encoding of these categories. The second encoding is
present only in the baseline diagnoses and not in the subsequent diagnoses,
and part of our modeling is to predict future decline, so we choose the first
encoding to remain consistent between diagnoses at different points of time.
"

# filter dataset for relevant subset
sub_df <- adni_df |> dplyr::filter(M == 0)

# first encoding
dxs <- sub_df$DX |>
    dplyr::recode('Dementia' = 'AD', 'CN' = 'CN', 'MCI' = 'MCI') |>
    as.factor()
dx_count <- data.frame(Diagnosis = levels(dxs), Count = table(dxs))

p1 <- ggplot(dx_count, aes(x = Diagnosis, y = Count.Freq)) +
    geom_bar(stat = "identity", fill = "skyblue", alpha = 0.7) +
    labs(title = "Baseline diagnoses (1st encoding)",
         x = "Diagnosis", y = "Frequency") +
    theme_minimal()

# second encoding
dx_bl_count <- data.frame(Diagnosis = levels(as.factor(sub_df$DX_bl)),
                          Count = table(as.factor(sub_df$DX_bl)))

p2 <- ggplot(dx_bl_count, aes(x = Diagnosis, y = Count.Freq)) +
    geom_bar(stat = "identity", fill = "lightgreen", alpha = 0.7) +
    labs(title = "Baseline diagnoses (2nd encoding)",
         x = "Diagnosis", y = "Frequency") +
    theme_minimal()


grid.arrange(p1, p2, ncol = 2)


# description predictors for initial diagnosis ====
description_predictors <- "
Important Predictors for Initial Diagnosis
From a histogram of baseline values for all predictors in our dataset
conditional on initial diagnosis, a few predictors stood out to us,
as being promising indicators of baseline diagnosis. These pertained
to examination scores, which doctors heavily rely on to make their initial
diagnoses. This influenced our decision to explore one such heavily
influential examination, the MMSE (Mini-Mental State Examination)
in particular.

The histograms below correspond to four cognitive tests that show a promising
separability between the three classes at baseline levels. We further explore
these four cognitive tests and three tests from the ADAS (Alzheimer's Disease
Assessment Scale) test in a pair plot that also shows promise in a combination
of cognitive tests being used to predict initial diagnoses.
"

color_dict <- c("CN" = "red", "AD" = "blue", "Dementia" = "blue",
                "LMCI" = "green", "MCI" = "green", "EMCI" = "orange",
                "SMC" = "purple")
name_dict <- c("CN" = "CN", "Dementia" = "AD", "MCI" = "MCI")

# baseline
baseline_df <- adni_df %>% filter(M == 0)

variables <- c("MMSE", "CDRSB", "RAVLT_immediate", "mPACCdigit")

plots <- lapply(variables, function(var) {
    df_copy <- baseline_df |>
        dplyr::select(M, DX, all_of(var)) |>
        drop_na()

    ggplot(df_copy, aes_string(x = var, fill = "DX")) +
        geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
        scale_fill_manual(values = color_dict) +
        labs(title = var, x = var, y = "Frequency") +
        theme_minimal() +
        theme(legend.title = element_blank())
})

grid.arrange(grobs = plots, ncol = 2)


library(GGally)

# description relationships bxt predictors ====
description_pairplot <- "
Exploring Relationships Between Predictors
This pair plot visualizes the relationships between baseline cognitive test
scores and their potential in distinguishing diagnoses (CN, MCI, AD).
Each point is color-coded based on the diagnosis category.
"

name_dict <- c("CN" = "CN", "Dementia" = "AD", "MCI" = "MCI")

# baseline data and map diagnoses
baseline_df <- adni_df |>
    dplyr::filter(M == 0) |>
    dplyr::mutate(DX_map = recode(DX, !!!name_dict))

cols <- c("CDRSB_bl", "MMSE_bl", "ADASQ4_bl", "ADAS11_bl", "ADAS13_bl",
          "RAVLT_immediate_bl", "mPACCdigit_bl")

baseline_df <- baseline_df |>
    dplyr::filter(DX != "Unknown", DX_bl != "Unknown")

ggpairs(
    baseline_df,
    columns = which(names(baseline_df) %in% cols), # Select relevant columns
    aes(color = DX_map, fill = DX_map, alpha = 0.6),
    upper = list(continuous = "points"),
    # upper = "blank",
    lower = list(continuous = "smooth"),
    diag = list(continuous = "densityDiag")
    # switch = "both"
) +
    scale_color_manual(values = c("CN" = "red", "MCI" = "green", "AD" = "blue")) +
    scale_fill_manual(values = c("CN" = "red", "MCI" = "green", "AD" = "blue")) +
    theme_gray(base_size = 14) +
    labs(title = "Pair-plot of baseline cognitive tests",
         subtitle = "Color-coded by diagnosis category",
         color = "Diagnosis",
         fill = "Diagnosis") +
    theme(strip.background = element_blank(), strip.placement = "outside")


# description dx over time ====
description_diagnoses_over_time <- "
Diagnoses Over Time
This section examines how the distribution of diagnoses changes over time,
focusing on transitions from Cognitively Normal (CN) to Mild Cognitive
Impairment (MCI) and from MCI to Alzheimer's Disease (AD).
A line plot visualizes the number of diagnoses across different time points.

In addition to being able to predict initial diagnosis, we also aim to
predict future decline i.e. patients who were initially diagnosed as
Cognitively Normal and later were diagnosed with Mild Cognitive Impairment,
or patients who were initially diagnosed with Mild Cognitive Impairment
and later became diagnosed with Alzheimer's Disease. In order to do this,
we have a few different visuals to get a sense of the distribution of
diagnoses over time and how MMSE scores (an important determiner of
diagnosis) vary over time for patients.
"

name_dict <- c("CN" = "Cognitively Normal", "Dementia" = "Alzheimer's Disease",
               "MCI" = "Mild Cognitive Impairment")
adni_df$DX <- recode(adni_df$DX, !!!name_dict)

kde_df <- adni_df |>
    dplyr::select(RID, M, DX)

# the group by time points (M) and count diagnoses
table_df <- kde_df |>
    dplyr::group_by(M) |>
    dplyr::count(DX) |>
    dplyr::filter(DX != "Unknown") |>
    tidyr::pivot_wider(names_from = DX, values_from = n, values_fill = 0)

ggplot(data = table_df, aes(x = M)) +
    geom_line(aes(y = `Cognitively Normal`, color = "CN"), size = 1.2) +
    geom_line(aes(y = `Mild Cognitive Impairment`, color = "MCI"), size = 1.2) +
    geom_line(aes(y = `Alzheimer's Disease`, color = "AD"), size = 1.2) +
    labs(
        title = "Number of diagnoses over time",
        x = "Time after baseline (Months)",
        y = "Number of diagnoses",
        color = "Diagnosis"
    ) +
    theme_classic(base_size = 15) +
    theme(
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "inside", legend.position.inside =  c(.9, .7),
        legend.box.background = element_rect(color = "red", size = 2),
    )

"We observe that the total number of diagnoses of each kind drops over time,
which may perhaps be a result of patients not carrying out a follow up study
for a number of different reasons."

# description mmse over time====
description_mmse <- "
Trajectories of MMSE Scores Over Time
This plot illustrates the trajectories of MMSE scores for each diagnosis
group (CN, MCI, AD) over time. It helps identify patterns in cognitive
decline across groups.

Next, we show how MMSE scores vary over time for participants. We plot
trajectories of MMSE scores for each of the three diagnosis types.
We anticipate there to be a gradual decline in MMSE scores for the AD
and MCI categories and perhaps a smaller decline for the CN category
but this plot helps discern between the different groups based on rate
of dropoff of MMSE scores.
"

color_dict <- c("CN" = "red", "MCI" = "green", "AD" = "blue")

name_dict <- c("CN" = "Cognitively Normal", "AD" = "Alzheimer's Disease",
               "MCI" = "Mild Cognitive Impairment")

traj_df <- adni_df |>
    dplyr::select(RID, M, MMSE, DX_bl, DX) |>
    drop_na() |>
    dplyr::filter(DX != "Unknown") |>
    dplyr::group_by(RID) |>
    summarise(
        Months = list(M),
        MMSE_scores = list(MMSE),
        DX_bl = first(DX_bl),
        DX = first(DX)
    ) |>
    dplyr::mutate("DX" = case_when(
        DX %in%  "Alzheimer's Disease" ~ "AD",
        DX %in%  "Mild Cognitive Impairment" ~ "MCI",
        DX %in%  "Cognitively Normal" ~ "CN"
        ))

ggplot(data = traj_df, aes(group = RID)) +
    geom_line(
        data = traj_df |> unnest(c(Months, MMSE_scores)) |>
            arrange(DX, RID, Months),
        aes(x = Months, y = MMSE_scores, color = DX),
        alpha = 0.7
    ) +
    scale_color_manual(values = color_dict, labels = name_dict) +
    facet_wrap(~ DX, labeller = as_labeller(name_dict), scales = "free_y") +
            labs(
        title = "Trajectories of MMSE scores over time",
        x = "Months after baseline",
        y = "MMSE levels",
        color = "Diagnosis group",
        subtitle = "Patterns of cognitive decline across diagnosis groups"
    ) +
    scale_y_continuous(limits = c(0, 31)) +
    theme_gray(base_size = 14) +
    theme(
        strip.text = element_text(size = 15, face = "bold"),
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 14)
    ) +
    theme(strip.background = element_blank(), strip.placement = "outside")


"As hypothesized, the dropoff of MMSE scores for the MCI group is faster
and larger than that for the CN group. It is striking to observe that many
patients initially diagnosed with AD have no observations after 3 years."

# description disease change between 1st and most recent visit ====
description_dx_change <- "
Change in Diagnoses Between First and Most Recent Visit
This plot shows the counts of initial diagnoses and most recent diagnoses
for patients who have multiple observations. It highlights transitions in
cognitive status over time.

Using our plots above, we turn to the question of predicting change
in diagnoses between first and most recent visit for a particular patient.
We subset our dataset to include only patients who have multiple observations
and show the number of diagnoses of each kind at baseline and at patients'
most recent visit. One extremely important caveat is that the most recent
visit can be 2 years later or 12 years later, and the time between first
and last visit is an extremely important predictor that we keep in mind
when modeling.
"

name_dict <- c("CN" = "Cognitively Normal", "Dementia" = "AD",
               "MCI" = "Mild Cognitive Impairment")

# identifying last visit for each patient
last_visits <- adni_df |>
    dplyr::group_by(PTID) |>
    dplyr::filter(M == max(M)) |>
    dplyr::ungroup()

# sub-setting initial and last visit diagnoses
initial_dx <- adni_df |>
    dplyr::filter(M == 0) |>
    dplyr::select(RID, DX) |>
    dplyr::rename(DX_initial = DX)

last_dx <- last_visits |>
    dplyr::select(RID, DX) |>
    dplyr::rename(DX_recent = DX)

# merging initial and last diagnoses
dx_change <- initial_dx |>
    dplyr::inner_join(last_dx, by = "RID") |>
    dplyr::mutate(
        DX_initial = recode(DX_initial, !!!name_dict),
        DX_recent = recode(DX_recent, !!!name_dict),
        Transition = paste(DX_initial, DX_recent, sep = " → ")
    )

# count diagnoses
dx_counts <- dx_change |>
    summarize(
        Initial = list(table(DX_initial)),
        Recent = list(table(DX_recent))
    ) |>
    tidyr::pivot_longer(cols = everything(),
                        names_to = "Type",
                        values_to = "Counts") |>
    unnest(Counts) |>
    dplyr::mutate(
        DX = names(Counts)
    ) |>
    dplyr::select(DX, Type, Counts)

ggplot(data = dx_counts, aes(x = DX, y = Counts, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7, width = 0.5) +
    scale_fill_manual(values = c("Initial" = "skyblue",
                                 "Recent" = "steelblue")) +
    theme_minimal(base_size = 14) +
    labs(
        title = "Initial and most recent diagnoses",
        x = "Diagnosis group",
        y = "Number of patients",
        fill = "Diagnosis timepoint",
        subtitle = "Comparison of diagnoses between Baseline and Last visit"
    ) +
    theme(
        axis.text.x = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
    )
    # scale_y_continuous(limits = c(0, 800))

"A general pattern seems to be that the total number of MCI and CN diagnoses
has decreased while the number of AD diagnoses has increased over time,
which is unsurprising.

We create a table, with each entry representing a participant's initial
and most recent diagnoses. A majority of participants start and end with
the same diagnoses, however a subset of participants do decline,
especially from MCI to AD. A non-trivial number of participants move
from MCI to CN, which is worth noting. Of the 1648 participants with
multiple observations, only 204 decline, which is about 12.3% of
our total sample. This is an important finding that suggests
that we should perform classification with balanced class weights
in our modeling phase.

The Sankey Diagram below illustrates changes in diagnoses and highlights
some major trends: that a majority of patients have the same initial
and most recent diagnoses."

library(networkD3)

# description transitions ====
description_transitions <- "
Analysis of Diagnosis Transitions
This section examines changes in participant diagnoses between baseline
and their most recent visits.
A Sankey Diagram highlights transitions, emphasizing that most patients
retain their original diagnosis,
while a smaller subset transitions, particularly from MCI to AD or MCI to CN.
"

name_dict <- c("Cognitively Normal" = "CN", "Alzheimer's Disease" = "AD",
               "Mild Cognitive Impairment" = "MCI")

dx_transition <- dx_change |>
    # dplyr::filter(DX_initial != "Unknown", DX_recent != "Unknown") |>
    dplyr::group_by(Transition) |>
    dplyr::summarise(Count = n()) |>
    dplyr::ungroup() |>
    separate(Transition, into = c("Source", "Target"), sep = " → ") |>
    dplyr::mutate(
        Source = recode(Source, !!!name_dict),
        Target = recode(Target, !!!name_dict)
    )

library(ggalluvial)
ggplot(dx_transition,
       aes(axis1 = Source, axis2 = Target, y = Count)) +
    geom_alluvium(aes(fill = Source), knot.pos = 0.3, width = 1/12, alpha = 0.7) +
    geom_stratum(aes(fill = Source), width = 1/12, fill = "grey39", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Source", "Target"), expand = c(.05, .05)) +
    theme(legend.position = "none") +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    ggtitle("Change in Diagnoses")

# creating node and link data for Sankey Diagram
nodes <- data.frame(
    name = unique(c(dx_transition$Source, dx_transition$Target))
)

links <- dx_transition |>
    dplyr::mutate(
        source = match(Source, nodes$name) - 1,
        target = match(Target, nodes$name) - 1
    ) |>
    dplyr::select(source, target, Count) |>
    dplyr::filter(if_all(source, ~ . != target))

# Sankey Diagram
sankey <- sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "Count",
    NodeID = "name",
    units = "Patients",
    fontSize = 14,
    nodeWidth = 30,
    sinksRight = TRUE
)

sankey


dx_change <- dx_change |>
    dplyr::mutate(
        DX_recent = ifelse(is.na(DX_recent) | DX_recent == "", "Unknown", DX_recent),
        DX_initial = ifelse(is.na(DX_initial) | DX_initial == "", "Unknown", DX_initial)
    )

transition_table <- dx_change |>
    dplyr::filter(DX_initial != "Unknown", DX_recent != "Unknown") |>
    dplyr::group_by(DX_initial, DX_recent) |>
    summarise(Count = n(), .groups = "drop") |>
    tidyr::pivot_wider(
        names_from = DX_recent,
        values_from = Count,
        values_fill = list(Count = 0)
    ) |>
    dplyr::rename(`Initial Diagnosis` = DX_initial)

library(ggalluvial)

