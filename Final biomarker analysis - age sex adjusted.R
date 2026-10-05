# ============================================================
# AGE- AND SEX-ADJUSTED BIOMARKER ANALYSIS
#
# THREE CLINICAL GROUPS:
#   1. Control
#   2. ICU (Non-pneumonia)
#   3. ICU (Pneumonia)
#
# PRIMARY ADJUSTED ANALYSIS:
#   All biomarkers except PSA/PZP:
#       log2(Biomarker) ~ Group + Age + Sex
#
#   Kallikrein 3/PSA (MALES ONLY):
#       log2(PSA) ~ Group + Age
#
#   PZP (FEMALES ONLY):
#       log2(PZP) ~ Group + Age
#
# WHY LOG2?
#   All 390 biomarker values in the uploaded dataset are > 0.
#   Modeling log2 concentrations helps with right-skew and makes
#   pairwise model estimates directly interpretable as fold changes:
#       Adjusted fold change = 2^(adjusted log2 difference)
#
# MULTIPLE-TESTING FAMILIES (same structure as prior script):
#   1. Omnibus adjusted Group tests:
#        BH across all biomarkers
#   2. ICU-v-Control adjusted pairwise comparisons:
#        NonP vs Control + Pneu vs Control pooled together
#   3. Pneu vs NonP adjusted pairwise comparisons:
#        separate BH family
#
# SIGNIFICANCE:
#   Omnibus_Group_FDR < 0.05 AND Pairwise_FDR < 0.05
#
# DESCRIPTIVES:
#   Raw N, median, Q1, Q3 are retained for interpretability.
#   Inference is from the adjusted regression/ANCOVA models.
# ============================================================


# ============================================================
# 1. INSTALL PACKAGES IF NEEDED
# ============================================================

packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "emmeans",
  "openxlsx"
)

new_packages <- packages[
  !(packages %in% installed.packages()[, "Package"])
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}


# ============================================================
# 2. LOAD PACKAGES
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(emmeans)
library(openxlsx)


# ============================================================
# 3. SELECT INPUT FILE + SET OUTPUT LOCATION
# ============================================================

cat(
  "\nPlease select the Excel file containing the biomarker data...\n"
)

# Opens a file browser so you can select the workbook
input_file <- file.choose()

# Get the folder where the selected workbook is located
input_folder <- dirname(input_file)

# Save the results into that same folder
output_file <- file.path(
  input_folder,
  "Final biomarker analysis - age sex adjusted.xlsx"
)

cat(
  "\nInput file:\n",
  input_file,
  "\n\nOutput file will be saved to:\n",
  output_file,
  "\n\n",
  sep = ""
)
# ============================================================
# 4. READ RAW DATA
# ============================================================

dat <- readxl::read_excel(
  input_file,
  sheet = "Final with clinical data"
)


# ============================================================
# 5. IDENTIFY BIOMARKERS
#    Every column after "Type of Bacteria"
# ============================================================

type_bacteria_col <- match(
  "Type of Bacteria",
  names(dat)
)

if (is.na(type_bacteria_col)) {
  stop("Could not find the column 'Type of Bacteria'.")
}

biomarker_vars <- names(dat)[
  (type_bacteria_col + 1):ncol(dat)
]

cat(
  "\nBiomarkers identified: ",
  length(biomarker_vars),
  "\n",
  sep = ""
)

if (length(biomarker_vars) != 390) {
  warning(
    paste0(
      "Expected 390 biomarkers but found ",
      length(biomarker_vars),
      ". Check workbook structure."
    )
  )
}


# ============================================================
# 6. VERIFY REQUIRED COLUMNS AND SPECIAL BIOMARKERS
# ============================================================

required_columns <- c(
  "Deidentified",
  "Group",
  "Age",
  "Sex"
)

missing_required <- required_columns[
  !required_columns %in% names(dat)
]

if (length(missing_required) > 0) {
  stop(
    paste0(
      "Missing required column(s): ",
      paste(missing_required, collapse = ", ")
    )
  )
}

special_biomarkers <- c(
  "Kallikrein 3/PSA",
  "PZP"
)

missing_special <- special_biomarkers[
  !special_biomarkers %in% biomarker_vars
]

if (length(missing_special) > 0) {
  stop(
    paste0(
      "Could not find the following biomarker(s): ",
      paste(missing_special, collapse = ", ")
    )
  )
}


# ============================================================
# 7. CLEAN GROUP, AGE, AND SEX
# ============================================================

group_order <- c(
  "Control",
  "ICU (Non-pneumonia)",
  "ICU (Pneumonia)"
)

dat <- dat %>%
  dplyr::mutate(
    Group = trimws(as.character(Group)),
    Age = suppressWarnings(as.numeric(Age)),
    Sex_raw = trimws(as.character(Sex)),
    Sex_clean = dplyr::case_when(
      tolower(Sex_raw) %in% c("male", "m") ~ "Male",
      tolower(Sex_raw) %in% c("female", "f") ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(
    Group %in% group_order
  )

dat$Group <- factor(
  dat$Group,
  levels = group_order
)

dat$Sex_clean <- factor(
  dat$Sex_clean,
  levels = c("Female", "Male")
)

cat(
  "\nAge missing: ",
  sum(is.na(dat$Age)),
  "\n",
  sep = ""
)

cat(
  "Sex missing: ",
  sum(is.na(dat$Sex_clean)),
  "\n",
  sep = ""
)

print(
  table(
    dat$Group,
    dat$Sex_clean,
    useNA = "ifany"
  )
)


# ============================================================
# 8. CONVERT BIOMARKERS TO NUMERIC
# ============================================================

dat <- dat %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(biomarker_vars),
      ~ suppressWarnings(
        as.numeric(.x)
      )
    )
  )


# ============================================================
# 9. VERIFY POSITIVE VALUES FOR LOG2 MODELING
# ============================================================

nonpositive_check <- dat %>%
  dplyr::select(
    dplyr::all_of(biomarker_vars)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Biomarker",
    values_to = "Value"
  ) %>%
  dplyr::filter(
    !is.na(Value),
    Value <= 0
  )

if (nrow(nonpositive_check) > 0) {
  stop(
    paste0(
      "At least one biomarker contains zero or negative values. ",
      "This script uses log2 biomarker values, so inspect those values first."
    )
  )
}


# ============================================================
# 10. CREATE LONG-FORM BIOMARKER DATA
#
# PSA: males only
# PZP: females only
# All others: both sexes
# ============================================================

long_dat <- dat %>%
  dplyr::select(
    Deidentified,
    Group,
    Age,
    Sex = Sex_clean,
    dplyr::all_of(biomarker_vars)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(biomarker_vars),
    names_to = "Biomarker",
    values_to = "Value"
  ) %>%
  dplyr::filter(
    !is.na(Value)
  ) %>%
  dplyr::filter(
    !(
      Biomarker == "Kallikrein 3/PSA" &
        Sex != "Male"
    ),
    !(
      Biomarker == "PZP" &
        Sex != "Female"
    )
  ) %>%
  dplyr::mutate(
    Log2_Value = log2(Value),
    
    Sex_restriction =
      dplyr::case_when(
        Biomarker ==
          "Kallikrein 3/PSA" ~
          "Male only",
        
        Biomarker ==
          "PZP" ~
          "Female only",
        
        TRUE ~
          "Both sexes"
      ),
    
    Adjustment =
      dplyr::case_when(
        Biomarker ==
          "Kallikrein 3/PSA" ~
          "Age adjusted; males only",
        
        Biomarker ==
          "PZP" ~
          "Age adjusted; females only",
        
        TRUE ~
          "Age + sex adjusted"
      )
  )


# ============================================================
# 11. SPECIAL BIOMARKER SAMPLE-SIZE CHECK
# ============================================================

special_check <- long_dat %>%
  dplyr::filter(
    Biomarker %in%
      special_biomarkers
  ) %>%
  dplyr::count(
    Biomarker,
    Sex_restriction,
    Group,
    Sex,
    name = "N"
  ) %>%
  dplyr::arrange(
    Biomarker,
    Group,
    Sex
  )

cat(
  "\n============================================\n",
  "SPECIAL BIOMARKER SEX CHECK\n",
  "============================================\n\n",
  sep = ""
)

print(
  special_check
)


# ============================================================
# 12. RAW GROUP DESCRIPTIVE STATISTICS
#
# These are NOT age/sex adjusted.
# They are included for descriptive reporting.
# ============================================================

group_descriptives <- long_dat %>%
  
  dplyr::group_by(
    Biomarker,
    Sex_restriction,
    Adjustment,
    Group
  ) %>%
  
  dplyr::summarise(
    
    N =
      sum(
        !is.na(Value)
      ),
    
    Median =
      median(
        Value,
        na.rm = TRUE
      ),
    
    Q1 =
      quantile(
        Value,
        probs = 0.25,
        na.rm = TRUE,
        names = FALSE
      ),
    
    Q3 =
      quantile(
        Value,
        probs = 0.75,
        na.rm = TRUE,
        names = FALSE
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    Biomarker,
    Group
  )


# ============================================================
# 13. CREATE WIDE DESCRIPTIVE TABLE
# ============================================================

descriptive_wide <- group_descriptives %>%
  
  dplyr::mutate(
    
    Group_short =
      dplyr::case_when(
        
        as.character(Group) ==
          "Control" ~
          "Control",
        
        as.character(Group) ==
          "ICU (Non-pneumonia)" ~
          "NonP",
        
        as.character(Group) ==
          "ICU (Pneumonia)" ~
          "Pneu"
      )
  ) %>%
  
  dplyr::select(
    Biomarker,
    Sex_restriction,
    Adjustment,
    Group_short,
    N,
    Median,
    Q1,
    Q3
  ) %>%
  
  tidyr::pivot_wider(
    
    names_from =
      Group_short,
    
    values_from =
      c(
        N,
        Median,
        Q1,
        Q3
      ),
    
    names_glue =
      "{.value}_{Group_short}"
  )


# ============================================================
# 14. ADJUSTED MODEL FUNCTION
#
# STANDARD BIOMARKERS:
#
#     log2(Value) ~ Group + Age + Sex
#
# PSA/PZP:
#
#     log2(Value) ~ Group + Age
#
# OMNIBUS GROUP TEST:
#
# drop1(model, test = "F")
#
# This tests whether Group explains significant variation in the
# biomarker after adjustment for the other covariates.
#
# PAIRWISE TESTS:
#
# Estimated marginal means from emmeans:
#
#   NonP vs Control
#   Pneu vs Control
#   Pneu vs NonP
#
# ============================================================

fit_one_biomarker <- function(
    biomarker_name
) {
  
  temp <- long_dat %>%
    
    dplyr::filter(
      Biomarker ==
        .env$biomarker_name
    ) %>%
    
    dplyr::filter(
      !is.na(Log2_Value),
      !is.na(Group),
      !is.na(Age)
    )
  
  
  sex_restriction_value <-
    unique(
      temp$Sex_restriction
    )[1]
  
  
  adjustment_value <-
    unique(
      temp$Adjustment
    )[1]
  
  
  # ----------------------------------------------------------
  # SPECIAL SEX-RESTRICTED BIOMARKERS
  # ----------------------------------------------------------
  
  if (
    biomarker_name %in%
    c(
      "Kallikrein 3/PSA",
      "PZP"
    )
  ) {
    
    model_formula <-
      Log2_Value ~
      Group +
      Age
    
  } else {
    
    # For all remaining biomarkers,
    # require valid sex information.
    
    temp <- temp %>%
      dplyr::filter(
        !is.na(Sex)
      )
    
    
    model_formula <-
      Log2_Value ~
      Group +
      Age +
      Sex
  }
  
  
  # ----------------------------------------------------------
  # CHECK THAT ANALYSIS CAN BE RUN
  # ----------------------------------------------------------
  
  if (
    length(
      unique(
        temp$Group
      )
    ) < 2 ||
    nrow(temp) < 4 ||
    length(
      unique(
        temp$Log2_Value
      )
    ) < 2
  ) {
    
    return(
      
      list(
        
        omnibus =
          data.frame(
            
            Biomarker =
              biomarker_name,
            
            Sex_restriction =
              sex_restriction_value,
            
            Adjustment =
              adjustment_value,
            
            Model_N =
              nrow(temp),
            
            Group_df =
              NA_real_,
            
            Residual_df =
              NA_real_,
            
            Group_F =
              NA_real_,
            
            Group_raw_p =
              NA_real_,
            
            Model_formula =
              paste(
                deparse(
                  model_formula
                ),
                collapse = " "
              ),
            
            stringsAsFactors =
              FALSE
          ),
        
        pairwise =
          NULL,
        
        emmeans =
          NULL
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # FIT MODEL
  # ----------------------------------------------------------
  
  model <- tryCatch(
    
    stats::lm(
      formula =
        model_formula,
      
      data =
        temp
    ),
    
    error =
      function(e)
        NULL
  )
  
  
  if (
    is.null(model)
  ) {
    
    return(
      
      list(
        
        omnibus =
          data.frame(
            
            Biomarker =
              biomarker_name,
            
            Sex_restriction =
              sex_restriction_value,
            
            Adjustment =
              adjustment_value,
            
            Model_N =
              nrow(temp),
            
            Group_df =
              NA_real_,
            
            Residual_df =
              NA_real_,
            
            Group_F =
              NA_real_,
            
            Group_raw_p =
              NA_real_,
            
            Model_formula =
              paste(
                deparse(
                  model_formula
                ),
                collapse = " "
              ),
            
            stringsAsFactors =
              FALSE
          ),
        
        pairwise =
          NULL,
        
        emmeans =
          NULL
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # OMNIBUS GROUP TEST
  #
  # Partial F test:
  # Does Group improve the model after Age / Sex are included?
  # ----------------------------------------------------------
  
  drop_table <- tryCatch(
    
    stats::drop1(
      model,
      test = "F"
    ),
    
    error =
      function(e)
        NULL
  )
  
  
  if (
    !is.null(
      drop_table
    ) &&
    "Group" %in%
    rownames(
      drop_table
    )
  ) {
    
    group_df <-
      as.numeric(
        drop_table[
          "Group",
          "Df"
        ]
      )
    
    
    group_F <-
      as.numeric(
        drop_table[
          "Group",
          "F value"
        ]
      )
    
    
    group_p <-
      as.numeric(
        drop_table[
          "Group",
          "Pr(>F)"
        ]
      )
    
  } else {
    
    group_df <-
      NA_real_
    
    group_F <-
      NA_real_
    
    group_p <-
      NA_real_
  }
  
  
  omnibus_row <-
    data.frame(
      
      Biomarker =
        biomarker_name,
      
      Sex_restriction =
        sex_restriction_value,
      
      Adjustment =
        adjustment_value,
      
      Model_N =
        stats::nobs(
          model
        ),
      
      Group_df =
        group_df,
      
      Residual_df =
        stats::df.residual(
          model
        ),
      
      Group_F =
        group_F,
      
      Group_raw_p =
        group_p,
      
      Model_formula =
        paste(
          deparse(
            stats::formula(
              model
            )
          ),
          collapse = " "
        ),
      
      stringsAsFactors =
        FALSE
    )
  
  
  # ----------------------------------------------------------
  # ESTIMATED MARGINAL MEANS
  # ----------------------------------------------------------
  
  emm <- tryCatch(
    
    emmeans::emmeans(
      model,
      specs =
        ~ Group
    ),
    
    error =
      function(e)
        NULL
  )
  
  
  if (
    is.null(
      emm
    )
  ) {
    
    return(
      
      list(
        
        omnibus =
          omnibus_row,
        
        pairwise =
          NULL,
        
        emmeans =
          NULL
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # DEFINE PAIRWISE CONTRASTS
  #
  # Group order:
  #
  # 1. Control
  # 2. ICU Non-pneumonia
  # 3. ICU Pneumonia
  #
  # Therefore:
  #
  # NonP - Control = -1, +1, 0
  # Pneu - Control = -1, 0, +1
  # Pneu - NonP    = 0, -1, +1
  # ----------------------------------------------------------
  
  contrast_list <-
    list(
      
      "NonP vs Control" =
        c(
          -1,
          1,
          0
        ),
      
      "Pneu vs Control" =
        c(
          -1,
          0,
          1
        ),
      
      "Pneu vs NonP" =
        c(
          0,
          -1,
          1
        )
    )
  
  
  cont <- tryCatch(
    
    emmeans::contrast(
      
      emm,
      
      method =
        contrast_list,
      
      adjust =
        "none"
    ),
    
    error =
      function(e)
        NULL
  )
  
  
  # ----------------------------------------------------------
  # PAIRWISE RESULTS
  # ----------------------------------------------------------
  
  if (
    is.null(
      cont
    )
  ) {
    
    pairwise_rows <-
      NULL
    
  } else {
    
    cont_df <-
      as.data.frame(
        
        summary(
          
          cont,
          
          infer =
            c(
              TRUE,
              TRUE
            ),
          
          level =
            0.95,
          
          adjust =
            "none"
        )
      )
    
    
    pairwise_rows <-
      cont_df %>%
      
      dplyr::transmute(
        
        Biomarker =
          biomarker_name,
        
        Sex_restriction =
          sex_restriction_value,
        
        Adjustment =
          adjustment_value,
        
        Comparison =
          as.character(
            contrast
          ),
        
        Adjusted_log2_difference =
          estimate,
        
        SE =
          SE,
        
        df =
          df,
        
        t_ratio =
          t.ratio,
        
        Raw_p =
          p.value,
        
        Log2_CI_low =
          lower.CL,
        
        Log2_CI_high =
          upper.CL,
        
        Adjusted_fold_change =
          2^estimate,
        
        Adjusted_FC_CI_low =
          2^lower.CL,
        
        Adjusted_FC_CI_high =
          2^upper.CL,
        
        Direction =
          dplyr::case_when(
            
            estimate > 0 ~
              "First group higher",
            
            estimate < 0 ~
              "First group lower",
            
            estimate == 0 ~
              "Equal adjusted means",
            
            TRUE ~
              NA_character_
          )
      )
  }
  
  
  # ----------------------------------------------------------
  # ADJUSTED GROUP MEANS
  # ----------------------------------------------------------
  
  emm_df <-
    as.data.frame(
      
      summary(
        
        emm,
        
        infer =
          c(
            TRUE,
            TRUE
          ),
        
        level =
          0.95
      )
    ) %>%
    
    dplyr::transmute(
      
      Biomarker =
        biomarker_name,
      
      Sex_restriction =
        sex_restriction_value,
      
      Adjustment =
        adjustment_value,
      
      Group =
        as.character(
          Group
        ),
      
      Adjusted_log2_mean =
        emmean,
      
      SE =
        SE,
      
      df =
        df,
      
      Log2_CI_low =
        lower.CL,
      
      Log2_CI_high =
        upper.CL,
      
      Adjusted_geometric_mean =
        2^emmean,
      
      Geometric_mean_CI_low =
        2^lower.CL,
      
      Geometric_mean_CI_high =
        2^upper.CL
    )
  
  
  # ----------------------------------------------------------
  # RETURN RESULTS
  # ----------------------------------------------------------
  
  list(
    
    omnibus =
      omnibus_row,
    
    pairwise =
      pairwise_rows,
    
    emmeans =
      emm_df
  )
}


# ============================================================
# 15. RUN ALL ADJUSTED MODELS
# ============================================================

cat(
  "\nRunning adjusted ANCOVA/regression models...\n"
)


model_results <- lapply(
  
  biomarker_vars,
  
  fit_one_biomarker
)


omnibus_results <-
  dplyr::bind_rows(
    
    lapply(
      model_results,
      `[[`,
      "omnibus"
    )
  )


pairwise_results <-
  dplyr::bind_rows(
    
    lapply(
      model_results,
      `[[`,
      "pairwise"
    )
  )


adjusted_means <-
  dplyr::bind_rows(
    
    lapply(
      model_results,
      `[[`,
      "emmeans"
    )
  )


# ============================================================
# 16. BH-FDR FOR OMNIBUS ADJUSTED GROUP TESTS
#
# One FDR family across all biomarkers
# ============================================================

omnibus_results <- omnibus_results %>%
  
  dplyr::mutate(
    
    Omnibus_Group_FDR =
      stats::p.adjust(
        Group_raw_p,
        method = "BH"
      ),
    
    Omnibus_Group_FDR_significant =
      !is.na(
        Omnibus_Group_FDR
      ) &
      Omnibus_Group_FDR <
      0.05
  ) %>%
  
  dplyr::arrange(
    Omnibus_Group_FDR,
    Group_raw_p
  )


# ============================================================
# 17. BH-FDR FOR PAIRWISE COMPARISONS
#
# Same FDR families as prior script:
#
# FAMILY 1:
#
#   NonP vs Control
#   Pneu vs Control
#
# pooled together
#
# FAMILY 2:
#
#   Pneu vs NonP
#
# separate correction
# ============================================================

icu_vs_control <-
  pairwise_results %>%
  
  dplyr::filter(
    
    Comparison %in%
      c(
        "NonP vs Control",
        "Pneu vs Control"
      )
  ) %>%
  
  dplyr::mutate(
    
    Pairwise_FDR =
      stats::p.adjust(
        Raw_p,
        method = "BH"
      )
  )


pneu_vs_nonp <-
  pairwise_results %>%
  
  dplyr::filter(
    
    Comparison ==
      "Pneu vs NonP"
  ) %>%
  
  dplyr::mutate(
    
    Pairwise_FDR =
      stats::p.adjust(
        Raw_p,
        method = "BH"
      )
  )


pairwise_results_fdr <-
  dplyr::bind_rows(
    
    icu_vs_control,
    
    pneu_vs_nonp
  ) %>%
  
  dplyr::mutate(
    
    Pairwise_FDR_significant =
      !is.na(
        Pairwise_FDR
      ) &
      Pairwise_FDR <
      0.05
  )


# ============================================================
# 18. JOIN RAW DESCRIPTIVES + OMNIBUS RESULTS
# ============================================================

pairwise_results_final <-
  pairwise_results_fdr %>%
  
  dplyr::left_join(
    
    descriptive_wide,
    
    by =
      c(
        "Biomarker",
        "Sex_restriction",
        "Adjustment"
      )
  ) %>%
  
  dplyr::left_join(
    
    omnibus_results %>%
      
      dplyr::select(
        
        Biomarker,
        
        Model_N,
        
        Group_df,
        
        Residual_df,
        
        Group_F,
        
        Group_raw_p,
        
        Omnibus_Group_FDR,
        
        Omnibus_Group_FDR_significant,
        
        Model_formula
      ),
    
    by =
      "Biomarker"
  ) %>%
  
  dplyr::mutate(
    
    Significant =
      Omnibus_Group_FDR_significant &
      Pairwise_FDR_significant
  )


# ============================================================
# 19. ORDER OUTPUT COLUMNS
# ============================================================

pairwise_results_final <-
  pairwise_results_final %>%
  
  dplyr::select(
    
    Biomarker,
    
    Sex_restriction,
    
    Adjustment,
    
    Comparison,
    
    
    N_Control,
    
    Median_Control,
    
    Q1_Control,
    
    Q3_Control,
    
    
    N_NonP,
    
    Median_NonP,
    
    Q1_NonP,
    
    Q3_NonP,
    
    
    N_Pneu,
    
    Median_Pneu,
    
    Q1_Pneu,
    
    Q3_Pneu,
    
    
    Adjusted_log2_difference,
    
    Log2_CI_low,
    
    Log2_CI_high,
    
    Adjusted_fold_change,
    
    Adjusted_FC_CI_low,
    
    Adjusted_FC_CI_high,
    
    Direction,
    
    
    SE,
    
    df,
    
    t_ratio,
    
    Raw_p,
    
    Pairwise_FDR,
    
    Pairwise_FDR_significant,
    
    
    Model_N,
    
    Group_df,
    
    Residual_df,
    
    Group_F,
    
    Group_raw_p,
    
    Omnibus_Group_FDR,
    
    Omnibus_Group_FDR_significant,
    
    
    Significant,
    
    Model_formula
  ) %>%
  
  dplyr::arrange(
    
    dplyr::desc(
      Significant
    ),
    
    Pairwise_FDR,
    
    Omnibus_Group_FDR
  )


# ============================================================
# 20. SIGNIFICANT RESULTS
# ============================================================

significant_results <-
  pairwise_results_final %>%
  
  dplyr::filter(
    Significant
  )


# ============================================================
# 21. COMPARISON-SPECIFIC TABLES
# ============================================================

nonp_vs_control_results <-
  pairwise_results_final %>%
  
  dplyr::filter(
    Comparison ==
      "NonP vs Control"
  )


pneu_vs_control_results <-
  pairwise_results_final %>%
  
  dplyr::filter(
    Comparison ==
      "Pneu vs Control"
  )


pneu_vs_nonp_results <-
  pairwise_results_final %>%
  
  dplyr::filter(
    Comparison ==
      "Pneu vs NonP"
  )


# ============================================================
# 22. BIOMARKER RESTRICTION / ADJUSTMENT TABLE
# ============================================================

restriction_table <-
  data.frame(
    
    Biomarker =
      biomarker_vars,
    
    Sex_restriction =
      dplyr::case_when(
        
        biomarker_vars ==
          "Kallikrein 3/PSA" ~
          "Male only",
        
        biomarker_vars ==
          "PZP" ~
          "Female only",
        
        TRUE ~
          "Both sexes"
      ),
    
    Adjustment =
      dplyr::case_when(
        
        biomarker_vars ==
          "Kallikrein 3/PSA" ~
          "Age adjusted; males only",
        
        biomarker_vars ==
          "PZP" ~
          "Age adjusted; females only",
        
        TRUE ~
          "Age + sex adjusted"
      ),
    
    stringsAsFactors =
      FALSE
  )


# ============================================================
# 23. SPECIAL INDIVIDUAL VALUES
#
# Allows verification that PSA/PZP used the correct participants.
# ============================================================

special_individual_values <-
  long_dat %>%
  
  dplyr::filter(
    
    Biomarker %in%
      special_biomarkers
  ) %>%
  
  dplyr::select(
    
    Deidentified,
    
    Group,
    
    Age,
    
    Sex,
    
    Biomarker,
    
    Sex_restriction,
    
    Adjustment,
    
    Value,
    
    Log2_Value
  ) %>%
  
  dplyr::arrange(
    
    Biomarker,
    
    Group,
    
    Deidentified
  )


# ============================================================
# 24. README
# ============================================================

readme <- data.frame(
  
  Item =
    c(
      
      "Analysis",
      
      "Clinical groups",
      
      "Number of biomarkers",
      
      "Outcome scale",
      
      "Standard biomarker model",
      
      "Kallikrein 3/PSA",
      
      "PZP",
      
      "Omnibus test",
      
      "Omnibus FDR",
      
      "Pairwise comparisons",
      
      "ICU vs Control FDR family",
      
      "Pneumonia vs Non-pneumonia FDR family",
      
      "Significant definition",
      
      "Adjusted fold change",
      
      "Descriptive statistics"
    ),
  
  
  Description =
    c(
      
      paste0(
        "Age- and sex-adjusted three-group biomarker analysis ",
        "using linear regression/ANCOVA on log2 biomarker ",
        "concentrations."
      ),
      
      
      paste(
        group_order,
        collapse = "; "
      ),
      
      
      as.character(
        length(
          biomarker_vars
        )
      ),
      
      
      paste0(
        "All biomarker concentrations are log2-transformed for ",
        "inferential models. Raw concentrations remain in ",
        "descriptive tables."
      ),
      
      
      "log2(Biomarker) ~ Group + Age + Sex.",
      
      
      paste0(
        "Males only; model is log2(PSA) ~ Group + Age. ",
        "Sex is not included because the analyzed sample ",
        "is restricted to males."
      ),
      
      
      paste0(
        "Females only; model is log2(PZP) ~ Group + Age. ",
        "Sex is not included because the analyzed sample ",
        "is restricted to females."
      ),
      
      
      paste0(
        "Partial F test for Group from drop1(model, test='F'), ",
        "adjusted for the other covariates in the model."
      ),
      
      
      paste0(
        "Benjamini-Hochberg correction across all 390 ",
        "omnibus Group P values."
      ),
      
      
      paste0(
        "Estimated marginal mean contrasts: NonP vs Control, ",
        "Pneu vs Control, and Pneu vs NonP."
      ),
      
      
      paste0(
        "NonP vs Control and Pneu vs Control raw pairwise P ",
        "values are corrected together as one BH family."
      ),
      
      
      paste0(
        "Pneu vs NonP raw pairwise P values are corrected ",
        "as a separate BH family."
      ),
      
      
      paste0(
        "Omnibus_Group_FDR < 0.05 AND Pairwise_FDR < 0.05."
      ),
      
      
      paste0(
        "Adjusted fold change = 2^(adjusted difference in log2 means). ",
        "Values >1 mean the first named group is higher; ",
        "values <1 mean the first named group is lower."
      ),
      
      
      paste0(
        "Raw N, median, Q1, and Q3 by clinical group are included ",
        "for descriptive reporting; these are not covariate-adjusted."
      )
    ),
  
  
  stringsAsFactors =
    FALSE
)


# ============================================================
# 25. CREATE EXCEL WORKBOOK
# ============================================================

wb <-
  openxlsx::createWorkbook()


# ============================================================
# 26. EXCEL STYLES
# ============================================================

header_style <-
  openxlsx::createStyle(
    
    fontColour =
      "#FFFFFF",
    
    fgFill =
      "#4472C4",
    
    textDecoration =
      "bold",
    
    halign =
      "center",
    
    valign =
      "center",
    
    border =
      "Bottom"
  )


significant_style <-
  openxlsx::createStyle(
    
    fgFill =
      "#E2F0D9",
    
    fontColour =
      "#006100"
  )


# ============================================================
# 27. HELPER FUNCTION FOR EXCEL SHEETS
# ============================================================

add_sheet <- function(
    sheet_name,
    data
) {
  
  
  openxlsx::addWorksheet(
    
    wb,
    
    sheet_name
  )
  
  
  openxlsx::writeData(
    
    wb,
    
    sheet =
      sheet_name,
    
    x =
      data,
    
    headerStyle =
      header_style,
    
    withFilter =
      TRUE
  )
  
  
  openxlsx::freezePane(
    
    wb,
    
    sheet =
      sheet_name,
    
    firstRow =
      TRUE
  )
  
  
  if (
    ncol(
      data
    ) > 0
  ) {
    
    openxlsx::setColWidths(
      
      wb,
      
      sheet =
        sheet_name,
      
      cols =
        seq_len(
          ncol(
            data
          )
        ),
      
      widths =
        "auto"
    )
  }
}


# ============================================================
# 28. ADD EXCEL SHEETS
# ============================================================

add_sheet(
  "Significant",
  significant_results
)


add_sheet(
  "All Pairwise Results",
  pairwise_results_final
)


add_sheet(
  "NonP vs Control",
  nonp_vs_control_results
)


add_sheet(
  "Pneu vs Control",
  pneu_vs_control_results
)


add_sheet(
  "Pneu vs NonP",
  pneu_vs_nonp_results
)


add_sheet(
  "Adjusted Omnibus",
  omnibus_results
)


add_sheet(
  "Adjusted Means",
  adjusted_means
)


add_sheet(
  "Group Descriptives",
  group_descriptives
)


add_sheet(
  "Descriptives Wide",
  descriptive_wide
)


add_sheet(
  "Restrictions",
  restriction_table
)


add_sheet(
  "Special Individual Values",
  special_individual_values
)


add_sheet(
  "Special Sex Check",
  special_check
)


add_sheet(
  "README",
  readme
)


# ============================================================
# 29. HIGHLIGHT SIGNIFICANT RESULTS
# ============================================================

if (
  nrow(
    significant_results
  ) > 0
) {
  
  openxlsx::addStyle(
    
    wb,
    
    sheet =
      "Significant",
    
    style =
      significant_style,
    
    rows =
      2:
      (
        nrow(
          significant_results
        ) + 1
      ),
    
    cols =
      seq_len(
        ncol(
          significant_results
        )
      ),
    
    gridExpand =
      TRUE
  )
}


# ============================================================
# 30. SAVE WORKBOOK
# ============================================================

openxlsx::saveWorkbook(
  
  wb,
  
  file =
    output_file,
  
  overwrite =
    TRUE
)


# ============================================================
# 31. CONSOLE SUMMARY
# ============================================================

cat(
  "\n\n============================================\n",
  "AGE/SEX-ADJUSTED BIOMARKER ANALYSIS COMPLETE\n",
  "============================================\n\n",
  sep = ""
)


cat(
  "Biomarkers analyzed: ",
  length(
    biomarker_vars
  ),
  "\n",
  sep = ""
)


cat(
  "Standard biomarkers: adjusted for AGE + SEX\n"
)


cat(
  "Kallikrein 3/PSA: MALES ONLY, adjusted for AGE\n"
)


cat(
  "PZP: FEMALES ONLY, adjusted for AGE\n\n"
)


cat(
  "Omnibus Group FDR-significant biomarkers: ",
  sum(
    omnibus_results$Omnibus_Group_FDR_significant,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Significant adjusted pairwise rows: ",
  nrow(
    significant_results
  ),
  "\n\n",
  sep = ""
)


cat(
  "Output saved here:\n",
  output_file,
  "\n",
  sep = ""
)