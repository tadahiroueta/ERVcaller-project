library(vcfR)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(readxl)
library(ggplot2)
library(tidytext)
library(reshape2)
library(data.table)
library(gridExtra)



AUSSIE_PATIENT_DATA_PATH <- paste0("Y:\\Lucas\\ERVcaller-project\\",
                                   "documentation\\aussie-file-tracker.xlsx")
BAYLOR_PATIENT_DATA_PATH <- paste0("Y:\\Lucas\\ERVcaller-project\\",
                                   "documentation\\baylor-patient-data.xlsx")
AUSSIE_OUTPUT_PATH <- paste0("Y:\\Lucas\\ERVcaller-project\\",
                             "ERVcaller-output\\aussie\\all-chr\\")
BAYLOR_OUTPUT_PATH <- paste0("Y:\\Lucas\\ERVcaller-project\\",
                             "ERVcaller-output\\baylor\\")
ERV_REFERENCE_PATH <- "Y:\\Lucas\\ERVcaller-project\\ERVs.xlsx"
DATA_FRAMES_PATH <- "Y:\\Lucas\\ERVcaller-project\\data-frames\\"

#' Parses vcf file
#' 
#' @param file_number Number id
#' @return A data frame with fixed data from vcf
read_aussie_vcf <- function(file_number) {
  path <- paste0(AUSSIE_OUTPUT_PATH, "SRR", as.character(file_number),
                 ".vcf")
  tryCatch(
    as.data.frame(read.vcfR(path)@fix),
    error = function(e) { NULL }
  )
}

#' Formats ERVcaller output file
#' 
#' @param insertions A data frame read from a vcf file
#' @param erv_reference df with ERV names
#' @param erv_only whether you should filter for HERVs
#' 
#' @return Data frame
#' \itemize{
#'  \item `chromosome` (character): chr_
#'  \item `position` (integer): start position of insertion
#'  \item `reference` (character): nucleotide
#'  \item `alternative` (character): insertion type
#'  \item `info` (character): misc
#'  \item `erv_name` (character): name of ERV
#' }
format_insertions <- function(insertions, erv_reference, erv_only = FALSE) {
  # file hasn't been made
  if (is.null(insertions)) {
    return(NULL)
  }
  
  insertions %>%
    select(CHROM, POS, REF, ALT, INFO) %>%
    rename(chromosome = CHROM,
           position = POS,
           reference = REF,
           alternative = ALT,
           info = INFO) %>%
    mutate(alternative = gsub("<INS_MEI:|>$", "", alternative)) %>%
    filter(!erv_only | alternative %in% c("HERV", "SVA"))
}

#' Finds insertions present in `insertion_1` but not in `insertion_2`
#' 
#' @param insertion_1 where insertions should be present
#' @param insertion_2 where insertions should not be present
#' 
#' @return list of insertions
exclusive_comparison <- function(insertion_1, insertion_2) {
  list(anti_join(insertion_1, insertion_2, by = c("chromosome", "position")))
}

#' Finds insertions present in both `insertion_1` and `insertion_2`
#' 
#' @param insertion_1 where insertions should be present
#' @param insertion_2 where insertions should also be present
#' 
#' @return list of insertions
overlap_comparison <- function(insertion_1, insertion_2) {
  list(inner_join(insertion_1, insertion_2,
                  by = c("chromosome", "position")) %>%
         select(chromosome, position, reference.x, alternative.x, info.x, info.y) %>%
         rename(reference = reference.x, alternative = alternative.x))
}

#' Finds insertions present in either `insertion_1` and `insertion_2`
#' 
#' @param insertion_1 where insertions could be present
#' @param insertion_2 where insertions could be present
#' 
#' @return list of insertions
inclusive_comparison <- function(insertion_1, insertion_2) {
  list(full_join(insertion_1, insertion_2, by = c("chromosome", "position")) %>%
         mutate(
           reference = case_when(!is.na(reference.x) ~ reference.x,
                                 !is.na(reference.y) ~ reference.y),
           alternative = case_when(!is.na(reference.x) ~ alternative.x,
                                   !is.na(reference.y) ~ alternative.y)) %>%
         select(chromosome, position, reference, alternative, info.x, info.y))
}

#' Fetches main ERVcaller output data from Aussie sample
#' 
#' @param .aussie_patient_data metadata about patients
#' @param erv_reference data frame with erv_names
#' @param erv_only whether you should filter for HERVs
#' 
#' @return data frame
#' \itemize{
#'  \item `id` (character): patient ID
#'  \item `race` (character): African or European
#'  \item `gleason` (integer): Gleason score
#'  \item `tumour_purity` (numeric): Tumour purity percentage
#'  \item `tumour_mutation_burden` (numeric): Tumour mutation burden
#'  \item `blood_insertions` (data frame): insertions found in blood
#'  \item `tumour_insertions` (data frame): insertions found in tumour
#'  \item `blood_exclusive_insertions` (data frame): insertions found in blood
#'  but not in tumour
#'  \item `tumour_exclusive_insertions` (data frame): insertions found in tumour
#'  but not in blood
#'  \item `overlapping_insertions` (data frame): insertions found both in blood
#'  and tumour
#'  \item `all_insertions` (data frame): insertions found in either blood or
#'  tumour
#'  }
fetch_aussie_patients <- function(.aussie_patient_data, erv_reference,
                                  erv_only = FALSE) {
  .aussie_patient_data %>%
    rowwise() %>%
    mutate(blood_insertions = list(format_insertions(
      read_aussie_vcf(blood_file_number), erv_reference, erv_only)),
      
      tumour_insertions = list(format_insertions(
        read_aussie_vcf(tumour_file_number), erv_reference, erv_only))) %>%
    
    # for missing VCFs
    filter(!is.null(blood_insertions) & !is.null(tumour_insertions)) %>%
    select(id, race, gleason, tumour_purity, tumour_mutation_burden,
           blood_insertions, tumour_insertions) %>%
    
    mutate(exclusive_blood_insertions =
             exclusive_comparison(blood_insertions, tumour_insertions),
           
           exclusive_tumour_insertions =
             exclusive_comparison(tumour_insertions, blood_insertions),
           
           overlapping_insertions =
             overlap_comparison(blood_insertions, tumour_insertions),
           
           all_insertions =
             inclusive_comparison(blood_insertions, tumour_insertions))
}

#' Fetches main ERVcaller output data from Baylor sample
#' 
#' @param .baylor_patient_data metadata about patients
#' @param erv_reference data frame with erv names
#' @param erv_only whether you should filter for HERVs
#' 
#' @return data frame
#' \itemize{
#'  \item `id` (character): patient ID
#'  \item `race` (character): African or European
#'  \item `gleason` (integer): Gleason score
#'  \item `tumour_stage` (character): Tumour date pT___
#'  \item `treatment` (vector): P({ surgery, radiation, ADT })
#'  \item `age` (integer): Age
#'  \item `time` (integer): In months
#'  \item `insertions` (data frame): insertions found
#'  }
fetch_baylor_patients <- function(.baylor_patient_data, erv_reference,
                                  erv_only = FALSE) {
  .baylor_patient_data %>%
    rowwise() %>%
    mutate(insertions = list(format_insertions(tryCatch(
      as.data.frame(read.vcfR(paste0(BAYLOR_OUTPUT_PATH, id, ".vcf"))@fix)
      , error = function(e) { NULL }), erv_reference, erv_only))) %>%
    filter(!is.null(insertions))
}

#' Generates grid of every insertion and whether each patient had it
#' 
#' @param aussie_patients data frame with sample data
#' @param insertions column with insertion data
#' 
#' @return data frame grid
#' \itemize{
#'  \item `chromosome` (character): chr_
#'  \item `position` (integer): start position of insertion
#'  \item `reference` (character): nucleotide
#'  \item `alternative` (character): insertion type
#'  \item `+` (1 | NA) one column for each patient; whether patient has it
#' }
get_grid <- function(aussie_patients, insertions) {
  grid <- data.frame(chromosome = character(),
                     position = character(),
                     reference = character(),
                     alternative = character())
  for (i in 1:nrow(aussie_patients)) {
    patient_id <- make.names(aussie_patients$id[i])
    grid <- full_join(grid, insertions[[i]],
                      by = c("chromosome", "position")) %>%
      mutate(!!patient_id := ifelse(!is.na(reference.y), 1, NA)) %>%
      mutate(reference = ifelse(is.na(reference.x), reference.y, reference.x),
             alternative = ifelse(is.na(alternative.x),
                                  alternative.y, alternative.x)) %>%
      select(-reference.x, -reference.y, -alternative.x, -alternative.y, -info)
  }
  grid <- grid %>%
    mutate(position = as.numeric(position)) %>%
    arrange(chromosome, position)
}

#' Generates table with whether each patient has each insertion
#' 
#' @param grid from get_grid
#' @param patient_data data frame without insertion data
#' @return data frame table
get_table <- function(grid, patient_data) {
  grid %>%
    pivot_longer(cols = -c(chromosome, position, reference, alternative),
                 names_to = "patient", values_to = "insertion") %>%
    mutate(patient = gsub("^X", "", patient)) %>%
    left_join(patient_data, by = c("patient" = "id")) %>%
    select(-intersect(names(grid), c("blood_file_number", 
                                     "tumour_file_number")))
}

#' Calculate and compare frequency of a insertion among races
#' 
#' @param table data frame with long tallies of insertions
#' @param insertion_value value of insertion column
#' @param patients data frame with main sample data
#' 
#' @return data frame grid
#' \itemize{
#'  \item `chromosome` (character): chr_
#'  \item `position` (integer): start position of insertion
#'  \item `reference` (character): nucleotide
#'  \item `alternative` (character): insertion type
#'  \item `african_count` (integer): count of African patients with
#'  the insertion
#'  \item `european_count` (integer): count of European patients
#'  with the insertion
#'  \item `african_proportion` (decimal): proportion of African patients
#'  with the insertion
#'  \item `european_proportion` (decimal): proportion of European
#'  patients with the insertion
#'  \item `proportion_difference` (decimal): difference between African
#'  and European proportions
#' }
count_insertions <- function(table, insertion_value, patients) {
  table %>%
    group_by(chromosome, position, reference, alternative, race) %>%
    summarise(count = sum(insertion == insertion_value, na.rm = TRUE)) %>%
    group_by(chromosome, position, reference, alternative) %>%
    summarise(african_count = first(count),
              european_count = last(count),
              .groups = "drop") %>%
    mutate(total_count = african_count + european_count,
           total_proportion = total_count / patients %>% nrow(),
           african_proportion = african_count / patients %>%
             filter(race == "African") %>% nrow(),
           european_proportion = european_count / patients %>%
             filter(race == "European") %>% nrow(),
           proportion_difference = african_proportion - european_proportion) %>%
    arrange(desc(proportion_difference))
}

#' Generates QQ plot
#' 
#' @param patients sample data frame to test
#' @param target_race race whose data will be tested
#' @param insertion_name column name of data to be tested
check_qq_plot <- function(patients, target_race, insertion_name) {
  patients_with_counts <- patients %>%
    mutate(count = nrow(!!ensym(insertion_name)))
  
  ggplot(patients_with_counts %>% filter(race == target_race),
         aes(sample = count)) +
    stat_qq() +
    stat_qq_line() +
    labs(title = paste(target_race, "QQ Plot"),
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal()
}

#' Generates box plot comparing insertion counts between races
#' 
#' @param patients sample data frame with insertion data
#' @param insertion_name type of insertion to count
#' @param title title for the graph
#' @return plot
make_box_plot <- function(patients, insertion_name, title) {
  patients_with_counts <- patients %>%
    mutate(count = nrow(!!ensym(insertion_name)))
  
  t_test <- t.test((patients_with_counts %>%
                      filter(race == "African"))$count,
                   (patients_with_counts %>%
                      filter(race == "European"))$count)
  print(t_test)
  
  ggplot(patients_with_counts, aes(x = race, y = count, fill = race)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .25) +
    annotate("text", x = .75,
             y = min((patients_with_counts %>%
                        filter(race == "African"))$count) - 5,
             label = paste("p =", round(t_test$p.value, digits = 6))) +
    scale_fill_manual(values = c("African" = "green", "European" = "orange")) +
    labs(x = "Race", y = "Insertion count", title = title) +
    guides(fill = "none") +
    theme_minimal() +
    theme(axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 1),
          axis.ticks.length = unit(7.5, "pt"))
}

#' Generates a line graph for proportions of each race's frequencies of a
#' insertion above a minimum proportion for a chromosome
#' 
#' @param insertion_count data frame with focus on insertions
#' @param chromosome_value chromosome to act as filter
#' @param minimum_proportion minimum proportion for either races 
#' @return plot
plot_proportions <- function(insertion_count, chromosome_value,
                             minimum_proportion) {
  long <- insertion_count %>%
    filter(chromosome == chromosome_value,
           african_proportion > minimum_proportion |
             european_proportion > minimum_proportion) %>%
    pivot_longer(cols = c(african_proportion, european_proportion),
                 names_to = "race", values_to = "proportion") %>%
    mutate(race = recode(race, "african_proportion" = "African",
                         "european_proportion" = "European"))
  
  ggplot(long, aes(x = position, y = proportion, colour = race)) +
    geom_line(linewidth = 0.1) +
    scale_colour_manual(values = c("African" = "green",
                                   "European" = "orange")) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(title = paste0("Heart Monitor, ", chromosome_value),
         x = "Insertion Start Position",
         y = "Proportion of Patients w/ Insertion")
}

#' Generates an area chart with the difference in proportions of frequency of
#' a insertion among races for a chromosome
#' 
#' @param insertion_count data frame with focus on insertions
#' @param chromosome_value chromosome to act as filter
#' @return plot
plot_difference_in_proportions <- function(insertion_count, chromosome_value) {
  ggplot(insertion_count %>% filter(chromosome == chromosome_value), 
         aes(x = position, y = proportion_difference)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Insertion Start Position",
         y = "Difference In Proportions")
}
#
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH

#
#' Generates a gradient tile visualization for proportions of each race's frequencies of
#' insertions above a minimum proportion for a chromosome
#' 
#' @param insertion_count data frame with focus on insertions
#' @param chromosome_value chromosome to act as filter
#' @param minimum_proportion minimum proportion for either races
#' @param length optional chromosome length for scaling the position axis
#' @param scale_min minimum value for color scale (default: 0)
#' @param scale_max maximum value for color scale (default: 1)
#' @return plot grid with both populations side by side
plot_proportions_gradient <- function(insertion_count, 
                                      chromosome_value,
                                      minimum_proportion,
                                      length = NULL,
                                      scale_min = 0,
                                      scale_max = 1,
                                      bin_size = NULL,
                                      debug = FALSE) {
  
  # Separate filtering for African and European proportions
  filtered_african <- insertion_count %>%
    filter(chromosome == chromosome_value,
           african_proportion > minimum_proportion)
  
  filtered_european <- insertion_count %>%
    filter(chromosome == chromosome_value,
           european_proportion > minimum_proportion)
  
  # Debug info
  if(debug) {
    cat("Filtered African data dimensions:", nrow(filtered_african), "rows x", ncol(filtered_african), "columns\n")
    cat("Filtered European data dimensions:", nrow(filtered_european), "rows x", ncol(filtered_european), "columns\n")
  }
  
  # Set y-axis limit
  y_min <- 0
  y_max <- ifelse(!is.null(length), length, max(insertion_count$position, na.rm = TRUE) * 1.01)
  
  # Determine bin size if not provided
  if(is.null(bin_size)) {
    range <- max(insertion_count$position, na.rm = TRUE) - min(insertion_count$position, na.rm = TRUE)
    bin_size <- ceiling(range / 1000)
    if(debug) cat("Auto-calculated bin size:", bin_size, "\n")
  }
  
  # Create binned data for African proportion
  binned_african <- filtered_african %>%
    mutate(position_bin = floor(position / bin_size) * bin_size) %>%
    group_by(position_bin) %>%
    summarize(
      african_proportion_mean = mean(african_proportion, na.rm = TRUE),
      count = n()
    ) %>%
    ungroup()
  
  # Create binned data for European proportion
  binned_european <- filtered_european %>%
    mutate(position_bin = floor(position / bin_size) * bin_size) %>%
    group_by(position_bin) %>%
    summarize(
      european_proportion_mean = mean(european_proportion, na.rm = TRUE),
      count = n()
    ) %>%
    ungroup()
  
  # Ensure values are within range to avoid warnings
  binned_african <- binned_african %>%
    mutate(african_proportion_mean = pmin(pmax(african_proportion_mean, scale_min), scale_max))
  
  binned_european <- binned_european %>%
    mutate(european_proportion_mean = pmin(pmax(european_proportion_mean, scale_min), scale_max))
  
  transform_color <- function(x) sqrt(x)

  
  # Create plot for African proportion
  p1 <- ggplot(binned_african, 
               aes(x = 0, y = position_bin, fill = transform_color(african_proportion_mean), height = bin_size)) +
    geom_tile(width = 0.8) +
    scale_fill_gradient(low = "white", high = "red", 
                        limits = c(transform_color(0), transform_color(0.7)),
                        na.value = "grey90") +  
    scale_y_continuous(limits = c(y_min, y_max)) +
    scale_x_continuous(breaks = 0, labels = "African") +
    labs(title = paste0("African Proportion, ", chromosome_value),
         x = "",
         y = "Position (Mb)",
         fill = "Proportion") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  # Create plot for European proportion  
  p2 <- ggplot(binned_european, 
               aes(x = 0, y = position_bin, fill = transform_color(european_proportion_mean), height = bin_size)) +
    geom_tile(width = 0.8) +
    scale_fill_gradient(low = "white", high = "blue", 
                        limits = c(transform_color(0), transform_color(0.7)),
                        na.value = "grey90") +  
    scale_y_continuous(limits = c(y_min, y_max),
                       labels = function(x) x/1000000) +  
    scale_x_continuous(breaks = 0, labels = "European") +
    labs(title = paste0("European Proportion, ", chromosome_value),
         x = "",
         y = "Position (Mb)",
         fill = "Proportion") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  # Combine plots
  grid.arrange(p1, p2, ncol = 2,
               top = grid::textGrob(paste0("Insertion Proportions, ", chromosome_value),
                                    gp = grid::gpar(fontsize = 14, font = 2)))
}

#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH
#     THIS IS A WORK IN PROGRESS PLEASE DONT BREAK INSHALLAH


# metadata
.aussie_patient_data <- read_xlsx(AUSSIE_PATIENT_DATA_PATH,
                                  range = "A3:G18") %>%
  rename(id = `Patient ID`,
         race = Race,
         gleason = Gleason,
         tumour_purity = `Tumour purity`,
         tumour_mutation_burden = `Tumour mutation burden`,
         blood_file_number = Blood,
         tumour_file_number = Tumour)

.baylor_patient_data <- read_xlsx(BAYLOR_PATIENT_DATA_PATH, 
                                  sheet = "clinical_data") %>%
  rename(id = `Assigned Name`,
         race = `RACE/ETHNICITY`,
         tumour_stage = `Tumor stage`,
         gleason_primary = `Gleason_primary`,
         gleason_secondary = `Gleason_secondary`,
         age = Age,
         time = time_rp_lastFU) %>%
  mutate(race = gsub("WA", "European", race),
         race = gsub("BA", "African", race),
         gleason = gleason_primary + gleason_secondary,
         treatment = gsub(" only", "", treatment),
         treatment = strsplit(treatment, ", "),
         treatment = trimws(treatment),
         time = gsub("months", "", time),
         time = gsub("month", "", time),
         time = trimws(time),
         time = as.integer(time)) %>%
  select(id, race, tumour_stage, gleason, treatment, age, time)

.erv_reference <- read_xlsx(
  ERV_REFERENCE_PATH, sheet = "Sheet1",
  col_types = c("text", "numeric", "numeric", "text")) %>%
  rename(chromosome = chr,
         erv_name = gene) %>%
  mutate(chromosome = paste0("chr", chromosome))

# main samples
aussie_patients <- fetch_aussie_patients(.aussie_patient_data, .erv_reference)
aussie_patients_herv <- fetch_aussie_patients(.aussie_patient_data,
                                              .erv_reference, erv_only = TRUE)

# Baylor samples showed almost no HERV insertions
baylor_patients <- fetch_baylor_patients(.baylor_patient_data, .erv_reference)

# insertion focus
.blood_grid <- get_grid(aussie_patients, aussie_patients$blood_insertions)
.tumour_grid <- get_grid(aussie_patients, aussie_patients$tumour_insertions)
.baylor_grid <- get_grid(baylor_patients, baylor_patients$insertions)

# long table to make plots
.blood_table <- get_table(.blood_grid, .aussie_patient_data)
.tumour_table <- get_table(.tumour_grid, .aussie_patient_data)
.baylor_table <- get_table(.baylor_grid, .baylor_patient_data)

.aussie_table <- full_join(
  .blood_table, .tumour_table,
  by = c("chromosome", "position", "patient", "race"),
  suffix = c(".blood", ".tumour")) %>%
  
  mutate(
    insertion = case_when(
      !is.na(insertion.blood) & !is.na(insertion.tumour) ~ "Tumour & Blood", 
      !is.na(insertion.blood) & is.na(insertion.tumour) ~ "Blood",
      is.na(insertion.blood) & !is.na(insertion.tumour) ~ "Tumour",
      TRUE ~ " "),
    reference = case_when(!is.na(reference.blood) ~ reference.blood,
                          !is.na(reference.tumour) ~ reference.tumour),
    alternative = case_when(
      !is.na(alternative.blood) ~ alternative.blood,
      !is.na(alternative.tumour) ~ alternative.tumour)) %>% 
  
  select(chromosome, position, patient, race,
         insertion, reference, alternative) %>%
  
  arrange(alternative, chromosome, position)

# frequency comparison for insertions
aussie_insertions <- count_insertions(.aussie_table, "Tumour",
                                      aussie_patients)
aussie_insertions_erv <- aussie_insertions %>%
  filter(alternative %in% c("HERV", "SVA"))
baylor_insertions <- count_insertions(.baylor_table, 1, baylor_patients)

# exporting tables
write_rds(aussie_patients, paste0(DATA_FRAMES_PATH, "aussie_patients.rds"))
write_rds(aussie_patients_herv, paste0(DATA_FRAMES_PATH,
                                       "aussie_patients_herv.rds"))
write_rds(baylor_patients, paste0(DATA_FRAMES_PATH, "baylor_patients.rds"))
write_rds(.blood_grid, paste0(DATA_FRAMES_PATH, "blood_grid.rds"))
write_rds(.tumour_grid, paste0(DATA_FRAMES_PATH, "tumour_grid.rds"))
write_rds(.baylor_grid, paste0(DATA_FRAMES_PATH, "baylor_grid.rds"))
write_rds(aussie_insertions, paste0(DATA_FRAMES_PATH,
                                    "aussie_insertions.rds"))
write_rds(aussie_insertions_erv, paste0(DATA_FRAMES_PATH,
                                        "aussie_insertions_erv.rds"))
write_rds(baylor_insertions, paste0(DATA_FRAMES_PATH,
                                    "baylor_insertions.rds"))

# Q-Q plot to check for normality - Aussie checks out... Baylor, not so much
# (400 x 300)
check_qq_plot(aussie_patients, "African", "exclusive_tumour_insertions")
check_qq_plot(aussie_patients, "European", "exclusive_tumour_insertions")
check_qq_plot(aussie_patients_herv, "African", "tumour_insertions")
check_qq_plot(aussie_patients_herv, "European", "tumour_insertions")
check_qq_plot(baylor_patients, "African", "insertions")
check_qq_plot(baylor_patients, "European", "insertions")


# box plot of ERV insertions over race (540 x 540)
make_box_plot(aussie_patients, "tumour_insertions", "Aussie, Tumour Insertions")
make_box_plot(aussie_patients, "exclusive_tumour_insertions",
              "Aussie, Exclusive Tumour Insertions")
make_box_plot(aussie_patients, "all_insertions", "Aussie, All Insertions")
make_box_plot(aussie_patients_herv, "tumour_insertions", "Aussie, Tumour HERV")
make_box_plot(baylor_patients, "insertions", "Baylor, All Insertions")

# TODO try with HERV
# grid of genes over Aussie patients (whether on blood and tumour)
# ggplot(.aussie_table %>% mutate(position = as.character(position)),
#       aes(patient, position, fill = insertion)) +
#  geom_tile() +
#  scale_fill_manual(values = c("Blood" = "red", 
#                               "Tumour" = "black",
#                               "Tumour & Blood" = "darkred",
#                               " " = "white")) +
#  theme_minimal() +
#  theme(axis.text.x = element_text(size = 6))

# heart monitor
plot_proportions(aussie_insertions, 'chr1', 0.2)

# chromo
plot_difference_in_proportions(aussie_insertions, "chr1")

# attempt on bar chromosomal position
plot_proportions_gradient(aussie_insertions, "chr1", 0.4, bin_size = 2500000)
