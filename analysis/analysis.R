
library(pacman)
p_load(
  tidyverse,
  here,
  ggplot2,
  lubridate,
  RColorBrewer,
  lme4,
  lmerTest,
  nlme,
  kableExtra,
  broom.mixed,
  see,
  blandr
)

#### Convenience functions (mostly used in reporting in the .rmd) and fancy violing plot function -----

# rounds to 2 or 3 decimal places
r2 <- function(num) {format(round(num, 2), nsmall = 2)}
r0 <- function(num) {format(round(num, 0), nsmall = 2)}
# formats to percentage
p <- function(num) {paste0((100*(num)) %>% r2, "%")}

# fancy split violin plot code
GeomSplitViolin <- ggplot2::ggproto(
  "GeomSplitViolin",
  ggplot2::GeomViolin,
  draw_group = function(self,
                        data,
                        ...,
                        # add the nudge here
                        nudge = 0,
                        draw_quantiles = NULL) {
    data <- transform(data,
                      xminv = x - violinwidth * (x - xmin),
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- dplyr::arrange(transform(data,
                                        x = if (grp %% 2 == 1) xminv else xmaxv),
                              if (grp %% 2 == 1) y else -y)
    newdata <- rbind(newdata[1, ],
                     newdata,
                     newdata[nrow(newdata), ],
                     newdata[1, ])
    newdata[c(1, nrow(newdata)-1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    # now nudge them apart
    newdata$x <- ifelse(newdata$group %% 2 == 1,
                        newdata$x - nudge,
                        newdata$x + nudge)
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data,
                                                           draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)),
                         setdiff(names(data), c("x", "y")),
                         drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin",
                       grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...),
                                      quantile_grob))
    }
    else {
      ggplot2:::ggname("geom_split_violin",
                       ggplot2::GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function(mapping = NULL,
                              data = NULL,
                              stat = "ydensity",
                              position = "identity",
                              # nudge param here
                              nudge = 0,
                              ...,
                              draw_quantiles = NULL,
                              trim = TRUE,
                              scale = "area",
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  ggplot2::layer(data = data,
                 mapping = mapping,
                 stat = stat,
                 geom = GeomSplitViolin,
                 position = position,
                 show.legend = show.legend,
                 inherit.aes = inherit.aes,
                 params = list(trim = trim,
                               scale = scale,
                               # don't forget the nudge
                               nudge = nudge,
                               draw_quantiles = draw_quantiles,
                               na.rm = na.rm,
                               ...))
}


#### create combined EMA dataset "dat" ------

miph_dat <- read_csv(here::here("data","clean_ema_data.csv" )) %>% 
  dplyr::select(-user_id, -timezone)

kiwi_dat <- read_csv(here("data","clean_ema_data_kiwi.csv" )) %>% 
  mutate(miph_id = as.character(miph_id)) %>% 
  dplyr::select(-user_id, -timezone)

# full dataset containing anyone who received ema pings
# "dat_full" contains all participants, even those with low response rates
dat_full <- bind_rows("miph" = miph_dat, "kiwi" = kiwi_dat, .id = "cohort") %>% 
  group_by(miph_id) %>% 
  dplyr::select(cohort:answered_date_local, daynumber,singingPeriod:phaseCompletionRate)

# ids to be excluded (response rate less than 50%)
exclusions <- dat_full %>% 
  filter(total_percentage_answered<.50) %>% 
  pull(miph_id) %>% 
  unique

# "dat" has low responders removed (final sample)
dat <- dat_full %>% 
  filter(!(miph_id %in% exclusions))

ids <- dat$miph_id %>% 
  unique

# mean number of pings answered
median_pings <- dat %>% 
  group_by(miph_id) %>% 
  mutate(n=n()) %>% 
  dplyr::select(miph_id, n) %>% 
  distinct() %>% 
  pull(n) %>% 
  median()

#### create combined demographic dataset "demo" ----

miph_demo <- read_csv(here("data","clean_survey_data.csv" )) %>% 
  filter(miph_id %in% ids) %>% # exclude low responders
  dplyr::select(miph_id:condition, demoStartDate:parentingSure) %>% 
  distinct() %>% 
  # make currentcountry = birthcountry if participant still lives where they were born
  mutate(currentCountry = ifelse(is.na(currentCountry), 
                                 birthCountry,
                                 currentCountry),
         # make numberChildren numeric
         numberChildren = ifelse(numberChildren == "5 or more", 
                                 5, 
                                 as.numeric(numberChildren))) 

kiwi_demo <- read_csv(here("data","clean_survey_data_kiwi.csv" )) %>% 
  filter(miph_id %in% ids) %>% # exclude low responders
  dplyr::select(miph_id, name=parentName, baby=babyName, condition, babyBirthday:parentingSure) %>% 
  distinct() %>% 
  mutate(currentCountry = ifelse(is.na(currentCountry), 
                                 birthCountry,
                                 currentCountry)) %>% 
  mutate(miph_id = as.character(miph_id))

demo <- bind_rows("miph" = miph_demo, "kiwi" = kiwi_demo, .id = "cohort") %>% 
  # merge in study start date from ema data
  left_join(dat[,c("miph_id","startDate")] , by = "miph_id") %>%  
  # calculate age at start of study in days
  mutate(age_at_start = as.numeric(startDate - babyBirthday),
         age_at_start_months = age_at_start / 30) %>% 
  distinct()

#### create combined repeated measures survey dataset "survey" ----

miph_survey <- read_csv(here("data","clean_survey_data.csv" )) %>% 
  dplyr::select(rm:parentingEasier)
kiwi_survey <- read_csv(here("data","clean_survey_data_kiwi.csv" )) %>% 
  dplyr::select(condition:parentingEasier) %>% 
  mutate(miph_id = as.character(miph_id))

survey <- bind_rows("miph"= miph_survey, "kiwi" = kiwi_survey, .id = "cohort") %>% 
  filter(miph_id %in% ids) %>% 
  mutate(childcare = factor(childcare, levels = c("I am not involved with child care (0%)",
                                                  "My partner/other adults, including daycare, are the primary caregivers of my baby and I assist (25%)",
                                                  "Caregiving is equally shared with my partner/other adults, including daycare (50%)",
                                                  "I am the primary caregiver and my partner/other adults (including daycare) assist (75%)",
                                                  "I am the sole primary caregiver (100%)")),
         childcare_num = as.numeric(childcare)) %>% 
  filter(miph_id %in% ids) %>%  # exclude low responders
  rowwise() %>% 
  group_by(miph_id) %>%
  mutate(caregiving_test = ifelse(any(childcare_num == 2), 0, 1)) 



#### store participant demographics for reporting ---- 

demographics <- list()

demographics[["total"]][["primary_caregiver"]] <- survey %>% 
  ungroup() %>% 
  dplyr::select(miph_id, caregiving_test) %>% 
  distinct() %>% 
  summarise(sum = sum(caregiving_test))

# full sample size
demographics[["miph"]][["n"]] <- demo[demo$cohort=="miph",] %>% nrow
demographics[["kiwi"]][["n"]] <- demo[demo$cohort=="kiwi",] %>% nrow
demographics[["total"]][["n"]] <- demo %>% nrow

# parent sex
demographics[["miph"]][["parentSex"]] <- demo[demo$cohort=="miph",]$parentSex %>% table
demographics[["kiwi"]][["parentSex"]] <- demo[demo$cohort=="kiwi",]$parentSex %>% table
demographics[["total"]][["parentSex"]] <- demo$parentSex %>% table

# parent age
demographics[["miph"]][["parentAge"]] <- demo[demo$cohort=="miph",]$parentAge %>%  mean(na.rm=T)
demographics[["kiwi"]][["parentAge"]] <- demo[demo$cohort=="kiwi",]$parentAge %>%  mean(na.rm=T)
demographics[["total"]][["parentAge"]] <- demo$parentAge %>%  mean(na.rm=T)

# current country of residence
demographics[["miph"]][["currentCountry"]] <- demo[demo$cohort=="miph",]$currentCountry %>% table
demographics[["kiwi"]][["currentCountry"]] <- demo[demo$cohort=="kiwi",]$currentCountry %>% table
demographics[["total"]][["currentCountry"]] <- demo$currentCountry %>% table

# parent race
demographics[["miph"]][["race"]] <- demo[demo$cohort=="miph",]$parentRace %>% table
demographics[["kiwi"]][["race"]] <- demo[demo$cohort=="kiwi",]$parentRace %>% table
demographics[["total"]][["race"]] <- demo$parentRace %>% table
demographics[["total"]][["hispanic"]] <- demo$parentHispanicLatino %>% table


# parent education
demographics[["miph"]][["education"]] <- demo[demo$cohort=="miph",]$parentEducation %>% table
demographics[["kiwi"]][["education"]] <- demo[demo$cohort=="kiwi",]$parentEducation %>% table
demographics[["total"]][["education"]] <- demo$parentEducation %>% table

# parent income
demographics[["miph"]][["income"]] <- demo[demo$cohort=="miph",]$parentIncome %>% table
demographics[["kiwi"]][["income"]] <- demo[demo$cohort=="kiwi",]$parentIncome %>% table
demographics[["total"]][["income"]] <- demo$parentIncome %>% table

# parent musical background
demographics[["miph"]][["musicbg"]] <- demo[demo$cohort=="miph",]$parentMusicTraining %>% table
demographics[["kiwi"]][["musicbg"]] <- demo[demo$cohort=="kiwi",]$parentMusicTraining %>% table
demographics[["total"]][["musicbg"]] <- demo$parentMusicTraining %>% table

# baby age
demographics[["miph"]][["babyAge"]][["min"]] <- demo[demo$cohort=="miph",]$age_at_start %>% min 
demographics[["miph"]][["babyAge"]][["max"]] <- demo[demo$cohort=="miph",]$age_at_start %>% max 
demographics[["miph"]][["babyAge"]][["mean"]] <- demo[demo$cohort=="miph",]$age_at_start %>%  mean(na.rm=T) 

demographics[["kiwi"]][["babyAge"]][["min"]] <- demo[demo$cohort=="kiwi",]$age_at_start %>% min 
demographics[["kiwi"]][["babyAge"]][["max"]] <- demo[demo$cohort=="kiwi",]$age_at_start %>% max 
demographics[["kiwi"]][["babyAge"]][["mean"]] <- demo[demo$cohort=="kiwi",]$age_at_start %>%mean(na.rm=T) 

demographics[["total"]][["babyAge"]][["min"]] <- demo$age_at_start %>% min 
demographics[["total"]][["babyAge"]][["max"]] <- demo$age_at_start %>% max 
demographics[["total"]][["babyAge"]][["mean"]] <- demo$age_at_start %>% mean(na.rm=T) 
demographics[["total"]][["babyAge"]][["iqr"]] <- demo$age_at_start %>% IQR(na.rm=T) 

# baby sex
demographics[["miph"]][["babySex"]] <- demo[demo$cohort=="miph",]$babySex %>% table
demographics[["kiwi"]][["babySex"]] <- demo[demo$cohort=="kiwi",]$babySex %>% table
demographics[["total"]][["babySex"]] <- demo$babySex %>% table


# how many babies were born preterm?
demographics[["miph"]][["preterm"]] <- demo[demo$cohort=="miph",]$babyDuedate %>% table
demographics[["kiwi"]][["preterm"]] <- demo[demo$cohort=="kiwi",]$babyDuedate %>% table
demographics[["total"]][["preterm"]] <- demo$babyDuedate %>% table

# expectancy effects
demographics[["total"]][["babyMoodPred"]] <- demo$babyMoodPred %>% table
demographics[["total"]][["parentMoodPred"]] <- demo$caregiverMoodPred %>% table

#### data for table: demographics ----

country_table <- demo %>% 
  mutate(total = n()) %>% 
  group_by(currentCountry, birthCountry) %>% 
  summarise(n = n(),
            p = ((n/total)*100) %>% round(.,1)) %>% 
  # mutate(`Recruitment Wave`= ifelse(cohort == "miph","First", "Second")) %>% 
  ungroup() %>% 
  dplyr::select(currentCountry, birthCountry,n,p) %>% 
  distinct() %>% 
  arrange(currentCountry, desc(n)) 

currentcountry_table <- demo %>% 
  mutate(total = n(),
         characteristic = "Country of residence") %>% 
  group_by(characteristic,currentCountry) %>% 
  summarise(n = n(),
            p = ((n/total)*100) %>% round(.,1)) %>% 
  # mutate(`Recruitment Wave`= ifelse(cohort == "miph","First", "Second")) %>% 
  ungroup() %>% 
  dplyr::select(characteristic, value = currentCountry,n,p) %>% 
  distinct() %>% 
  arrange(desc(n)) 

birthcountry_table <- demo %>% 
  mutate(total = n(),
         characteristic = "Parent's country of birth") %>% 
  group_by(characteristic, birthCountry) %>% 
  summarise(n = n(),
            p = ((n/total)*100) %>% round(.,1)) %>% 
  # mutate(`Recruitment Wave`= ifelse(cohort == "miph","First", "Second")) %>% 
  ungroup() %>% 
  dplyr::select(characteristic, value = birthCountry,n,p) %>% 
  distinct() %>% 
  arrange(desc(n)) 

race_table <- demo %>% 
  mutate(total = n(),
         characteristic = "Parent race/ethnicity") %>% 
  group_by(parentRace,characteristic, total) %>% 
  summarise(freq = n()) %>% 
  mutate(value = case_when(parentRace %in% c("White", "European/NZ European") ~ "White/European/New Zealand European",
                           TRUE ~ parentRace)) %>% 
  group_by(characteristic,value) %>% 
  summarise(n = sum(freq),
            p = ((n/total)*100) %>% round(.,1)) %>% distinct() %>% 
  mutate(value = factor(value, 
                        levels= c("White/European/New Zealand European",
                                  "Asian",
                                  "Black or African American",
                                  "Māori",
                                  "More than one race",
                                  "I'd prefer not to say"))) %>% 
  arrange(value) %>% 
  dplyr::select(characteristic, value, n, p)


education_table <- demo %>% 
  mutate(parentEducation = case_when(parentEducation=="Some college/university (still studying or dropped out)" ~ "Some college/university",
                                     parentEducation=="College/university graduate (4 year)" ~ "College/university graduate",
                                     TRUE ~ parentEducation
  ),
  parentEducation = factor(parentEducation, levels = c("High school or equivalent",
                                                       "Vocational/technical school (2 year)",
                                                       "Some college/university",
                                                       "College/university graduate",
                                                       "Master's degree (MA or equivalent)",
                                                       "Doctoral degree (PhD or equivalent)",
                                                       "Professional degree (MD, JD, etc)")),
  total = n(),
  characteristic = "Parent's highest level of education") %>% 
  group_by(characteristic, parentEducation) %>% 
  summarise(n = n(),
            p = ((n/total)*100) %>% round(.,1)) %>% 
  ungroup() %>% 
  dplyr::select(characteristic, value= parentEducation,n,p) %>% 
  distinct() %>% 
  arrange(value) 

# simplified income table (full version below)
income_table <- demo %>% 
  mutate(parentIncome = case_when(cohort == "kiwi" &  parentIncome=="$40,000 to $49,999" ~ "$20,000 to $29,999",
                                  cohort == "kiwi" &  parentIncome=="$50,000 to $74,999" ~ "$30,000 to $39,999",
                                  cohort == "kiwi" &  parentIncome=="$75,000 to $99,999" ~ "$50,000 to $74,999",
                                  cohort == "kiwi" &  parentIncome=="$100,000 to $150,000" ~ "$75,000 to $99,999",
                                  cohort == "kiwi" &  parentIncome=="Over $150,000" ~ "$100,000 to $150,000",
                                  TRUE ~ parentIncome),
         parentIncome_recoded = case_when(parentIncome %in% c("$40,000 to $49,999",
                                                              "$30,000 to $39,999",
                                                              "$20,000 to $29,999",
                                                              "Under $10,000") ~ "Below $50,000",
                                          TRUE ~ parentIncome),
         parentIncome_recoded = factor(parentIncome_recoded, levels= c("Over $150,000",
                                                                       "$100,000 to $150,000" ,
                                                                       "$75,000 to $99,999",
                                                                       "$50,000 to $74,999",
                                                                       "Below $50,000",
                                                                       "I'd prefer not to say")), 
         total = n(),
         characteristic = "Current household income (USD)") %>% 
  group_by(characteristic, parentIncome_recoded) %>% 
  reframe(n = n(),
          p = ((n/total)*100) %>% round(.,1)) %>% 
  ungroup() %>% 
  dplyr::select(characteristic, value = parentIncome_recoded,n,p) %>% 
  distinct() %>% 
  arrange(value) 

table1 <- bind_rows(currentcountry_table, 
                    birthcountry_table,
                    race_table,
                    education_table,
                    income_table)


#### data for table: musical background ----

music_table <- demo %>% 
  dplyr::select(miph_id, parentMusicTraining,parentMusicTraining_7_TEXT) %>% 
  mutate(`No formal music training`= ifelse(grepl("No formal music training", parentMusicTraining),1,0),
         `Lessons or classes before or during elementary school years`= ifelse(grepl("Lessons or classes before or during elementary school years", parentMusicTraining),1,0),
         `Lessons or classes during middle and high school years`= ifelse(grepl("Lessons or classes during middle and high school years", parentMusicTraining),1,0),
         `Lessons or classes in adulthood`= ifelse(grepl("Lessons or classes in adulthood", parentMusicTraining),1,0),
         `Community-based music groups (e.g., church choir, community ensembles)`= ifelse(grepl("Community-based music groups", parentMusicTraining),1,0),
         `Majored in music`= ifelse(grepl("Majored in music", parentMusicTraining),1,0),
         `Other`= ifelse(grepl("Other", parentMusicTraining),1,0))

simple_music_table <- demo %>% 
  dplyr::select(miph_id, parentMusicTraining,parentMusicTraining_7_TEXT) %>% 
  mutate(`Lessons or classes before or during elementary school years`= ifelse(grepl("Lessons or classes before or during elementary school years", parentMusicTraining),1,0),
         `Lessons or classes during middle and high school years`= ifelse(grepl("Lessons or classes during middle and high school years", parentMusicTraining),1,0),
         `Lessons or classes in adulthood`= ifelse(grepl("Lessons or classes in adulthood", parentMusicTraining),1,0),
         `Community-based music groups (e.g., church choir, community ensembles)`= ifelse(grepl("Community-based music groups", parentMusicTraining),1,0),
         `Majored in music`= ifelse(grepl("Majored in music", parentMusicTraining),10,0),
         `Other`= ifelse(grepl("Other", parentMusicTraining),10,0)) %>% 
  rowwise() %>% 
  mutate(sum = sum(`Lessons or classes before or during elementary school years`,
                   `Lessons or classes during middle and high school years`,
                   `Lessons or classes in adulthood`,
                   `Community-based music groups (e.g., church choir, community ensembles)`,
                   `Majored in music`,
                   `Other`),
         bg = case_when(sum == 0 ~ "None",
                        sum <= 2 ~ "Some",
                        sum <= 4 ~ "Intermediate",
                        TRUE ~ "Advanced"),
         bg = factor(bg, levels = c("Advanced",
                                    "Intermediate",
                                    "Some",
                                    "None"))) %>% 
  dplyr::select(miph_id, bg ) %>% 
  group_by(bg) %>% 
  summarise(n = n(),
            p = ((n/110)*100) %>% r0)
  
         
#### data for table: EMA studies with parent-infant dyads ------

EMAreview <- data.frame(
  Study = c("Allen et al (2018)", "von Stumm & Latham (2018)", "Adams et al (2019); Paul et al (2014)", "Franchak (2019)", "Omowale et al (2022); Mendez et al (2019)", "Suga et al (2022)", "Corpuz et al (2023)", "de Barbaro et al (2023)", "Wenze et al (2023)", "Imhoff et al (2023); Hoffmann et al (2023)", "Franchak et al (2024)", "Malachowski et al (2023)"),
  Population = c("Pregnant women with a recent history of smoking", "Mother-infant dyads", "Primiparous mother-infant dyads", "Caregivers of infants", "Pregnant (and postpartum) women", "Mother-infant dyads", "First time fathers", "Mother-infant dyads", "Postpartum mothers", "Expecting (and postpartum) parents", "Caregivers of infants", "Mother-infant dyads"),
  Infant_Age = c("Followed from birth through 12 weeks", "2-15 weeks", "Followed from birth through 3 years", "3, 6, 9, & 12 months", "Followed from birth through 1 year", "2-8 months", "9-12 months", "Mage = 3.9 months", "6-12 & 18-24 weeks", "Followed from birth through 6 months", "Followed from 10 to 13 months","4-6 months"),
  Sample_Size = c(46, 53, 157, 95, 296, 129, 194, 53, 130, 297, 62, 60),
  Context_of_Assessment = c("Maternal mood and smoking-related symptomatology", "Mother/infant dietary, sleeping patterns, maternal wellbeing", "Infant fussiness and soothing techniques in daily life", "Infant body position", "Maternal weight and cardiometabolic health", "Mother-infant interaction and maternal wellbeing", "Paternal care of infants", "Maternal mental health in relation to infant crying", "Mother-infant bonding, stress, and sleep", "Parents' gendered stereotypes for infants", "Infants' everyday experiences with restraint, body position, and objects", "Infant placements and their exposure to adult speech"),
  Study_Duration = c("16 weeks (36 weeks gestation - 12 weeks postpartum)", "3 weeks", "3 years", "7 days", "Mean of 15 months (18-24 weeks' gestation - 1 year postpartum)", "1 month", "13 months (third trimester - 10 months postpartum)", "7 days", "7 days", "Pregnancy - 6 months postpartum", "4 months", "1 month"),
  Study_Methods = c("Baseline visit + EMA + 5 clinic visits", "Online surveys + EMA", "parenting intervention (4 nurse home visits) + annual clinic visit + EMA", "EMA", "3 lab visits + postpartum assessments + EMA", "Pre-questionnaire + EMA + post-questionnaire", "3 home visits + EMA", "Introductory session + EMA + audio sampling", "Baseline survey + EMA", "Psychological measurements + EMA", "Introductory call + EMA + exit survey", "Introductory session + EMA + audio sampling"),
  EMA_Frequency = c("1/day over 16 weeks", "4/week over 3 weeks", "5/day over 5-8 days at two separate times (3 & 8 weeks)", "5/day over 7 days", "3/day on average over ~15 months", "1-3/day over 1 month", "8/day at 6 separate days over 12 weeks", "6/day over 7 days", "4/day over 7 days", "1/day for 2 weeks then 1/week for 4 weeks", "10/day over 4 days at four separate time points", "3 separate days within a 1-month period"),
  EMA_Items = c("Mood (24) daily, postpartum depression (10) weekly", "Parent and infant wellness (2), mother diet (4), mother sleep (7), infant diet (4), infant sleep (1) 3/week; mother wellbeing (3), support (1) 1/week", "Infant fussiness and soothing techniques used (7)", "Infant positions (3)", "Sleep, diet, physical activity, stress, mood, support, breastfeeding, depression, experience with racism", "Feelings & heart rate measurement", "Interactions at the moment (1)", "Depression (2), anxiety (2), negative affect (5) 1/day; (2); negative affect (5) 6/day", "Sleep-related (14) 1/day; mood (20), stress(1), fatigue(1), % connections with infants(1) 4/day", "Emotional well-being (12), breastfeeding (3), wound healing (6), infant well-being and behavior (10)", "Infants' restraint, body position & object holding (6)", "Infant use of independence-supporting placements (2)"),
  Compliance_Rate = c("80.91% ± 2.28% (postpartum)", "Not specified", "Not specified", "95.60%", "69.5% (1 month postpartum), 66.2% (2 months), 64.2% (3 months), 61.8% (4 months), 62% (5 months)", "Not specified", "91%", "~78% (average 33 out of 42 pings)", "~75% (average of 21.1 out of 28)", "Not specified", "~71% (average of 28.4 out of 40)", "~91.6%"),
  Retention_Rate = c("87%", "93% (7% attrition)", "~93% (completed)", "~97% (completed)", "~15% attrition (terminated data collection before publication)", "~75% (completed)", "~93% (completed)", "~90% (completed)", "Unclear (221 baseline, 130 completed)", "Not specified", "~85% (completed)", "Not specified"),
  Latency = c("Not specified", "EMA expired at midnight", "EMA data collected within 1 hour", "Median = 2.3 min", "Not specified", "Not specified", "EMA data collected within 30 min", "Not specified", "Median = 1.68 min", "Not specified", "Median = 0.05 min", "Not specified")
)



#### can we merge the two cohorts? ----

# empty list to put cohort comparison test results in
cohort_compare <- list()


# compare infant age at start of study in two cohorts
cohort_compare[["age"]] <- wilcox.test(age_at_start ~ cohort, 
                                       data = demo, 
                                       paired = F) %>% broom::tidy()


# data with compliance rates by study phase 
compliance_by_phase_dat <- dat %>% 
  dplyr::select(miph_id, cohort, studyPhase, singingPeriod, total_percentage_answered, phaseCompletionRate) %>% 
  distinct()

compliance_by_group_dat <- dat %>% 
  dplyr::select(miph_id, cohort,condition, studyPhase, singingPeriod, total_percentage_answered, phaseCompletionRate) %>% 
  distinct()

# data with compliance rates across whole study
compliance_total_dat <- dat %>%
  left_join(demo, by=c("miph_id","cohort","condition")) %>% 
  dplyr::select(miph_id, cohort,condition,  total_percentage_answered,age_at_start) %>% 
  distinct() 

# compliance was not related to infant age at the start of the study
cohort_compare[["compliance_by_age"]] <- lm(total_percentage_answered~age_at_start, compliance_total_dat) %>% tidy

cohort_compare[["compliance_by_condition"]] <- lm(total_percentage_answered~condition, compliance_total_dat) %>% tidy
cohort_compare[["compliance_by_condition_pretest"]] <- lm(phaseCompletionRate~condition, compliance_by_group_dat %>% filter(studyPhase == 1)) %>% tidy
cohort_compare[["compliance_by_condition_intervention"]] <- lm(phaseCompletionRate~condition, compliance_by_group_dat %>% filter(studyPhase == 2)) %>% tidy
cohort_compare[["compliance_by_condition_posttest"]] <- lm(phaseCompletionRate~condition, compliance_by_group_dat %>% filter(studyPhase == 3)) %>% tidy
cohort_compare[["compliance_by_condition_7-10"]] <- lm(phaseCompletionRate~condition, compliance_by_group_dat %>% filter(studyPhase == 4)) %>% tidy

#### singing frequency analyses ----

singing_analyses <- list()

# data with weekly proportion of times participant reported having sung/played music/etc in last 2-3 hours
sing_check_dat <-dat %>% 
  filter(EMA_wereYouWithBaby == "Yes" & EMA_sickBaby == "No") %>% 
  group_by(cohort, condition, miph_id, week) %>% 
  summarise(sangPercentage = mean(EMA_sungToBabyLastHour == "Yes", na.rm = TRUE),
            musicForSelfPercentage = mean(EMA_musicForOwnEnjoyment == "Yes", na.rm = TRUE),
            musicRecordedPercentage = mean(EMA_playedRecordedMusic == "Yes", na.rm = TRUE)) %>% 
  pivot_wider(id_cols = c("condition","miph_id","cohort"), names_from = week, values_from = c(sangPercentage, musicForSelfPercentage,musicRecordedPercentage))

# compare conditions 
singing_analyses[["hourly_singing_check"]][["compare_conditions_pretest"]] <- t.test(sing_check_dat[sing_check_dat$condition == "SingFirst",]$sangPercentage_1, 
                                                                                     sing_check_dat[sing_check_dat$condition == "SingSecond",]$sangPercentage_1, 
                                                                                     paired = F) %>% tidy()

singing_analyses[["hourly_singing_check"]][["compare_conditions_posttest"]] <- t.test(sing_check_dat[sing_check_dat$condition == "SingFirst",]$sangPercentage_6, 
                                                                                     sing_check_dat[sing_check_dat$condition == "SingSecond",]$sangPercentage_6, 
                                                                                     paired = F) %>% tidy()


sing_check_dat_long <- dat %>% 
  filter(week <= 6) %>% 
  filter(EMA_wereYouWithBaby == "Yes" & EMA_sickBaby == "No") %>% 
  group_by(cohort, condition, miph_id, daynumber) %>% 
  summarise(sangPercentage = mean(EMA_sungToBabyLastHour == "Yes", na.rm = TRUE) * 100,
            musicForSelfPercentage = mean(EMA_musicForOwnEnjoyment == "Yes", na.rm = TRUE) * 100,
            musicRecordedPercentage = mean(EMA_playedRecordedMusic == "Yes", na.rm = TRUE) * 100) %>% 
  mutate(condition = factor(condition, 
                            levels=c("SingSecond","SingFirst")))
singing_analyses[["hourly_singing_check"]][["autoregression"]] <- lme(sangPercentage ~ 1 + condition*daynumber,
                                                                          random = ~1|miph_id,
                                                                          correlation = corAR1(form = ~ daynumber|miph_id),
                                                                          data = sing_check_dat_long %>% filter) %>% tidy(.,"fixed")

singing_analyses[["hourly_music_check"]][["self"]][["autoregression"]] <- lme(musicForSelfPercentage ~ 1 + condition*daynumber,
                                                                                  random = ~1|miph_id,
                                                                                  correlation = corAR1(form = ~ daynumber|miph_id),
                                                                                  data = sing_check_dat_long)  %>% tidy(.,"fixed")
singing_analyses[["hourly_music_check"]][["recorded"]][["autoregression"]] <- lme(musicRecordedPercentage ~ 1 + condition*daynumber,
                                                                                      random = ~1|miph_id,
                                                                                      correlation = corAR1(form = ~ daynumber|miph_id),
                                                                                      data = sing_check_dat_long) %>% tidy(.,"fixed")


# calculate mean percentage of having sung in past 2-3 hours for each week
for (condition in c("SingFirst","SingSecond")) {
  for (week in c(1:10)) {
    
    var = paste0("sangPercentage_", week)
    
    singing_analyses[["hourly_singing_check"]][[condition]][["mean"]][[week]] <- 
      sing_check_dat[sing_check_dat$condition==condition,][[var]] %>% mean(na.rm=T) 
    singing_analyses[["hourly_singing_check"]][[condition]][["sd"]][[week]] <- 
      sing_check_dat[sing_check_dat$condition==condition,][[var]] %>% sd(na.rm=T) 
    
  }
}


# data with weekly average of singing frequency
sing_freq_dat <-  dat %>% 
  group_by(condition,cohort, week, miph_id) %>% 
  mutate(weeklySing = mean(EMA_babySingHowManyTimes, na.rm=T)) %>% 
  ungroup() %>% 
  dplyr::select(condition, cohort,week, studyPhase, miph_id, weeklySing) %>% 
  distinct() %>% 
  pivot_wider(id_cols = c("condition","miph_id","cohort"), names_from = week, values_from = weeklySing) %>% 
  mutate(d = `6`- `1`)


# calculate mean singing frequency for each week
for (condition in c("SingFirst","SingSecond", "total")) {
  for (week in c(1:10)) {
    
    var = as.character(week)
    
    if (condition %in% c("SingFirst","SingSecond")){ 
      singing_analyses[["singing_freq"]][[condition]][[week]] <- 
        sing_freq_dat[sing_freq_dat$condition==condition,][[var]] %>% mean(na.rm=T) 
    } else {
      singing_analyses[["singing_freq"]][[condition]][[week]] <- 
        sing_freq_dat[[var]] %>% mean(na.rm=T) 
    }
    
  }
}



# compare weekly singing frequency between conditions

singing_analyses[["singing_freq"]][["compare_conditions_pretest"]] <- t.test(sing_freq_dat[sing_freq_dat$condition == "SingFirst",]$`1`, 
                                                                             sing_freq_dat[sing_freq_dat$condition == "SingSecond",]$`1`, 
                                                                             paired = F)


singing_analyses[["singing_freq"]][["compare_conditions_posttest"]] <- t.test(sing_freq_dat[sing_freq_dat$condition == "SingFirst",]$`6`, 
                                                                              sing_freq_dat[sing_freq_dat$condition == "SingSecond",]$`6`, 
                                                                              paired = F)


sing_freq_dat_long <- dat %>% 
  filter(week <= 6) %>% 
  group_by(condition, cohort,daynumber, miph_id) %>% 
  mutate(dailySing = mean(EMA_babySingHowManyTimes, na.rm=T)) %>% 
  ungroup() %>% 
  dplyr::select(condition,cohort, daynumber, studyPhase, miph_id, dailySing) %>% 
  distinct() %>% 
  drop_na(dailySing) %>% 
  mutate(condition = factor(condition, 
                            levels=c("SingSecond","SingFirst")))


singing_analyses[["singing_freq"]][["autoregression "]] <- lme(dailySing ~ 1 + condition*daynumber,
                                                               random = ~1|miph_id,
                                                               correlation = corAR1(form = ~ daynumber|miph_id),
                                                               data = sing_freq_dat_long) %>% tidy(.,"fixed")

lme(dailySing ~ 1 + condition*daynumber,
    random = ~1|miph_id,
    correlation = corAR1(form = ~ daynumber|miph_id),
    data = sing_freq_dat_long %>% filter(cohort == "kiwi")) %>% summary


#### soothing methods analyses ------

soothing_analyses <- list()


soothing_method_dat_long <- dat %>% 
  filter(week<=6) %>% 
  dplyr::select(miph_id, condition, cohort, week, daynumber, studyPhase,EMA_sickBaby, EMA_babyFussy, starts_with("EMA_babyHowCalmDown.")) %>% 
  filter(EMA_babyFussy == "Yes"  & EMA_sickBaby == "No") %>% 
  dplyr::select(-EMA_babyFussy, -EMA_sickBaby, -studyPhase, -daynumber) %>% 
  group_by(condition,cohort,week, miph_id) %>% 
  mutate(total = n())  %>% 
  group_by(condition,cohort,week, miph_id, total) %>% 
  summarise_all(sum, na.rm = TRUE) %>% 
  ungroup()  %>% 
  mutate(across(starts_with("EMA_babyHowCalmDown."), ~./total)) %>% 
  ungroup() 

soothing_method_dat <- soothing_method_dat_long%>% 
  group_by(condition, week) %>% 
  summarise(across(starts_with("EMA_babyHowCalmDown."), mean, na.rm = TRUE)) %>% 
  pivot_longer(cols = starts_with("EMA_babyHowCalmDown."), names_to = "method", values_to = "proportion") %>% 
  mutate(method = str_replace_all(method, "EMA_babyHowCalmDown.", ""))

# soothing behaviours in order of frequency (overall)
soothing_method_dat_long %>% 
  summarise(across(starts_with("EMA_babyHowCalmDown."), mean, na.rm = TRUE)) %>% 
  pivot_longer(cols = starts_with("EMA_babyHowCalmDown."), names_to = "method", values_to = "proportion") %>% 
  mutate(method = str_replace_all(method, "EMA_babyHowCalmDown.", "")) %>% 
  arrange(desc(proportion))

soothing_plot1_data <- soothing_method_dat %>% 
  filter(condition == "SingFirst" & week <= 6) %>%
  filter(method %in% c("Rub.or.pat", "Sing","Feed","Bounce.rock.walk.swing", "Shush.white.noise","Play.recorded.music")) %>% 
  mutate(method = case_when( method == "Rub.or.pat" ~ "Rub or pat",
                             method == "Bounce.rock.walk.swing" ~ "Bounce/Rock/Walk/Swing",
                             method == "Shush.white.noise" ~ "Shush or use white noise",
                             method == "Play.recorded.music" ~ "Play recorded music",
                             TRUE ~ method))

soothing_plot2_data <- soothing_method_dat_long %>% 
  filter(week <=6) %>% 
  mutate(phase = case_when(week ==1 ~ "Pre-test\n",
                           week %in% c(2:5) ~ "Intervention \n(avg.)",
                           week == 6 ~ "Post-test\n"),
         phase = factor(phase, levels = c("Pre-test\n", "Intervention \n(avg.)", "Post-test\n")),
         EMA_babyHowCalmDown.Sing = EMA_babyHowCalmDown.Sing) %>%
  dplyr::select(miph_id,condition, phase, EMA_babyHowCalmDown.Sing) %>% 
  group_by(condition, phase) %>% 
  summarise(mean = mean(EMA_babyHowCalmDown.Sing, na.rm = TRUE),
            sd = sd(EMA_babyHowCalmDown.Sing, na.rm = TRUE),
            n = n(),
            se = sd/sqrt(n)) 
  
soothing_plot2_data_weekly <- soothing_method_dat_long %>% 
  filter(week <=6) %>% 
  mutate(EMA_babyHowCalmDown.Sing = EMA_babyHowCalmDown.Sing*100) %>%
  dplyr::select(miph_id,condition, week, EMA_babyHowCalmDown.Sing) %>% 
  group_by(condition, week) %>% 
  summarise(mean = mean(EMA_babyHowCalmDown.Sing, na.rm = TRUE),
            sd = sd(EMA_babyHowCalmDown.Sing, na.rm = TRUE),
            n = n(),
            se = sd/sqrt(n)) 
  


soothing_analyses[["singing_change"]][["manipulation"]] <- t.test(soothing_method_dat_long[soothing_method_dat_long$condition == "SingFirst" & soothing_method_dat_long$week == "6",]$`EMA_babyHowCalmDown.Sing`, 
                                                                  soothing_method_dat_long[soothing_method_dat_long$condition == "SingFirst" & soothing_method_dat_long$week == "1",]$`EMA_babyHowCalmDown.Sing`
                                                                  ) %>%  tidy

soothing_analyses[["singing_change"]][["control"]] <- t.test(soothing_method_dat_long[soothing_method_dat_long$condition == "SingSecond" & soothing_method_dat_long$week == "6",]$`EMA_babyHowCalmDown.Sing`, 
                                                                  soothing_method_dat_long[soothing_method_dat_long$condition == "SingSecond" & soothing_method_dat_long$week == "1",]$`EMA_babyHowCalmDown.Sing`
                                                             ) %>%  tidy

soothing_analyses[["singing_change"]][["compare_groups"]] <- t.test(soothing_method_dat_long[soothing_method_dat_long$condition == "SingFirst" & soothing_method_dat_long$week == "6",]$`EMA_babyHowCalmDown.Sing`, 
                                                                    soothing_method_dat_long[soothing_method_dat_long$condition == "SingSecond"  & soothing_method_dat_long$week == "6",]$`EMA_babyHowCalmDown.Sing`, 
                                                                    paired = F) %>% tidy


soothing_analyses[["recorded_music_change"]][["compare_groups"]] <- t.test(soothing_method_dat_long[soothing_method_dat_long$condition == "SingFirst" & soothing_method_dat_long$week == "6",]$`EMA_babyHowCalmDown.Play.recorded.music`, 
                                                                          soothing_method_dat_long[soothing_method_dat_long$condition == "SingSecond" & soothing_method_dat_long$week == "6",]$`EMA_babyHowCalmDown.Play.recorded.music`)


##### baby mood analyses ----

# z-score mood measures within participant
dat <- dat %>% 
  group_by(miph_id) %>% 
  mutate(across(contains(c("Positive", "Passive", "Energetic")), 
                scale, 
                .names = "{.col}_z")) %>% 
  ungroup()

baby_mood_dat <- dat %>% 
  # drop sick days and empty data points 
  filter(EMA_wereYouWithBaby=="Yes" & EMA_sickBaby == "No") %>% 
  group_by(condition, cohort,week, miph_id) %>% 
  mutate(weeklyPos = mean(EMA_babyMoodPositiveNegative, na.rm=T),
         weeklyEne = mean(EMA_babyMoodPassiveAroused, na.rm=T),
         weeklyPos_scaled = mean(EMA_babyMoodPositiveNegative_z, na.rm=T),
         weeklyEne_scaled = mean(EMA_babyMoodPassiveAroused_z, na.rm=T),
         weeklyParentMood_scaled = mean(EMA_adultHowPositive_z, na.rm=T),
  ) %>% 
  ungroup() %>% 
  dplyr::select(condition, cohort,week, miph_id, starts_with("weekly")) %>% 
  distinct()

baby_mood_dat_wide <- baby_mood_dat %>% 
  pivot_wider(id_cols = c("condition","cohort", "miph_id"), names_from = week, values_from = starts_with("weekly")) 

baby_tests <- list()


# across conditions: compare week 1 to week 6 
baby_tests[["mood"]][["all_1_6"]] <- t.test(baby_mood_dat_wide$weeklyPos_scaled_6,
                                                baby_mood_dat_wide$weeklyPos_scaled_1,
                                                alternative = "greater",
                                                paired=T)


# between condition: compare manipulation and control groups at pre-test/week 1 (nonsig)
baby_tests[["mood"]][["compare_pretest"]] <- t.test(baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingFirst", ]$weeklyPos_scaled_1,
                                                    baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingSecond", ]$weeklyPos_scaled_1)

# between condition: compare manipulation and control groups at post-test/week 6 (sig)
baby_tests[["mood"]][["compare_posttest"]] <- t.test(baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingFirst", ]$weeklyPos_scaled_6,
                                                     baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingSecond", ]$weeklyPos_scaled_6,
                                                     alternative = "greater")


# calculate mean baby mood score for each week
for (condition in c("SingFirst","SingSecond")) {
  for (week in c(1:10)) {
    
    var = paste0("weeklyPos_", week)
    
    baby_tests[["mood"]][[condition]][[week]] <- 
      baby_mood_dat_wide[baby_mood_dat_wide$condition==condition,][[var]] %>% mean(na.rm=T) 
    
  }
}


# replicate results with mixed model with auto regression
mm_dat <- dat %>% 
  filter(week <= 6 & EMA_sickBaby == "No") %>% 
  drop_na(EMA_babyMoodPositiveNegative) %>% # drop missing data
  dplyr::select(miph_id, cohort,condition, week, daynumber, EMA_babyMoodPositiveNegative, EMA_adultHowPositive, EMA_babyMoodPassiveAroused, EMA_sickBaby, EMA_sungToBabyLastHour, EMA_playedRecordedMusic, EMA_musicForOwnEnjoyment) %>% 
  group_by(miph_id, daynumber, condition,cohort) %>% 
  summarise(mood = mean(EMA_babyMoodPositiveNegative),
            energy = mean(EMA_babyMoodPassiveAroused),
            parent_mood = mean(EMA_adultHowPositive)) %>% 
  mutate(condition = factor(condition,
                            levels= c("SingSecond","SingFirst")))

# significant condition*day interaction: manipulation babies' mood increased over time
baby_tests[["mood"]][["autoregression"]] <- lme(mood ~ 1 + condition*daynumber,
                                                random = ~1|miph_id,
                                                correlation = corAR1(form = ~ daynumber|miph_id),
                                                data = mm_dat) %>% tidy(.,"fixed")

baby_tests[["mood"]][["sd_increase"]] <- (1.56)/sd(mm_dat$mood) #0.08


baby_mood_dat_plot <- dat %>% 
  dplyr::select(condition,miph_id,week, starts_with("EMA_babyMood")) %>% 
  filter(week %in% c(1,6)) %>% 
  mutate(week = factor(ifelse(week == 1, "Pre-test", "Post-test"), levels = c("Pre-test", "Post-test"))) %>% 
  drop_na() %>% 
  group_by(condition,week, miph_id) %>%
  summarise(weeklyPos=mean(EMA_babyMoodPositiveNegative),
            weeklyPos_scaled=mean(EMA_babyMoodPositiveNegative_z),
            weeklyEne=mean(EMA_babyMoodPassiveAroused))

##### parent mood analyses ----

# nothing particularly interesting here
parent_mood_dat <- dat %>% 
  filter(EMA_wereYouWithBaby=="Yes" & EMA_sickBaby == "No") %>% 
  group_by(condition, week, miph_id) %>% 
  mutate(weeklyPos = mean(EMA_adultHowPositive, na.rm=T),
         weeklyEne = mean(EMA_adultHowEnergetic, na.rm=T),
         weeklyPos_scaled = mean(EMA_adultHowPositive_z, na.rm=T)) %>% 
  ungroup() %>% 
  dplyr::select(condition, week, miph_id, starts_with("weekly")) %>% 
  distinct()

parent_mood_dat_wide <- parent_mood_dat %>% 
  pivot_wider(id_cols = c("condition","miph_id"), names_from = week, values_from = starts_with("weekly")) 

parent_tests <- list()


parent_tests[["compare_baseline_mood_scaled"]] <-  t.test(parent_mood_dat_wide[parent_mood_dat_wide$condition=="SingFirst", ]$weeklyPos_scaled_1,
                                                          parent_mood_dat_wide[parent_mood_dat_wide$condition=="SingSecond", ]$weeklyPos_scaled_1)

parent_tests[["compare_posttest_mood_scaled"]] <-  t.test(parent_mood_dat_wide[parent_mood_dat_wide$condition=="SingFirst", ]$weeklyPos_scaled_6,
                                                          parent_mood_dat_wide[parent_mood_dat_wide$condition=="SingSecond", ]$weeklyPos_scaled_6,
                                                          alternative = "greater")



parent_tests[["mood"]][["autoregression"]] <- lme(parent_mood ~ 1 + condition*daynumber,
                                                random = ~1|miph_id,
                                                correlation = corAR1(form = ~ daynumber|miph_id),
                                                data = mm_dat) %>% tidy(.,"fixed")

# adding parent mood to the condition*daynumber model removes the interaction effect

baby_tests[["mood"]][["autoregression_with_parents"]] <- lme(mood ~ 1 + condition*daynumber*parent_mood,
    random = ~1|miph_id,
    correlation = corAR1(form = ~ daynumber|miph_id),
    data = mm_dat) %>% tidy(.,"fixed")


parent_mood_dat_plot <- dat %>% 
  dplyr::select(condition,miph_id,week, EMA_adultHowPositive_z) %>% 
  filter(week %in% c(1,6)) %>% 
  mutate(week = factor(ifelse(week == 1, "Pre-test", "Post-test"), levels = c("Pre-test", "Post-test"))) %>% 
  drop_na() %>% 
  group_by(condition, week, miph_id) %>%
  summarise(weeklyPos_scaled=mean(EMA_adultHowPositive_z))

##### were parent and baby mood correlated?

# yes, correlated
parent_tests[["baby_parent_corr"]] <- cor.test(dat$EMA_babyMoodPositiveNegative_z, dat$EMA_adultHowPositive_z, 
                                               method = "spearman") %>% tidy()

parent_tests[["baby_parent_corr_compare"]] <- lm(EMA_babyMoodPositiveNegative_z ~ EMA_adultHowPositive_z*condition, dat %>% filter(daynumber<42)) %>% tidy

ggplot(dat, aes(x = EMA_babyMoodPositiveNegative_z, y = EMA_adultHowPositive_z)) + 
  geom_point(alpha = .2)+
  geom_smooth()

# how to prove that these two measures are tapping different constructs?

## parent and infant mood differentially predict the parent's social life

parent_v_baby <- dat %>% 
  filter(EMA_wereYouWithBaby=="Yes" & EMA_sickBaby == "No") %>% 
  group_by(condition, daynumber, miph_id) %>%
  mutate(dailyParentingStress_num = as.numeric(factor(EMA_adultHowStressful, 
                                                      levels = c("Not stressful at all",
                                                                 "A bit stressful",
                                                                 "Somewhat stressful",
                                                                 "Very stressful"
                                                      )))) %>% 
  summarise(dailyParentingStress_num=mean(dailyParentingStress_num),
         EMA_adultHowSocial = mean(EMA_adultHowSocial),
         EMA_adultHowPositive_z = mean(EMA_adultHowPositive_z),
         EMA_babyMoodPositiveNegative_z = mean(EMA_babyMoodPositiveNegative_z),
         EMA_adultHowPositive = mean(EMA_adultHowPositive),
         EMA_babyMoodPositiveNegative = mean(EMA_babyMoodPositiveNegative)
         ) 
  
parent_tests[["baby_parent_social"]] <-  lme(EMA_adultHowSocial ~ 1 + EMA_babyMoodPositiveNegative*EMA_adultHowPositive,
                                             random = ~1|miph_id,
                                             correlation = corAR1(form = ~ daynumber|miph_id),
                                             data = parent_v_baby) %>% tidy(.,"fixed")

parent_tests[["baby_parent_stressful"]] <-  lme(dailyParentingStress_num ~ 1 + EMA_babyMoodPositiveNegative*EMA_adultHowPositive,
                                             random = ~1|miph_id,
                                             correlation = corAR1(form = ~ daynumber|miph_id),
                                             data = parent_v_baby) %>% tidy(.,"fixed")


 ggplot(parent_v_baby, aes(x=EMA_adultHowPositive_z , y = EMA_adultHowSocial,col = miph_id)) + 
  theme_classic() +
  geom_smooth(method = "lm",se = F, alpha =.2, size = .1, span=.75) +
  theme(legend.position = "none") +
  labs(x ="Parent: Z-scored Mood", y = "Parent: Connectedness")

ggplot(parent_v_baby, aes(x=EMA_babyMoodPositiveNegative_z , y = EMA_adultHowSocial, col=miph_id)) + 
  theme_classic() +
  geom_smooth(method = "lm",se = F, alpha =.2, size = .1) +
  theme(legend.position = "none") +
  labs(x = "Baby: Z-scored Mood" , y = "Parent: Connectedness")
 

## parent and infant mood differentially predict the parent's stress levels

ggplot(parent_v_baby, aes(x=EMA_adultHowPositive_z , y = dailyParentingStress_num, col=miph_id)) + 
  geom_smooth(method = "lm", se = F, alpha =.2, size = .1) +
  theme_classic() +
  theme(legend.position = "none")+
  labs(x ="Parent: Z-scored Mood", y = "Parenting Stress") +
  xlim(-6,4)

ggplot(parent_v_baby, aes(x=EMA_babyMoodPositiveNegative_z , y = dailyParentingStress_num, col=miph_id)) + 
  geom_smooth(method = "lm", se = F, alpha =.2, size = .1) +
  theme_classic() +
  theme(legend.position = "none")+
  labs(x ="Baby: Z-scored Mood", y = "Parenting Stress")+
  xlim(-6,4)


#### latency analyses -----

latency_data <- kiwi_dat %>% 
  filter(miph_id %in% ids) %>% 
  left_join(demo, by=c("miph_id","condition")) %>% 
  distinct() %>% 
  mutate(week = case_when(daynumber < 7 ~ 1,
                          daynumber < 14 ~ 2,
                          daynumber < 21 ~ 3,
                          daynumber < 28 ~ 4,
                          daynumber < 35 ~ 5,
                          daynumber < 42 ~ 6,
                          daynumber < 49 ~ 7,
                          daynumber < 56 ~ 8,
                          daynumber < 63 ~ 9,
                          daynumber < 70 ~10)) %>%
  dplyr::select(miph_id, age_at_start, daynumber, week, pingtimestamp, answered_timestamp_utc) %>% 
  mutate(current_age = age_at_start + daynumber,
         pingtimestamp_num = force_tz(as.POSIXct(pingtimestamp, format = "%a, %d %b %Y %H:%M:%S GMT")),
         answered_timestamp_utc =force_tz(answered_timestamp_utc),
         latency = difftime(answered_timestamp_utc,pingtimestamp_num, unit="mins")) 
  

latency <- median(latency_data$latency)
hifreq_latency <- median(latency_data[latency_data$week %in% c(1,6),]$latency)
lowfreq_latency <- median(latency_data[!(latency_data$week %in% c(1,6)),]$latency)

latency_model <- lme(mean_latency ~ 1 + current_age,
                     random = ~1|miph_id,
                     correlation = corAR1(form = ~ daynumber|miph_id),
                     data = latency_data %>% group_by(miph_id, daynumber, current_age) %>% 
                       summarise(mean_latency = mean(latency))) %>% 
  tidy(.,"fixed")

#### replicate mood results in miph ----

baby_tests[["mood"]][["miph"]][["compare_pretest"]] <- t.test(baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingFirst" & baby_mood_dat_wide$cohort=="miph", ]$weeklyPos_scaled_1,
                                                    baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingSecond" & baby_mood_dat_wide$cohort=="miph", ]$weeklyPos_scaled_1)%>% tidy()

# between condition: compare manipulation and control groups at post-test/week 6 (sig)
baby_tests[["mood"]][["miph"]][["compare_posttest"]] <- t.test(baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingFirst" & baby_mood_dat_wide$cohort=="miph", ]$weeklyPos_scaled_6,
                                                     baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingSecond" & baby_mood_dat_wide$cohort=="miph", ]$weeklyPos_scaled_6,
                                                     alternative = "greater")%>% tidy()


baby_tests[["mood"]][["miph"]][["autoregression"]] <- lme(mood ~ 1 + condition*daynumber,
    random = ~1|miph_id,
    correlation = corAR1(form = ~ daynumber|miph_id),
    data = mm_dat %>% filter(cohort == "miph")) %>% tidy(.,fixed=T)

#### replicate mood results in kiwi ----


baby_tests[["mood"]][["kiwi"]][["compare_pretest"]] <- t.test(baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingFirst" & baby_mood_dat_wide$cohort=="kiwi", ]$weeklyPos_scaled_1,
                                                              baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingSecond" & baby_mood_dat_wide$cohort=="kiwi", ]$weeklyPos_scaled_1) %>% tidy()


baby_tests[["mood"]][["kiwi"]][["compare_posttest"]] <- t.test(baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingFirst" & baby_mood_dat_wide$cohort=="kiwi", ]$weeklyPos_scaled_6,
                                                               baby_mood_dat_wide[baby_mood_dat_wide$condition=="SingSecond" & baby_mood_dat_wide$cohort=="kiwi", ]$weeklyPos_scaled_6,
                                                               alternative = "greater") %>% tidy()

baby_tests[["mood"]][["kiwi"]][["autoregression"]] <- lme(mood ~ 1 + condition*daynumber,
    random = ~1|miph_id,
    correlation = corAR1(form = ~ daynumber|miph_id),
    data = mm_dat %>% filter(cohort == "kiwi")) %>% tidy(.,fixed=T)

#### missing data -----

missing_dat <- dat %>% 
  dplyr::select(miph_id, daynumber) %>% 
  mutate(responded = 1) %>% 
  group_by(miph_id,daynumber) %>% 
  mutate(id = row_number())

length(ids)

missing_dat2 <- tibble(miph_id = rep(ids, each=98), 
                       daynumber = rep(c(rep(c(0:6), each = 3), c(7:34), rep(c(35:41), each = 3), c(42:69)), 110)) %>% 
  group_by(miph_id, daynumber) %>% 
  mutate(id = row_number())

tmp <- missing_dat %>% 
  full_join(missing_dat2, by= c("miph_id", "daynumber","id")) %>% 
  mutate(responded = replace_na(responded, 0)) %>% 
  arrange(miph_id, daynumber, id)

#### save it alllll ----

save(list = ls(), file = here("results", "analysis.RData"))

