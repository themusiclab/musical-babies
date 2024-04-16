
# note: the dev version of bib2df may be necessary as the CRAN version
# currently has an issue with parsing curly braces
# if you do run into this issue, you can load the github latest version 
# like this:
# devtools::install_github("ropensci/bib2df")


# libraries ---------------------------------------------------------------


library(tidyverse)
library(here)
library(bib2df)


# parameters --------------------------------------------------------------


manuscript_file_path <- here("writing", "manuscript.Rmd")
bib_file_path <- here("writing", "themusiclab_long.bib")
output_file_path <- here("writing", "themusiclab.bib")


# extract cite keys from manuscript ---------------------------------------


# load manuscript
Rmd <- readChar(manuscript_file_path, nchars=1e9)

# extract out all citekeys
cites <- str_extract_all(Rmd, "@[a-zA-Z0-9-]*(?=(\\s)|(;)|(])|(\\.)|($)|(\\\\;))")[[1]] %>%
  unique() %>%
  tibble() %>%
  filter(str_detect(., "\\d|(inpress)")) %>%
  mutate(cite = str_remove(., "@")) %>%
  # WARNING: remove a few outliers (assuming that all citekeys start with capital letter...)
  filter(str_detect(., "[A-Z]")) %>%
  pull(cite)


# filter .bib file to just those keys -------------------------------------


bib_out <- bib2df(here("writing", "themusiclab_full.bib")) %>%
  filter(BIBTEXKEY %in% all_of(cites)) %>%
  select(!contains("."))

# bib_out <- bib_out |>
#   mutate(across(where(is.character), ~ str_remove_all(.x, "\\{|\\}")))

df2bib(bib_out, file = here("writing", "themusiclab.bib"))
