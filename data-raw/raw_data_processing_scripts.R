


# Set here path
here::here("README.md")

# Import true diagnostic imprintome designs
design_scores = list()
design_scores$Design_Score_Summary <- read.csv(here::here("data-raw", "TruDx_AllDesigns_Design_Score_Summary.csv"))
design_scores$Pass_Score_Threshold <- read.csv(here::here("data-raw", "TruDx_AllDesigns_Pass_Score_Threshold.csv"))
design_scores$Fail_Score_Threshold <- read.csv(here::here("data-raw", "TruDx_AllDesigns_Fail_Score_Threshold.csv"))
design_scores$Failed_Designs <- read.csv(here::here("data-raw", "TruDx_AllDesigns_Failed_Designs.csv"))

usethis::use_data(design_scores)

save(design_scores, file = "design_scores.rda")
