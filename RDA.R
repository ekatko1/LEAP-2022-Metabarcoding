library(readr)
library(dplyr)

## Read data
df = readr::read_tsv("18S_combined_taxa_aggregated.tsv")
View(df)

## Convert to long format
dfl = df |> tidyr::pivot_longer(cols=where(is.numeric), 
                          names_to = c("Mesocosm", "Time", "Replicate"), 
                          names_sep = "_", 
                          values_to = "Count") |> 
  # Remove non-mesocosm data points for this analysis
  filter(Mesocosm != "Resv", Mesocosm != "Cont")

## Read tratment data
trt = readr::read_csv(file = "Treatments.csv")

## Explore experimental time range
dfl |> pull(Time) |> unique()

# Focus on timepoint 5 for rest of analysis (for now).

## Check which mesocosms we don't have data for
trt |> anti_join(
  dfl |> filter(Time == "T5") |> group_by(Mesocosm) |> summarise()
)

T5 = dfl |> filter(Time == "T5")

## Prepare data to fit into a Site x Species matrix
T5_wide = T5 |> 
  tidyr::pivot_wider(values_from = Count, names_from = Domain:species) 

## Match environmental data frame to fit to species matrix
env_T5 = T5_wide |> select(Mesocosm) |> 
  left_join(trt) |> select(-Mesocosm) |> 
  mutate(Node=as.numeric(Node))

## Create the Site x Species matrix. 
matrix_T5 = T5_wide |> 
  mutate(ID = glue::glue("{Mesocosm}_{Time}_{Replicate})"), .keep="unused", .before=1) |>
  tibble::column_to_rownames("ID") |>
  as.matrix()



## Fit RDA to check if there are significant effects
library(vegan)
## TODO: hellinger transform or some other standardization that makes sense for both environment and species counts
fit = rda(formula = matrix_T5 ~ Movement*NutrientP*Node, data = env_T5)
fit
anova.cca(fit, step = 1000, by="term")


# Type 1 scaling
ordiplot(fit, scaling = 1, type = "text")
# Type 2 scaling
ordiplot(fit, scaling = 2, type = "text")
## TODO: Simplify data to a higher order of biological organization to make these plots readbale
            