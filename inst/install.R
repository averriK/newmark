file.remove("NAMESPACE") |> suppressWarnings()
# source("data-raw/DATASET data is not .R")
usethis::use_proprietary_license(copyright_holder = "Alejandro Verri Kozlowski")
devtools::check(document = TRUE,cran = TRUE)
remove.packages("dsra") |> suppressWarnings()
devtools::install()

## PUSH MAIN
# remotes::install_github("averriK/dsra")

