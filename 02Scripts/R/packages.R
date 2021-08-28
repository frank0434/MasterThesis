library(data.table)
library(readxl)
library(magrittr)
library(forcats)
library(lubridate)
library(DBI)
library(RSQLite)
library(inspectdf)
library(openxlsx)
library(knitr)
library(ggplot2)
library(cowplot)
# library(zoo)
library(broom)
library(dplyr)
library(tidyr)
library(kableExtra)
library(tabulizer)
# library(bookdown)
library(here)

# Workflow control
library(targets)
library(tarchetypes)
library(future)
library(future.callr)

# Stats
library(hydroGOF)
# Customised package
library(autoapsimx)
library(mcp)
library(tidybayes)
potentialPKGS <- c("autoapsimx","clipr"
,"dplyr     ","drake     "
,"DT        ","ggthemes  "
,"ggvis     ","inspectdf "
,"lubridate ","pool      "
,"purrr     ","readr     "
,"reticulate","RPostgreSQL"
,"shiny     ","shinyjs   "
,"tidyr     ","visNetwork")
