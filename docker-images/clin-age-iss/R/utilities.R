suppressPackageStartupMessages(library("synapseClient"))

## Synapse utility functions

get.challenge.files <- function() {
    chal_data_table <- synTableQuery('select id,name from syn9763946')
    chal_data_df <- chal_data_table@values
    chal_data_df
}

