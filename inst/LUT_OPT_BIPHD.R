Opt_algorithm = c('OCI_Hu12', 'OCI_Hu12', 'TC2_clean', 'TC2_clean', 'TBA_Gil10', 
                  'BR_Git11', 'OCI_Hu12', 'Gons08', 'TC2_clean', 
                  'TBA_Git11', 'TBA_Git11', 'TC2_turbid', 'TBA_Git11', 
                  'TC2_turbid', 'TBA_Gil10', 'Bloom', 'Bloom')
names(Opt_algorithm) <- paste("Cluster", 1:length(Opt_algorithm))

Remove_algorithm <- list(
  c('TBA_Gil10', 'TBA_Git11', 'TC2_turbid', 'Bloom'),
  c('TBA_Gil10', 'TBA_Git11', 'TC2_clean', 'Bloom'),
  c('BR_Git11', 'TBA_Git11', 'TC2_turbid', 'Bloom'),
  c('TBA_Git11', 'OCI_Hu12', 'TC2_turbid', 'Bloom'),
  c('TBA_Git11', 'OCI_Hu12', 'TC2_turbid', 'Bloom'),
  c('OCI_Hu12', 'Gons08', 'TC2_clean', 'Bloom'),
  c('BR_Git11', 'TBA_Git11', 'TC2_clean', 'Bloom'),
  c('TBA_Git11', 'OCI_Hu12', 'TC2_turbid', 'Bloom'),
  c('BR_Git11', 'TBA_Gil10', 'TC2_turbid', 'Bloom'),
  c('BR_Git11', 'TBA_Gil10', 'OCI_Hu12', 'Bloom'),
  c('OCI_Hu12', 'Gons08', 'TC2_clean', 'Bloom'),
  c('BR_Git11', 'OCI_Hu12', 'TC2_clean', 'Bloom'),
  c('OCI_Hu12', 'Gons08', 'TC2_clean', 'Bloom'),
  c('BR_Git11', 'TBA_Gil10', 'Gons08', 'Bloom'),
  c('BR_Git11', 'OCI_Hu12', 'TC2_clean', 'Bloom'),
  c('TBA_Git11', 'OCI_Hu12', 'TC2_clean', 'TC2_turbid'),
  c('BR_Git11', 'OCI_Hu12', 'Gons08', 'TC2_turbid')
)

LUT_OPT_BIPHD <- list(
  Opt_algorithm = Opt_algorithm,
  Remove_algorithm = Remove_algorithm
)

save(LUT_OPT_BIPHD, file = "./data/LUT_OPT_BIPHD.rda", compression_level = 9)
