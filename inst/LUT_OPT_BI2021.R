Opt_algorithm = c('TC2_clean', 'BR_Gil10', 'Gons08', 'TC2_clean', 'TC2_clean', 'BR_Git11', 'TBA_Git11', 
                  'TBA_Gil10', 'TC2_turbid', 'TBA_Gil10', 'TC2_turbid', 'TC2_turbid', 'Gons08', 'Bloom')
names(Opt_algorithm) <- paste("Cluster", 1:length(Opt_algorithm))

Remove_algorithm <- list(
  c('Bloom', 'BR_Git11', 'TBA_Git11', 'TC2_turbid'),
  c('Bloom', 'TBA_Git11', 'Gons08', 'TC2_turbid'),
  c('TBA_Gil10', 'Bloom', 'TBA_Git11', 'TC2_turbid'),
  c('Bloom', 'TBA_Git11', 'Gons08', 'TC2_turbid'),
  c('BR_Gil10', 'TBA_Gil10', 'Bloom', 'BR_Git11'),
  c('TBA_Gil10', 'Bloom', 'TC2_clean', 'TC2_turbid'),
  c('BR_Gil10', 'TBA_Gil10', 'Bloom', 'BR_Git11'),
  c('BR_Gil10', 'Bloom', 'BR_Git11', 'TC2_clean'),
  c('BR_Gil10', 'TBA_Gil10', 'Bloom', 'TBA_Git11'),
  c('BR_Gil10', 'Bloom', 'BR_Git11', 'TBA_Git11'),
  c('BR_Gil10', 'Bloom', 'Gons08', 'TC2_clean'),
  c('BR_Gil10', 'TBA_Gil10', 'Bloom', 'BR_Git11'),
  c('TBA_Gil10', 'Bloom', 'TBA_Git11', 'TC2_clean'),
  c('BR_Gil10', 'Gons08', 'TC2_clean', 'TC2_turbid')
)

LUT_OPT <- list(
  Opt_algorithm = Opt_algorithm,
  Remove_algorithm = Remove_algorithm
)