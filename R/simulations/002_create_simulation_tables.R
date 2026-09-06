
################################################################################
# 2. Read power from RData
################################################################################

read_power <- function(file) {
  
  env <- new.env()
  
  objects <- load(
    file,
    envir = env
  )
  
  power_object <- grep(
    "^potencia_",
    objects,
    value = TRUE
  )
  
  if (length(power_object) != 1) {
    stop(
      "Expected one 'potencia_' object in: ",
      basename(file)
    )
  }
  
  as.numeric(
    env[[power_object]]
  )
}


################################################################################
# 3. Read one scenario
################################################################################

read_scenario <- function(scenario) {
  
  files <- list.files(
    scenario_paths[[scenario]],
    pattern = paste0(scenario, "\\.RData$"),
    full.names = TRUE
  )
  
  
  ############################################################################
  # No RData files
  ############################################################################
  
  if (length(files) == 0) {
    
    message(
      "\tScenario ", scenario,
      ": 0 RData files found - creating empty table"
    )
    
    return(
      data.frame(
        statistic = character(0),
        method    = character(0),
        n1        = numeric(0),
        n2        = numeric(0),
        power     = I(list()),
        sample    = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  
  ############################################################################
  # Available RData files
  ############################################################################
  
  message(
    "\tScenario ", scenario,
    ": using ", length(files),
    " RData files"
  )
  
  
  ############################################################################
  # Read files
  ############################################################################
  
  results <- do.call(
    rbind,
    lapply(
      files,
      function(file) {
        
        name <- sub(
          "\\.RData$",
          "",
          basename(file)
        )
        
        parts <- strsplit(
          name,
          "_"
        )[[1]]
        
        # Expected examples:
        # t2_montecarlo_100_50
        # l1std_boot_100_50_1A
        # l2std_perms_200_100_1A
        #
        # The possible scenario suffix is ignored.
        
        if (length(parts) < 4) {
          stop(
            "Invalid filename: ",
            basename(file)
          )
        }
        
        power <- read_power(file)
        
        data.frame(
          statistic = parts[1],
          method    = parts[2],
          n1        = as.numeric(parts[3]),
          n2        = as.numeric(parts[4]),
          power     = I(list(power)),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  
  ############################################################################
  # Create sample identifier
  ############################################################################
  
  results$sample <- paste(
    results$n1,
    results$n2,
    sep = "_"
  )
  
  
  ############################################################################
  # Define ordering
  ############################################################################
  
  results$sample <- factor(
    results$sample,
    levels = sample_order[[scenario]]
  )
  
  results$statistic <- factor(
    results$statistic,
    levels = stat_order
  )
  
  results$method <- factor(
    results$method,
    levels = method_order
  )
  
  
  ############################################################################
  # Validate filenames
  ############################################################################
  
  if (
    any(is.na(results$sample)) ||
    any(is.na(results$statistic)) ||
    any(is.na(results$method))
  ) {
    stop(
      "Unrecognized filename in Scenario ",
      scenario
    )
  }
  
  
  ############################################################################
  # Order results
  ############################################################################
  
  results <- results[
    order(
      results$sample,
      results$statistic,
      results$method
    ),
  ]
  
  rownames(results) <- NULL
  
  
  ############################################################################
  # Return results
  ############################################################################
  
  results
}


################################################################################
# 4. Format numbers
################################################################################

format_power <- function(x) {
  sprintf("%.3f", x)
}


################################################################################
# 5. Create one LaTeX row
################################################################################

create_row <- function(statistic, method, power, scenario_type) {
  
  statistic_text <- stat_latex[[statistic]]
  method_text    <- method_latex[[method]]
  
  power <- format_power(power)
  
  # Scenario 1: delta = 0, 0.3, 0.5, 0.7, 0.9
  if (scenario_type == "1") {
    
    values <- power
    
  }
  
  # Scenario 2:
  # T2     -> 0, 0.02, 0.03, 0.05, -, -, -
  # L1/L2  -> 0, -, -, -, 0.3, 0.5, 0.6
  if (scenario_type == "2") {
    
    if (statistic == "t2") {
      
      values <- c(
        power,
        "-",
        "-",
        "-"
      )
      
    } else {
      
      values <- c(
        power[1],
        "-",
        "-",
        "-",
        power[2:4]
      )
    }
  }
  
  paste0(
    " & ",
    statistic_text,
    " & ",
    method_text,
    " & ",
    paste(
      values,
      collapse = " & "
    ),
    " \\\\"
  )
}


################################################################################
# 6. Create one scenario block
################################################################################

create_scenario_block <- function(data, scenario, index) {
  
  scenario_type <- substr(
    scenario,
    1,
    1
  )
  
  n1 <- unique(data$n1)
  n2 <- unique(data$n2)
  
  n_columns <- ifelse(
    scenario_type == "1",
    8,
    10
  )
  
  lines <- c(
    paste0(
      "\\multicolumn{",
      n_columns,
      "}{@{}l}{\\textbf{Scenario ",
      scenario,
      ".",
      index,
      ": } $n_1=",
      n1,
      ",\\, n_2=",
      n2,
      "$} \\\\"
    ),
    "\\midrule"
  )
  
  for (stat in stat_order) {
    
    temp <- data[
      as.character(data$statistic) == stat,
    ]
    
    for (j in seq_len(nrow(temp))) {
      
      row <- create_row(
        statistic = stat,
        method = as.character(temp$method[j]),
        power = temp$power[[j]],
        scenario_type = scenario_type
      )
      
      # Statistic name only in first method row
      if (j > 1) {
        
        row <- sub(
          paste0(
            " & ",
            stat_latex[[stat]],
            " & "
          ),
          " &  & ",
          row,
          fixed = TRUE
        )
      }
      
      lines <- c(
        lines,
        row
      )
    }
  }
  
  lines
}


################################################################################
# 7. Create complete LaTeX table
################################################################################

create_latex_table <- function(scenario) {
  
  results <- read_scenario(
    scenario
  )
  
  scenario_type <- substr(
    scenario,
    1,
    1
  )
  
  samples <- sample_order[[scenario]]
  
  blocks <- list()
  
  if (nrow(results) > 0) {
    
    blocks <- lapply(
      seq_along(samples),
      function(i) {
        
        data_sample <- results[
          results$sample == samples[i],
        ]
        
        if (nrow(data_sample) > 0) {
          
          create_scenario_block(
            data_sample,
            scenario,
            i
          )
          
        } else {
          
          NULL
        }
      }
    )
    
    blocks <- Filter(
      Negate(is.null),
      blocks
    )
  }
  
  ##########################################################################
  # Header - Scenarios 1A and 1B
  ##########################################################################
  
  if (scenario_type == "1") {
    
    header <- c(
      
      "\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}",
      "l",
      ">{\\raggedright\\arraybackslash}p{3.2cm}",
      "c",
      "*{5}{S[table-format=1.3]}",
      "@{}}",
      
      "\\toprule",
      
      paste0(
        " & \\textbf{Statistic}",
        " & \\textbf{Method}",
        " & \\multicolumn{5}{c}{$\\boldsymbol{\\delta}$} \\\\"
      ),
      
      "\\cmidrule(l){4-8}",
      
      paste0(
        " & & & {0}",
        " & {0.3}",
        " & {0.5}",
        " & {0.7}",
        " & {0.9} \\\\"
      ),
      
      "\\midrule"
    )
  }
  
  
  ##########################################################################
  # Header - Scenarios 2A and 2B
  ##########################################################################
  
  if (scenario_type == "2") {
    
    header <- c(
      
      "\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}",
      "l",
      ">{\\raggedright\\arraybackslash}p{3.2cm}",
      "c",
      "*{7}{S[table-format=1.3]}",
      "@{}}",
      
      "\\toprule",
      
      paste0(
        " & \\textbf{Statistic}",
        " & \\textbf{Method}",
        " & \\multicolumn{4}{c}{$\\boldsymbol{\\eta}$ (Hotelling $T^{2}$)}",
        " & \\multicolumn{3}{c}{$\\boldsymbol{\\eta}$ ($L^{1}$/$L^{2}$)} \\\\"
      ),
      
      "\\cmidrule(lr){4-7}",
      "\\cmidrule(l){8-10}",
      
      paste0(
        " & & & {0}",
        " & {0.02}",
        " & {0.03}",
        " & {0.05}",
        " & {0.3}",
        " & {0.5}",
        " & {0.6} \\\\"
      ),
      
      "\\midrule"
    )
  }
  
  
  ##########################################################################
  # Join everything
  ##########################################################################
  
  latex <- header
  
  for (i in seq_along(blocks)) {
    
    if (i > 1) {
      latex <- c(
        latex,
        "\\midrule"
      )
    }
    
    latex <- c(
      latex,
      blocks[[i]]
    )
  }
  
  c(
    latex,
    "\\bottomrule",
    "\\end{tabular*}"
  )
}



################################################################################
# SIMULATION TABLES
# Scenarios 1A, 1B, 2A and 2B
################################################################################
source("R/001_helper_functions.R")
message("\tCreating simulation tables")

################################################################################
# 1. Configuration
################################################################################

scenario_paths <- c(
  # "1A" = "results/simulations/Scenario1A",
  "1A" = "results/simulations",
  "1B" = "results/simulations",
  "2A" = "results/simulations/",
  "2B" = "results/simulations/"
)

output_path <- paste0("results/simulations/tables")

sample_order <- list(
  "1A" = c("100_50", "150_100", "200_100"),
  "1B" = c("100_50", "150_100", "200_100"),
  "2A" = c("40_20", "60_40", "80_40"),
  "2B" = c("40_20", "60_40", "80_40")
)

stat_order <- c(
  "t2",
  "l1std",
  "l2std"
)

method_order <- c(
  "montecarlo",
  "boot",
  "perms"
)

stat_latex <- c(
  t2    = "Hotelling $T^{2}$",
  l1std = "$L^{1}$ std",
  l2std = "$L^{2}$ std"
)

method_latex <- c(
  montecarlo = "M",
  boot  = "B",
  perms      = "P"
)



################################################################################
# 2. Generate all four tables
################################################################################

scenarios <- c(
  "1A",
  "1B",
  "2A",
  "2B"
)

output_files <- c(
  "1A" = "res_1A.tex",
  "1B" = "res_1B.tex",
  "2A" = "res_2A.tex",
  "2B" = "res_2B.tex"
)

for (scenario in scenarios) {
  
  latex <- create_latex_table(
    scenario
  )
  
  con <- file(
    file.path(
      output_path,
      output_files[[scenario]]
    ),
    open = "w",
    encoding = "UTF-8"
  )
  
  writeLines(
    latex,
    con = con
  )
  
  close(con)
  
  message(
    "\tCreated: ",
    output_files[[scenario]]
  )
}
