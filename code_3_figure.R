# ============================================================================
# REDE OFICIAL GINsim: DDR / MALAT1 / miR-204-5p / SIRT1 / GSDMD / PIROPTOSE
# Output de parada do ciclo: CELL_CYCLE_ARREST
# Estados estaveis exatos com e sem perturbacoes, sem Monte Carlo
# ============================================================================

setwd("~/Downloads/ARTIGO ENTRE ALUNOS/")

required_packages <- c("BoolNet", "ggplot2", "xml2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Pacotes ausentes: ", paste(missing_packages, collapse = ", "),
    ". Instale-os com: install.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(BoolNet)
  library(ggplot2)
  library(xml2)
})

# ------------------------------- CAMINHOS -----------------------------------

command_line_arguments <- commandArgs(trailingOnly = TRUE)
MODEL_FILE <- if (length(command_line_arguments) >= 1L) {
  command_line_arguments[[1L]]
} else {
  "model_pyroptosis_mir206.zginml"
}

if (!file.exists(MODEL_FILE)) {
  candidates <- list.files(
    path = ".",
    pattern = "model_pyroptosis_mir206.*[.]zginml$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(candidates) == 1L) {
    MODEL_FILE <- candidates[[1L]]
  } else {
    stop(
      "Arquivo atualizado model_pyroptosis_mir206(2).zginml nao encontrado.  ",
      "Coloque o modelo na mesma pasta do script ou informe seu caminho ",
      "como primeiro argumento do Rscript. Se houver varias versoes na pasta, ",
      "informe explicitamente o arquivo que contem CELL_CYCLE_ARREST."
    )
  }
}

MODEL_FILE <- normalizePath(MODEL_FILE, mustWork = TRUE)
OUTPUT_DIR <- file.path(getwd(), "GINsim_DDR_ON_exact_results")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------- EXTRAIR O MODELO GINSIM --------------------------

extract_regulatory_graph <- function(model_file) {
  archive_contents <- tryCatch(
    utils::unzip(model_file, list = TRUE),
    error = function(e) NULL
  )

  if (!is.null(archive_contents)) {
    graph_entries <- archive_contents$Name[
      grepl("(^|/)regulatoryGraph[.]ginml$", archive_contents$Name)
    ]

    if (length(graph_entries) != 1L) {
      stop(
        "O arquivo .zginml deveria conter exatamente um ",
        "regulatoryGraph.ginml; encontrados: ", length(graph_entries), "."
      )
    }

    extraction_directory <- tempfile("ginsim_model_")
    dir.create(extraction_directory)
    on.exit(unlink(extraction_directory, recursive = TRUE), add = TRUE)

    utils::unzip(
      zipfile = model_file,
      files = graph_entries,
      exdir = extraction_directory,
      overwrite = TRUE
    )

    extracted_file <- file.path(extraction_directory, graph_entries)
    return(xml2::read_xml(
      extracted_file,
      options = c("RECOVER", "NOERROR", "NOWARNING", "NONET")
    ))
  }

  # Compatibilidade com arquivos GINML fornecidos diretamente como XML.
  xml2::read_xml(
    model_file,
    options = c("RECOVER", "NOERROR", "NOWARNING", "NONET")
  )
}

model_document <- extract_regulatory_graph(MODEL_FILE)
graph_node <- xml2::xml_find_first(
  model_document,
  ".//graph[@class='regulatory']"
)

if (inherits(graph_node, "xml_missing")) {
  stop("Nenhum grafo regulatorio foi encontrado no arquivo .zginml.")
}

graph_id <- xml2::xml_attr(graph_node, "id")
node_order_text <- xml2::xml_attr(graph_node, "nodeorder")
node_order <- strsplit(trimws(node_order_text), "[[:space:]]+")[[1L]]
node_elements <- xml2::xml_find_all(graph_node, "./node")
node_ids <- xml2::xml_attr(node_elements, "id")

if (!setequal(node_order, node_ids)) {
  stop(
    "A lista nodeorder do GINsim nao corresponde aos nos do grafo. ",
    "Ausentes: ", paste(setdiff(node_order, node_ids), collapse = ", "),
    "; extras: ", paste(setdiff(node_ids, node_order), collapse = ", ")
  )
}

max_values <- suppressWarnings(as.integer(xml2::xml_attr(node_elements, "maxvalue")))
if (any(is.na(max_values)) || any(max_values != 1L)) {
  stop(
    "Esta versao do script aceita somente redes booleanas com maxvalue = 1."
  )
}

input_nodes <- node_ids[
  !is.na(xml2::xml_attr(node_elements, "input")) &
    xml2::xml_attr(node_elements, "input") == "true"
]

extract_rule <- function(node_element) {
  node_id <- xml2::xml_attr(node_element, "id")
  expression_node <- xml2::xml_find_first(
    node_element,
    "./value[@val='1']/exp"
  )

  if (inherits(expression_node, "xml_missing")) {
    if (node_id %in% input_nodes) {
      return(node_id)
    }
    stop("Regra booleana ausente para o no: ", node_id)
  }

  rule <- xml2::xml_attr(expression_node, "str")
  if (is.na(rule) || !nzchar(trimws(rule))) {
    stop("Expressao booleana vazia para o no: ", node_id)
  }
  trimws(rule)
}

rules_by_node <- vapply(node_elements, extract_rule, character(1))
names(rules_by_node) <- node_ids
rules_by_node <- rules_by_node[node_order]

required_nodes <- c(
  "DDR", "lncRNA_MALAT1", "miR_204_5p", "SIRT1",
  "PROLIFERATION", "RESISTANCE", "CELL_CYCLE_ARREST",
  "PYROPTOSIS", "APOPTOSIS"
)

if (!all(required_nodes %in% node_order)) {
  stop(
    "O modelo nao possui todos os nos esperados. Ausentes: ",
    paste(setdiff(required_nodes, node_order), collapse = ", ")
  )
}

if (grepl("mir204", basename(MODEL_FILE), ignore.case = TRUE) &&
    "miR_204_5p" %in% node_order) {
  message(
    "ATENCAO: o nome do arquivo menciona mir206, mas o no existente no ",
    "modelo e miR_204_5p. A analise usara o conteudo interno oficial."
  )
}

rule_table <- data.frame(
  Node = node_order,
  Boolean_rule = unname(rules_by_node),
  Is_input = node_order %in% input_nodes,
  stringsAsFactors = FALSE
)

write.csv(
  rule_table,
  file.path(OUTPUT_DIR, "01_rules_extracted_from_GINsim.csv"),
  row.names = FALSE
)

boolnet_network_file <- file.path(OUTPUT_DIR, "02_network_imported_by_BoolNet.csv")
writeLines(
  c(
    "targets, factors",
    paste(node_order, unname(rules_by_node), sep = ", ")
  ),
  boolnet_network_file
)

network <- BoolNet::loadNetwork(boolnet_network_file)

if (!setequal(network$genes, node_order)) {
  stop(
    "Os nos importados pelo BoolNet diferem dos nos do GINsim. Ausentes: ",
    paste(setdiff(node_order, network$genes), collapse = ", "),
    "; extras: ", paste(setdiff(network$genes, node_order), collapse = ", ")
  )
}

saveRDS(network, file.path(OUTPUT_DIR, "03_imported_BoolNet_network.rds"))

message(
  "Modelo oficial importado: ", graph_id,
  " | ", length(node_order), " nos | input(s): ",
  paste(input_nodes, collapse = ", ")
)

# ------------------------------ CONDICOES -----------------------------------

conditions <- list(
  "DDR ON | No perturbation" = c(DDR = 1L),
  "DDR ON | MALAT1 OE" = c(DDR = 1L, lncRNA_MALAT1 = 1L),
  "DDR ON | miR-204-5p OE" = c(DDR = 1L, miR_204_5p = 1L),
  "DDR ON | miR-204-5p KO" = c(DDR = 1L, miR_204_5p = 0L),
  "DDR ON | SIRT1 OE" = c(DDR = 1L, SIRT1 = 1L),
  "DDR ON | SIRT1 KO" = c(DDR = 1L, SIRT1 = 0L)
)

perturbation_table <- data.frame(
  Condition = names(conditions)[-1L],
  Target = c(
    "lncRNA_MALAT1", "miR_204_5p", "miR_204_5p", "SIRT1", "SIRT1"
  ),
  Type = c("OE", "OE", "KO", "OE", "KO"),
  Fixed_value = c(1L, 1L, 0L, 1L, 0L),
  Symbol = c("\u2191", "\u2191", "\u2193", "\u2191", "\u2193"),
  stringsAsFactors = FALSE
)

validate_condition <- function(condition_name, fixed_values) {
  if (is.null(names(fixed_values)) || any(names(fixed_values) == "")) {
    stop("Valores fixos sem nomes na condicao: ", condition_name)
  }
  if (!all(names(fixed_values) %in% node_order)) {
    stop(
      "No desconhecido na condicao ", condition_name, ": ",
      paste(setdiff(names(fixed_values), node_order), collapse = ", ")
    )
  }
  if (any(!fixed_values %in% c(0L, 1L))) {
    stop("Somente os valores 0 e 1 sao permitidos: ", condition_name)
  }
  if (!"DDR" %in% names(fixed_values) || fixed_values[["DDR"]] != 1L) {
    stop("DDR precisa estar fixado em 1 na condicao: ", condition_name)
  }
  invisible(TRUE)
}

invisible(Map(validate_condition, names(conditions), conditions))

output_nodes <- c(
  "PROLIFERATION", "RESISTANCE", "CELL_CYCLE_ARREST",
  "PYROPTOSIS", "APOPTOSIS"
)

output_labels <- c(
  PROLIFERATION = "Proliferation",
  RESISTANCE = "Resistance",
  CELL_CYCLE_ARREST = "Cell Cycle Arrest",
  PYROPTOSIS = "Pyroptosis",
  APOPTOSIS = "Apoptosis"
)

classify_phenotype <- function(state) {
  active_outputs <- output_nodes[state[output_nodes] == 1L]
  if (length(active_outputs) == 0L) {
    return("No active phenotype")
  }
  paste(unname(output_labels[active_outputs]), collapse = " + ")
}

# ----------------------- PONTOS FIXOS EXATOS -------------------------------

decode_fixed_points <- function(attractor_information, condition_name) {
  number_of_attractors <- length(attractor_information$attractors)
  if (number_of_attractors == 0L) {
    return(data.frame())
  }

  decoded_rows <- vector("list", number_of_attractors)

  for (attractor_index in seq_len(number_of_attractors)) {
    attractor_sequence <- as.data.frame(
      BoolNet::getAttractorSequence(attractor_information, attractor_index),
      check.names = FALSE
    )

    if (nrow(attractor_sequence) != 1L) {
      stop(
        "Foi encontrado um atrator com mais de um estado, apesar de a busca ",
        "estar restrita a pontos fixos."
      )
    }

    if (!all(network$genes %in% colnames(attractor_sequence))) {
      stop("A sequencia do atrator nao contem todos os genes da rede.")
    }

    state <- as.integer(unlist(
      attractor_sequence[1L, network$genes, drop = FALSE],
      use.names = FALSE
    ))
    names(state) <- network$genes

    decoded_rows[[attractor_index]] <- data.frame(
      Condition = condition_name,
      State_key = paste0(state[node_order], collapse = ""),
      Phenotype = classify_phenotype(state),
      as.list(state[node_order]),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  decoded <- unique(do.call(rbind, decoded_rows))

  phenotype_order <- c(
    "Proliferation + Resistance", "Proliferation", "Resistance",
    "Cell Cycle Arrest", "Pyroptosis", "Apoptosis", "No active phenotype"
  )
  phenotype_rank <- match(decoded$Phenotype, phenotype_order)
  phenotype_rank[is.na(phenotype_rank)] <- length(phenotype_order) + 1L
  decoded <- decoded[order(phenotype_rank, decoded$State_key), , drop = FALSE]
  decoded$Stable_state <- paste0("SS", seq_len(nrow(decoded)))
  decoded <- decoded[, c(
    "Condition", "Stable_state", "State_key", "Phenotype", node_order
  )]
  row.names(decoded) <- NULL
  decoded
}

analyze_condition <- function(condition_name, fixed_values) {
  validate_condition(condition_name, fixed_values)

  perturbed_network <- BoolNet::fixGenes(
    network,
    fixIndices = names(fixed_values),
    values = as.numeric(fixed_values)
  )

  attractor_information <- BoolNet::getAttractors(
    perturbed_network,
    type = "synchronous",
    method = "sat.restricted",
    maxAttractorLength = 1,
    canonical = TRUE,
    returnTable = TRUE
  )

  decode_fixed_points(attractor_information, condition_name)
}

message("Iniciando enumeracao exata de pontos fixos; Monte Carlo nao sera usado...")

state_tables <- vector("list", length(conditions))
names(state_tables) <- names(conditions)

for (condition_index in seq_along(conditions)) {
  condition_name <- names(conditions)[[condition_index]]
  message("[", condition_index, "/", length(conditions), "] ", condition_name)
  state_tables[[condition_name]] <- analyze_condition(
    condition_name,
    conditions[[condition_name]]
  )
}

if (all(vapply(state_tables, nrow, integer(1)) == 0L)) {
  stop("Nenhum ponto fixo foi encontrado.")
}

all_states <- do.call(
  rbind,
  state_tables[vapply(state_tables, nrow, integer(1)) > 0L]
)
row.names(all_states) <- NULL

all_states$Condition <- factor(
  all_states$Condition,
  levels = names(conditions)
)
all_states <- all_states[order(all_states$Condition), , drop = FALSE]
all_states$Condition <- as.character(all_states$Condition)
all_states$Row_label <- paste(
  all_states$Condition,
  all_states$Stable_state,
  all_states$Phenotype,
  sep = " | "
)

write.csv(
  all_states[, c(
    "Condition", "Stable_state", "Phenotype", "State_key", node_order
  )],
  file.path(OUTPUT_DIR, "04_exact_fixed_points_all_conditions.csv"),
  row.names = FALSE
)

condition_summary <- do.call(
  rbind,
  lapply(names(conditions), function(condition_name) {
    condition_states <- state_tables[[condition_name]]
    data.frame(
      Condition = condition_name,
      Fixed_nodes = paste(
        paste0(names(conditions[[condition_name]]), "=", conditions[[condition_name]]),
        collapse = "; "
      ),
      Number_of_fixed_points = nrow(condition_states),
      Phenotypes = if (nrow(condition_states) == 0L) {
        "None"
      } else {
        paste(unique(condition_states$Phenotype), collapse = "; ")
      },
      Monte_Carlo = "Not used",
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  condition_summary,
  file.path(OUTPUT_DIR, "05_condition_summary.csv"),
  row.names = FALSE
)

# Controle de reproducibilidade para esta versao oficial da rede.
expected_fixed_point_counts <- c(4L, 1L, 3L, 1L, 1L, 3L)
observed_fixed_point_counts <- condition_summary$Number_of_fixed_points

if (!identical(observed_fixed_point_counts, expected_fixed_point_counts)) {
  warning(
    "A quantidade de pontos fixos diferiu da auditoria esperada. Esperado: ",
    paste(expected_fixed_point_counts, collapse = ", "),
    "; observado: ", paste(observed_fixed_point_counts, collapse = ", "),
    ". Confira a versao do arquivo .zginml e do BoolNet."
  )
}

# A nova regra de CELL_CYCLE_ARREST torna os tres destinos antitumorais
# mutuamente exclusivos nos pontos fixos desta analise.
expected_phenotypes <- list(
  "DDR ON | No perturbation" = c(
    "Resistance", "Cell Cycle Arrest", "Pyroptosis", "Apoptosis"
  ),
  "DDR ON | MALAT1 OE" = "Resistance",
  "DDR ON | miR-204-5p OE" = c(
    "Cell Cycle Arrest", "Pyroptosis", "Apoptosis"
  ),
  "DDR ON | miR-204-5p KO" = "Resistance",
  "DDR ON | SIRT1 OE" = "Resistance",
  "DDR ON | SIRT1 KO" = c(
    "Cell Cycle Arrest", "Pyroptosis", "Apoptosis"
  )
)

for (condition_name in names(expected_phenotypes)) {
  observed <- sort(unique(state_tables[[condition_name]]$Phenotype))
  expected <- sort(expected_phenotypes[[condition_name]])
  if (!identical(observed, expected)) {
    warning(
      "Fenotipos inesperados em ", condition_name,
      ". Esperado: ", paste(expected, collapse = "; "),
      "; observado: ", paste(observed, collapse = "; "), "."
    )
  }
}

# Alcancabilidade = presenca nos pontos fixos, e nao probabilidade.
phenotype_reachability <- do.call(
  rbind,
  lapply(names(conditions), function(condition_name) {
    condition_states <- state_tables[[condition_name]]
    total_states <- nrow(condition_states)

    do.call(
      rbind,
      lapply(output_nodes, function(output_node) {
        active_states <- if (total_states == 0L) {
          0L
        } else {
          sum(condition_states[[output_node]] == 1L)
        }

        status <- if (active_states == 0L) {
          "Absent"
        } else if (active_states == total_states) {
          "Present in all fixed points"
        } else {
          "Present in some fixed points"
        }

        data.frame(
          Condition = condition_name,
          Output = unname(output_labels[[output_node]]),
          Active_fixed_points = active_states,
          Total_fixed_points = total_states,
          Status = status,
          stringsAsFactors = FALSE
        )
      })
    )
  })
)

write.csv(
  phenotype_reachability,
  file.path(OUTPUT_DIR, "06_phenotype_reachability.csv"),
  row.names = FALSE
)

# ------------------------------- HEATMAPS -----------------------------------

make_long_state_table <- function(states) {
  do.call(
    rbind,
    lapply(seq_len(nrow(states)), function(row_index) {
      values <- as.integer(unlist(
        states[row_index, node_order, drop = FALSE],
        use.names = FALSE
      ))

      fill_group <- ifelse(values == 1L, "Active node", "Inactive node")
      active_outputs <- node_order %in% output_nodes & values == 1L
      fill_group[active_outputs] <- unname(output_labels[node_order[active_outputs]])

      data.frame(
        Condition = states$Condition[[row_index]],
        Row_label = states$Row_label[[row_index]],
        Node = node_order,
        Node_index = seq_along(node_order),
        Value = values,
        Fill_group = fill_group,
        stringsAsFactors = FALSE
      )
    })
  )
}

build_stable_state_heatmap <- function(states, title, filename_prefix) {
  if (nrow(states) == 0L) {
    warning("Heatmap ignorado por ausencia de estados: ", title)
    return(invisible(NULL))
  }

  state_long <- make_long_state_table(states)
  displayed_rows <- unique(states$Row_label)
  state_long$Row_label <- factor(
    state_long$Row_label,
    levels = rev(displayed_rows)
  )

  perturbation_overlay <- merge(
    unique(states[, c("Condition", "Row_label")]),
    perturbation_table,
    by = "Condition",
    all = FALSE,
    sort = FALSE
  )
  perturbation_overlay$Node_index <- match(
    perturbation_overlay$Target,
    node_order
  )
  perturbation_overlay$Row_label <- factor(
    perturbation_overlay$Row_label,
    levels = rev(displayed_rows)
  )

  input_position <- match("DDR", node_order)
  output_positions <- match(output_nodes, node_order)
  number_of_rows <- length(displayed_rows)

  heatmap_plot <- ggplot(
    state_long,
    aes(x = Node_index, y = Row_label, fill = Fill_group)
  ) +
    geom_tile(color = "white", linewidth = 0.45, width = 0.96, height = 0.90) +
    annotate(
      "rect",
      xmin = input_position - 0.52,
      xmax = input_position + 0.52,
      ymin = 0.45,
      ymax = number_of_rows + 0.55,
      fill = NA,
      color = "#E9A23B",
      linewidth = 1.1
    ) +
    annotate(
      "rect",
      xmin = min(output_positions) - 0.52,
      xmax = max(output_positions) + 0.52,
      ymin = 0.45,
      ymax = number_of_rows + 0.55,
      fill = NA,
      color = "#2A9D8F",
      linewidth = 1.1
    ) +
    scale_fill_manual(
      values = c(
        "Inactive node" = "#F1F4F7",
        "Active node" = "#2F5D8C",
        "Proliferation" = "#A65A3A",
        "Resistance" = "#5B4B8A",
        "Cell Cycle Arrest" = "#D4A017",
        "Pyroptosis" = "#168C80",
        "Apoptosis" = "#D65A4A"
      ),
      name = "Boolean state / active output",
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = seq_along(node_order),
      labels = node_order,
      expand = expansion(add = 0.55)
    ) +
    scale_y_discrete(drop = FALSE) +
    labs(
      title = title,
      subtitle = paste0(
        "Official GINsim rules; DDR fixed at ON; exact SAT fixed points; ",
        "no Monte Carlo sampling"
      ),
      x = NULL,
      y = NULL,
      caption = paste0(
        "Light cells = 0; blue cells = 1; colored output cells identify ",
        "active phenotypes. Arrows mark fixed KO/OE nodes."
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 58, hjust = 1, vjust = 1, size = 7.7, color = "#1F2933"
      ),
      axis.text.y = element_text(size = 7.5, color = "#1F2933"),
      plot.title = element_text(face = "bold", size = 16, color = "#17202A"),
      plot.subtitle = element_text(color = "#52616B"),
      plot.caption = element_text(hjust = 0, color = "#52616B"),
      legend.position = "bottom",
      plot.margin = margin(12, 18, 12, 12)
    )

  if (nrow(perturbation_overlay) > 0L) {
    heatmap_plot <- heatmap_plot +
      geom_tile(
        data = perturbation_overlay,
        aes(x = Node_index, y = Row_label),
        inherit.aes = FALSE,
        width = 0.82,
        height = 0.76,
        fill = NA,
        color = "#C1121F",
        linewidth = 1.0
      ) +
      geom_text(
        data = perturbation_overlay,
        aes(x = Node_index, y = Row_label, label = Symbol),
        inherit.aes = FALSE,
        color = "#F28E2B",
        size = 5.0,
        fontface = "bold"
      )
  }

  figure_height <- max(5.8, 0.36 * nrow(states) + 3.8)
  png_file <- file.path(OUTPUT_DIR, paste0(filename_prefix, ".png"))
  pdf_file <- file.path(OUTPUT_DIR, paste0(filename_prefix, ".pdf"))

  ggsave(
    png_file,
    heatmap_plot,
    width = 18,
    height = figure_height,
    units = "in",
    dpi = 400,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    pdf_file,
    heatmap_plot,
    width = 18,
    height = figure_height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )

  invisible(heatmap_plot)
}

unperturbed_states <- all_states[
  all_states$Condition == "DDR ON | No perturbation",
  ,
  drop = FALSE
]
perturbed_states <- all_states[
  all_states$Condition != "DDR ON | No perturbation",
  ,
  drop = FALSE
]

build_stable_state_heatmap(
  unperturbed_states,
  "DDR ON stable states without perturbation",
  "07_heatmap_without_perturbation"
)

build_stable_state_heatmap(
  perturbed_states,
  "DDR ON stable states under the requested perturbations",
  "08_heatmap_with_perturbations"
)

build_stable_state_heatmap(
  all_states,
  "DDR ON stable-state comparison: unperturbed and perturbed conditions",
  "09_heatmap_combined_comparison"
)

# Heatmap resumido de alcance dos outputs. Os numeros sao contagens de pontos
# fixos, nao percentuais nem probabilidades de convergencia.
phenotype_reachability$Condition <- factor(
  phenotype_reachability$Condition,
  levels = rev(names(conditions))
)
phenotype_reachability$Output <- factor(
  phenotype_reachability$Output,
  levels = unname(output_labels[output_nodes])
)
phenotype_reachability$Status <- factor(
  phenotype_reachability$Status,
  levels = c(
    "Absent", "Present in some fixed points", "Present in all fixed points"
  )
)

reachability_plot <- ggplot(
  phenotype_reachability,
  aes(x = Output, y = Condition, fill = Status)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = paste0(Active_fixed_points, "/", Total_fixed_points)),
    fontface = "bold",
    size = 3.5,
    color = "#17202A"
  ) +
  scale_fill_manual(
    values = c(
      "Absent" = "#EDF1F4",
      "Present in some fixed points" = "#F2C14E",
      "Present in all fixed points" = "#4C956C"
    ),
    name = "Output reachability",
    drop = FALSE
  ) +
  labs(
    title = "Phenotype reachability under DDR ON",
    subtitle = "Cell labels: fixed points with active output / total fixed points",
    x = NULL,
    y = NULL,
    caption = "Counts of exact fixed points are not probabilities or cell frequencies."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#52616B"),
    plot.caption = element_text(hjust = 0, color = "#52616B"),
    legend.position = "bottom"
  )

ggsave(
  file.path(OUTPUT_DIR, "10_phenotype_reachability_heatmap.png"),
  reachability_plot,
  width = 12.5,
  height = 7.8,
  units = "in",
  dpi = 400,
  bg = "white"
)
ggsave(
  file.path(OUTPUT_DIR, "10_phenotype_reachability_heatmap.pdf"),
  reachability_plot,
  width = 12.5,
  height = 7.8,
  units = "in",
  bg = "white"
)


# =============================================================================
# FIGURE 3 - MODEL-PREDICTED STABLE STATES AND LITERATURE CORRESPONDENCE IN HCC
# =============================================================================
# This panel is built directly from the exact fixed points already stored in
# `all_states`. No Boolean state is manually reconstructed for the figure.
#
# Interpretation:
#   - E1 = Boolean node fixed at 1 (gain-of-function constraint)
#   - KO = Boolean node fixed at 0 (loss-of-function constraint)
#   - Fixed points indicate qualitative stable-state reachability, not
#     frequencies, probabilities, or basin sizes.
#   - "Literature-supported" indicates qualitative correspondence with the
#     cited literature and should not be read as direct experimental validation
#     of the complete Boolean state.
# =============================================================================

figure_condition_order <- c(
  "DDR ON | MALAT1 OE",
  "DDR ON | miR-204-5p OE",
  "DDR ON | miR-204-5p KO",
  "DDR ON | SIRT1 OE",
  "DDR ON | SIRT1 KO"
)

condition_display <- c(
  "DDR ON | MALAT1 OE"      = "MALAT1 E1",
  "DDR ON | miR-204-5p OE" = "miR-204-5p E1",
  "DDR ON | miR-204-5p KO" = "miR-204-5p KO",
  "DDR ON | SIRT1 OE"       = "SIRT1 E1",
  "DDR ON | SIRT1 KO"       = "SIRT1 KO"
)

# -----------------------------------------------------------------------------
# Literature correspondence shown in the figure.
# Keep this table synchronized with the final manuscript references.
# -----------------------------------------------------------------------------
literature_map <- data.frame(
  Condition = c(
    "DDR ON | MALAT1 OE",
    "DDR ON | miR-204-5p OE",
    "DDR ON | miR-204-5p OE",
    "DDR ON | miR-204-5p OE",
    "DDR ON | miR-204-5p KO",
    "DDR ON | SIRT1 OE",
    "DDR ON | SIRT1 KO",
    "DDR ON | SIRT1 KO",
    "DDR ON | SIRT1 KO"
  ),
  Phenotype = c(
    "Resistance",
    "Cell Cycle Arrest",
    "Pyroptosis",
    "Apoptosis",
    "Resistance",
    "Resistance",
    "Cell Cycle Arrest",
    "Pyroptosis",
    "Apoptosis"
  ),
  Reference = c(
    "Hou et al., 2017",
    "Jiang et al., 2016",
    NA,
    "Jiang et al., 2016",
    "Hou et al., 2017",
    "Wang et al., 2021",
    "Jiang et al., 2016",
    "Chen et al., 2020",
    "Jiang et al., 2016"
  ),
  Evidence = c(
    "Literature-supported",
    "Literature-supported",
    "Model-derived prediction",
    "Literature-supported",
    "Literature-supported",
    "Literature-supported",
    "Literature-supported",
    "Literature-supported",
    "Literature-supported"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Extract only the requested perturbed exact fixed points.
# -----------------------------------------------------------------------------
figure_states <- all_states[
  all_states$Condition %in% figure_condition_order,
  ,
  drop = FALSE
]

if (nrow(figure_states) == 0L) {
  stop("Figure 3: no perturbed fixed points were found in all_states.")
}

figure_states$Condition <- factor(
  figure_states$Condition,
  levels = figure_condition_order
)

phenotype_rank <- c(
  "Proliferation"       = 1L,
  "Resistance"          = 2L,
  "Cell Cycle Arrest"   = 3L,
  "Pyroptosis"          = 4L,
  "Apoptosis"           = 5L,
  "No active phenotype" = 6L
)

figure_states$Phenotype_rank <- unname(
  phenotype_rank[figure_states$Phenotype]
)
figure_states$Phenotype_rank[is.na(figure_states$Phenotype_rank)] <- 99L

figure_states <- figure_states[
  order(
    figure_states$Condition,
    figure_states$Phenotype_rank,
    figure_states$State_key
  ),
  ,
  drop = FALSE
]
row.names(figure_states) <- NULL

# Highest Y value is rendered at the top.
figure_states$Plot_row <- seq_len(nrow(figure_states))
figure_states$Y <- nrow(figure_states) - figure_states$Plot_row + 1L
figure_states$Perturbation <- unname(
  condition_display[as.character(figure_states$Condition)]
)

# -----------------------------------------------------------------------------
# Join literature annotations line-by-line by condition and phenotype.
# -----------------------------------------------------------------------------
literature_key <- paste(
  literature_map$Condition,
  literature_map$Phenotype,
  sep = " || "
)
state_key <- paste(
  as.character(figure_states$Condition),
  figure_states$Phenotype,
  sep = " || "
)
match_index <- match(state_key, literature_key)

figure_states$Reference <- literature_map$Reference[match_index]
figure_states$Evidence <- literature_map$Evidence[match_index]
figure_states$Reference_display <- ifelse(
  is.na(figure_states$Reference),
  "\u2014",
  figure_states$Reference
)

if (any(is.na(match_index))) {
  warning(
    paste0(
      "Figure 3: literature annotation missing for: ",
      paste(state_key[is.na(match_index)], collapse = "; ")
    )
  )
}

# -----------------------------------------------------------------------------
# Convert exact Boolean states to a long table for ggplot2.
# -----------------------------------------------------------------------------
figure_long <- do.call(
  rbind,
  lapply(seq_len(nrow(figure_states)), function(i) {
    values <- as.integer(unlist(
      figure_states[i, node_order, drop = FALSE],
      use.names = FALSE
    ))

    fill_group <- ifelse(values == 1L, "Active node", "Inactive node")

    active_outputs <- node_order %in% output_nodes & values == 1L
    fill_group[active_outputs] <- unname(
      output_labels[node_order[active_outputs]]
    )

    data.frame(
      Y = figure_states$Y[[i]],
      X = seq_along(node_order),
      Node = node_order,
      Value = values,
      Fill_group = fill_group,
      stringsAsFactors = FALSE
    )
  })
)

fill_levels <- c(
  "Inactive node",
  "Active node",
  "Proliferation",
  "Resistance",
  "Cell Cycle Arrest",
  "Pyroptosis",
  "Apoptosis"
)

figure_long$Fill_group <- factor(
  figure_long$Fill_group,
  levels = fill_levels
)

# -----------------------------------------------------------------------------
# Perturbation arrows are inherited from perturbation_table.
# -----------------------------------------------------------------------------
overlay_rows <- unique(
  figure_states[, c("Condition", "Y"), drop = FALSE]
)
overlay_rows$Condition <- as.character(overlay_rows$Condition)

overlay_match <- match(
  overlay_rows$Condition,
  perturbation_table$Condition
)
overlay_rows$Target <- perturbation_table$Target[overlay_match]
overlay_rows$Symbol <- perturbation_table$Symbol[overlay_match]
overlay_rows$X <- match(overlay_rows$Target, node_order)
overlay_rows <- overlay_rows[!is.na(overlay_rows$X), , drop = FALSE]

# -----------------------------------------------------------------------------
# Group geometry used for perturbation labels and separators.
# -----------------------------------------------------------------------------
group_info <- do.call(
  rbind,
  lapply(figure_condition_order, function(cond) {
    y_values <- figure_states$Y[
      as.character(figure_states$Condition) == cond
    ]

    if (length(y_values) == 0L) {
      return(NULL)
    }

    data.frame(
      Condition = cond,
      Perturbation = unname(condition_display[[cond]]),
      Y_center = mean(range(y_values)),
      Y_min = min(y_values),
      Y_max = max(y_values),
      stringsAsFactors = FALSE
    )
  })
)

separator_y <- if (nrow(group_info) > 1L) {
  group_info$Y_min[-nrow(group_info)] - 0.5
} else {
  numeric(0)
}

# -----------------------------------------------------------------------------
# Publication labels for nodes.
# -----------------------------------------------------------------------------
node_display <- node_order
node_display[node_order == "miR_204_5p"] <- "miR-204-5p"
node_display[node_order == "CELL_CYCLE_ARREST"] <- "CELL-CYCLE ARREST"

input_position <- match("DDR", node_order)
output_positions <- match(output_nodes, node_order)

if (is.na(input_position)) {
  stop("DDR was not found in node_order.")
}
if (any(is.na(output_positions))) {
  stop("One or more phenotype output nodes were not found in node_order.")
}

n_nodes <- length(node_order)
n_rows <- nrow(figure_states)

x_perturbation <- n_nodes + 1.10
x_reference <- n_nodes + 5.30
x_evidence <- n_nodes + 8.75
x_right <- n_nodes + 15.25

supported_rows <- figure_states[
  !is.na(figure_states$Evidence) &
    figure_states$Evidence == "Literature-supported",
  ,
  drop = FALSE
]

prediction_rows <- figure_states[
  !is.na(figure_states$Evidence) &
    figure_states$Evidence == "Model-derived prediction",
  ,
  drop = FALSE
]

# -----------------------------------------------------------------------------
# Final figure.
# -----------------------------------------------------------------------------
p_validation <- ggplot2::ggplot() +
  ggplot2::geom_tile(
    data = figure_long,
    ggplot2::aes(x = X, y = Y, fill = Fill_group),
    width = 0.96,
    height = 0.90,
    color = "white",
    linewidth = 0.42
  ) +
  ggplot2::annotate(
    "rect",
    xmin = input_position - 0.52,
    xmax = input_position + 0.52,
    ymin = 0.45,
    ymax = n_rows + 0.55,
    fill = NA,
    color = "#E69F00",
    linewidth = 1.05
  ) +
  ggplot2::annotate(
    "rect",
    xmin = min(output_positions) - 0.52,
    xmax = max(output_positions) + 0.52,
    ymin = 0.45,
    ymax = n_rows + 0.55,
    fill = NA,
    color = "#009E73",
    linewidth = 1.05
  ) +
  ggplot2::geom_tile(
    data = overlay_rows,
    ggplot2::aes(x = X, y = Y),
    inherit.aes = FALSE,
    width = 0.80,
    height = 0.74,
    fill = NA,
    color = "#D55E00",
    linewidth = 0.95
  ) +
  ggplot2::geom_text(
    data = overlay_rows,
    ggplot2::aes(x = X, y = Y, label = Symbol),
    inherit.aes = FALSE,
    color = "#E69F00",
    size = 5.0,
    fontface = "bold"
  ) +
  ggplot2::geom_hline(
    yintercept = separator_y,
    color = "#D7DEE7",
    linewidth = 0.50
  ) +
  ggplot2::geom_text(
    data = group_info,
    ggplot2::aes(x = x_perturbation, y = Y_center, label = Perturbation),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 4.0,
    color = "#17202A"
  ) +
  ggplot2::geom_text(
    data = figure_states,
    ggplot2::aes(x = x_reference, y = Y, label = Reference_display),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.45,
    fontface = "italic",
    color = "#17202A"
  ) +
  ggplot2::geom_label(
    data = supported_rows,
    ggplot2::aes(x = x_evidence, y = Y, label = "LITERATURE-SUPPORTED"),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.9,
    label.size = 0.55,
    label.padding = grid::unit(0.18, "lines"),
    label.r = grid::unit(0.18, "lines"),
    fill = "#E9F8F2",
    color = "#007A5A"
  ) +
  ggplot2::geom_label(
    data = prediction_rows,
    ggplot2::aes(
      x = x_evidence,
      y = Y,
      label = paste0(
        "MODEL-DERIVED PREDICTION\n",
        "Direct experimental validation required"
      )
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.65,
    lineheight = 0.92,
    label.size = 0.65,
    label.padding = grid::unit(0.18, "lines"),
    label.r = grid::unit(0.18, "lines"),
    fill = "#FFF4E5",
    color = "#C65D00"
  ) +
  ggplot2::annotate(
    "segment",
    x = -0.85,
    xend = -0.85,
    y = 0.55,
    yend = n_rows + 0.45,
    linewidth = 0.8,
    color = "#17202A"
  ) +
  ggplot2::annotate(
    "segment",
    x = -0.85,
    xend = -0.50,
    y = 0.55,
    yend = 0.55,
    linewidth = 0.8,
    color = "#17202A"
  ) +
  ggplot2::annotate(
    "segment",
    x = -0.85,
    xend = -0.50,
    y = n_rows + 0.45,
    yend = n_rows + 0.45,
    linewidth = 0.8,
    color = "#17202A"
  ) +
  ggplot2::annotate(
    "text",
    x = -1.32,
    y = (n_rows + 1) / 2,
    label = "Stable states",
    angle = 90,
    size = 4.0,
    fontface = "bold",
    color = "#17202A"
  ) +
  ggplot2::annotate(
    "text",
    x = input_position,
    y = n_rows + 1.20,
    label = "Input\n(DDR = 1)",
    size = 3.9,
    fontface = "bold",
    lineheight = 0.95
  ) +
  ggplot2::annotate(
    "text",
    x = mean(range(output_positions)),
    y = n_rows + 1.20,
    label = "Outputs\n(phenotypes)",
    size = 3.9,
    fontface = "bold",
    lineheight = 0.95
  ) +
  ggplot2::annotate(
    "text",
    x = x_perturbation,
    y = n_rows + 1.20,
    label = "Perturbation",
    hjust = 0,
    size = 3.9,
    fontface = "bold"
  ) +
  ggplot2::annotate(
    "text",
    x = x_reference,
    y = n_rows + 1.20,
    label = "Reference",
    hjust = 0,
    size = 3.9,
    fontface = "bold"
  ) +
  ggplot2::annotate(
    "text",
    x = x_evidence,
    y = n_rows + 1.20,
    label = "Evidence",
    hjust = 0,
    size = 3.9,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Inactive node" = "#F2F6FA",
      "Active node" = "#2F5D8C",
      "Proliferation" = "#A65A3A",
      "Resistance" = "#6A3D9A",
      "Cell Cycle Arrest" = "#E69F00",
      "Pyroptosis" = "#009E73",
      "Apoptosis" = "#E84A5F"
    ),
    breaks = c(
      "Active node",
      "Inactive node",
      "Proliferation",
      "Resistance",
      "Cell Cycle Arrest",
      "Pyroptosis",
      "Apoptosis"
    ),
    labels = c(
      "Active node",
      "Inactive node",
      "Proliferation",
      "Resistance",
      "Cell-cycle arrest",
      "Pyroptosis",
      "Apoptosis"
    ),
    drop = FALSE,
    name = NULL
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq_along(node_order),
    labels = node_display,
    limits = c(-1.65, x_right),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(-0.10, n_rows + 1.75),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = paste0(
      "Validation of Model Predictions against Experimental Observations in HCC"
    ),
    subtitle = paste0(
      "DDR fixed at 1; exact SAT-restricted fixed points. ",
      "Literature correspondence indicates qualitative concordance."
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "E1 = Boolean node fixed at 1; KO = node fixed at 0. ",
      "Fixed points indicate qualitative stable-state reachability, ",
      "not frequencies or probabilities. Literature support does not ",
      "constitute direct validation of the complete Boolean state."
    )
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 52,
      hjust = 1,
      vjust = 1,
      size = 8.0,
      color = "#1F2933"
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 18,
      hjust = 0.5,
      color = "#111827",
      margin = ggplot2::margin(b = 5)
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10.5,
      hjust = 0.5,
      color = "#4B5563",
      margin = ggplot2::margin(b = 12)
    ),
    plot.caption = ggplot2::element_text(
      size = 8.8,
      hjust = 0,
      color = "#4B5563",
      margin = ggplot2::margin(t = 10)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = ggplot2::element_text(size = 9.3),
    legend.key.width = grid::unit(1.0, "lines"),
    legend.key.height = grid::unit(1.0, "lines"),
    plot.margin = ggplot2::margin(18, 28, 16, 30)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(color = NA)
    )
  )

# Export publication-quality versions.
ggplot2::ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_3_model_predictions_vs_literature.png"
  ),
  plot = p_validation,
  width = 18,
  height = 10.125,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_3_model_predictions_vs_literature.pdf"
  ),
  plot = p_validation,
  width = 18,
  height = 10.125,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

print(p_validation)

# --------------------------- RELATORIO METODOLOGICO -------------------------

methodological_report <- c(
  paste0("Source model: ", MODEL_FILE),
  paste0("GINsim graph id: ", graph_id),
  paste0("Number of nodes: ", length(node_order)),
  paste0("Input nodes: ", paste(input_nodes, collapse = ", ")),
  "DDR constraint: DDR = 1 in the unperturbed reference and all perturbations",
  paste0("Conditions: ", paste(names(conditions), collapse = "; ")),
  "Attractor method: exact SAT-restricted enumeration",
  "Maximum attractor length: 1 (fixed points only)",
  "Monte Carlo sampling: not used",
  "Convergence probabilities: not calculated",
  paste0(
    "Fixed-point counts: ",
    paste(
      paste0(condition_summary$Condition, "=", condition_summary$Number_of_fixed_points),
      collapse = "; "
    )
  ),
  paste0(
    "Important naming audit: archive=", basename(MODEL_FILE),
    "; internal microRNA node=miR_204_5p"
  ),
  paste0(
    "Cell-cycle-arrest output: CELL_CYCLE_ARREST = ",
    rules_by_node[["CELL_CYCLE_ARREST"]]
  ),
  "Interpretation: exact fixed points are stable under synchronous and asynchronous update semantics.",
  "Limitation: complex asynchronous attractors and their reachability were not assessed."
)

writeLines(
  methodological_report,
  file.path(OUTPUT_DIR, "11_methodological_report.txt")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "12_sessionInfo.txt")
)

message("Analise concluida.")
message("Resultados salvos em: ", normalizePath(OUTPUT_DIR, mustWork = TRUE))
print(condition_summary)

