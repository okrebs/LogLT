#' Convert standardized IO Table to simulation inputs
#'
#' \code{convert_io} takes an input output table in standardized long format and
#' converts it to the necessary model inputs for simulations.
#'
#' @details Our simulation requires data on country-sector revenues, trade
#' deficits, sectoral import shares by use category (final or sectoral
#' intermediate), country-sector consumption shares, labor cost shares and cost
#' shares of material inputs, all of which can be derived from an input output
#' table.
#'
#' @param iot An input output table in standardized long format (see the package
#' 'iotr') with the columns, \code{origin}, \code{sector}, \code{destination},
#' \code{use} and \code{flow}. Alternative column names can be specified with
#' the following parameters
#' @param orig_col alternative name for 'origin' column
#' @param sec_col alternative name for 'sector' column
#' @param dest_col alternative name for 'destination' column
#' @param use_col alternative name for 'use' column
#' @param flow_col alternative name for 'flow' column
#' @return Returns a list of
#'   \describe{
#'   \item{location_id}{vector of location ids}
#'   \item{sector_id}{vector of sector ids}
#'   \item{R}{matrix of location-sector revenues with
#'            \code{N = length(location_id)} rows and
#'            \code{J = length(sector_id) columns}}
#'   \item{D}{vector of locations' deficit transfers, i.e. positive for a trade
#'            deficit and negative for a trade surplus}
#'   \item{pi}{vector of import shares across locations for each
#'             sector-destination-use combination, where the import share of
#'             origin \code{o} in sector \code{s} products used in destination
#'             \code{d} and use category \code{u} is at position
#'             \code{o + (d-1)*N + (s-1)*N*N + (u-1)*N*N*J}}
#'   \item{alpha}{matrix of sectoral consumption shares in each location with
#'                \code{N} rows and \code{J} columns}
#'   \item{gamma_jrs}{vector of intermediate cost shares across sectrs for each
#'                    destination-use combination, where the cost share of
#'                    sector \code{s} intermediates used in destination
#'                    \code{d} and use category \code{u} is at position
#'                    \code{d + (s-1)*N + (u-1)*N*J}}
#'   \item{gamma_js}{matrix of location-sector labor cost shares with
#'                \code{N} rows and \code{J} columns}}
#' @export convert_io
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by ungroup mutate distinct n_distinct arrange
#' @importFrom dplyr summarise

convert_io <- function(iot,
                       orig_col = "origin",
                       sec_col = "sector",
                       dest_col = "destination",
                       use_col = "use",
                       flow_col = "flow") {

  iot <- iot %>%
    rename(origin = all_of(orig_col),
           sector = all_of(sec_col),
           destination = all_of(dest_col),
           use = all_of(use_col),
           flow = all_of(flow_col))

  J <- n_distinct(iot$origin)
  S <- n_distinct(iot$sector)

  locations <- sort(unique(iot$origin))
  sectors <- sort(unique(iot$sector))
  fc_name <- unique(iot$use)[!(unique(iot$use) %in% sectors)]

  # revenue matrix
  iot <- iot %>%
    group_by(origin, sector) %>%
    mutate(revenue = sum(flow)) %>%
    ungroup()

  R <- iot %>%
    distinct(origin, sector, revenue) %>%
    arrange(sector, origin) %>%
    .$revenue %>%
    matrix(nrow=J, ncol=S)

  # deficit transfers (positive for trade deficits)
  iot <- iot %>%
    group_by(destination) %>%
    mutate(agg_use = sum(flow)) %>%
    group_by(origin) %>%
    mutate(agg_revenue = sum(flow),
           agg_use = unique(agg_use[origin == destination]),
           deficit_transfer = agg_use - agg_revenue) %>%
    ungroup()

  D <- iot %>%
    distinct(origin, deficit_transfer) %>%
    arrange(origin) %>%
    .$deficit_transfer

  # sector-by-sector and sector-by-final-demand import shares
  iot <- iot %>%
    group_by(destination, sector, use) %>%
    mutate(imp_shr = flow / ifelse(sum(flow) == 0, 1, sum(flow))) %>%
    ungroup()

  pi <- iot %>%
    arrange(use, sector, destination, origin) %>%
    .$imp_shr

  # sectoral consumption shares (inclusive of potential tariffs)
  tmp <- iot %>%
    filter(use == fc_name) %>%
    group_by(destination, sector) %>%
    summarise(sec_consumption = sum(flow), .groups = "drop_last") %>%
    mutate(cons_sec_shr = sec_consumption / sum(sec_consumption)) %>%
    ungroup()

  alpha <- tmp %>%
    arrange(sector, destination) %>%
    .$cons_sec_shr %>%
    matrix(nrow = J, ncol = S)

  # intersectoral intermediate shares
  tmp <- iot %>%
    filter(use != fc_name) %>%
    group_by(destination, use) %>%
    mutate(imp_revenue = revenue[origin == destination & sector == use]) %>%
    group_by(destination, sector, use) %>%
    summarise(sec_interm = sum(flow), imp_revenue = unique(imp_revenue),
              .groups = "drop") %>%
    mutate(sec_use_shr = ifelse(imp_revenue == 0, 0, sec_interm / imp_revenue))

  gamma_jrs <- tmp %>%
    arrange(use, sector, destination) %>%
    .$sec_use_shr %>%
    array(dim = c(J, S, S))

  # vad shares
  gamma_js <- 1 - apply(gamma_jrs, 3, rowSums)
  gamma_jrs <- as.vector(gamma_jrs)

  return(list("location_id" = locations,
              "sector_id" = sectors,
              "R" = R,
              "D" = D,
              "pi" = pi,
              "alpha" = alpha,
              "gamma_jrs" = gamma_jrs,
              "gamma_js" = gamma_js))
}
