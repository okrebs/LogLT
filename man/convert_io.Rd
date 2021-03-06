% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_io.R
\name{convert_io}
\alias{convert_io}
\title{Convert standardized IO Table to simulation inputs}
\usage{
convert_io(
  iot,
  orig_col = "origin",
  sec_col = "sector",
  dest_col = "destination",
  use_col = "use",
  flow_col = "flow"
)
}
\arguments{
\item{iot}{An input output table in standardized long format (see the package
'iotr') with the columns, \code{origin}, \code{sector}, \code{destination},
\code{use} and \code{flow}. Alternative column names can be specified with
the following parameters}

\item{orig_col}{alternative name for 'origin' column}

\item{sec_col}{alternative name for 'sector' column}

\item{dest_col}{alternative name for 'destination' column}

\item{use_col}{alternative name for 'use' column}

\item{flow_col}{alternative name for 'flow' column}
}
\value{
Returns a list of
\describe{
\item{location_id}{vector of location ids}
\item{sector_id}{vector of sector ids}
\item{R}{matrix of location-sector revenues with
\code{N = length(location_id)} rows and
\code{J = length(sector_id) columns}}
\item{D}{vector of locations' deficit transfers, i.e. positive for a trade
deficit and negative for a trade surplus}
\item{pi}{vector of import shares across locations for each
sector-destination-use combination, where the import share of
origin \code{o} in sector \code{s} products used in destination
\code{d} and use category \code{u} is at position
\code{o + (d-1)*N + (s-1)*N*N + (u-1)*N*N*J}}
\item{alpha}{matrix of sectoral consumption shares in each location with
\code{N} rows and \code{J} columns}
\item{gamma_jrs}{vector of intermediate cost shares across sectrs for each
destination-use combination, where the cost share of
sector \code{s} intermediates used in destination
\code{d} and use category \code{u} is at position
\code{d + (s-1)*N + (u-1)*N*J}}
\item{gamma_js}{matrix of location-sector labor cost shares with
\code{N} rows and \code{J} columns}}
}
\description{
\code{convert_io} takes an input output table in standardized long format and
converts it to the necessary model inputs for simulations.
}
\details{
Our simulation requires data on country-sector revenues, trade
deficits, sectoral import shares by use category (final or sectoral
intermediate), country-sector consumption shares, labor cost shares and cost
shares of material inputs, all of which can be derived from an input output
table.
}
