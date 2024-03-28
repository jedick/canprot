\name{NEWS}
\title{News for Package 'canprot'}
\encoding{UTF-8}

\section{Changes in canprot version 1.1.2-42 (2024-03-28)}{

    \itemize{

      \item Objects in \code{cplab} are now \code{quote()}-ed instead of
      \code{expression()}s for easier handling by other functions.

      \item Functions for chemical metrics (\code{ZCAA()}, \code{O2AA()},
      \code{H2OAA()}, and \code{MWAA()}) have been renamed to \code{Zc()},
      \code{nO2()}, \code{nH2O()}, and \code{MW()}.

      \item Add \code{terminal_H2O} argument to \code{nH2O()}, used to specify
      how terminal -H and -OH groups are counted (they are not by default). 

      \item Add \code{...} argument to functions to permit using the same
      arguments for all metrics in \code{get()} or \code{do.call()}
      constructions.

      \item Add functions for metrics: \code{plength()} (protein length),
      \code{V0()} and \code{S0()} (standard molal volume and entropy per
      residue), \code{pV0()} and \code{pS0()} (volume and entropy per protein),
      \code{V0g()} and \code{S0g()} (specific volume and entropy),
      \code{nO2g()} and \code{nH2Og()} (specific oxidation and hydration state),
      \code{Density()} (mass density), \code{SV()} (entropy density),
      \code{nC()} and \code{pnC()} (C atoms per residue and protein),
      \code{Cost()} (metabolic cost of amino acids from
      \href{https://doi.org/10.1073/pnas.062526999}{Akashi and Gojobori,
        2002}), \code{RespiratoryCost()} and \code{FermentativeCost()}
      (respiratory and fermentative costs from
      \href{https://doi.org/10.1093/molbev/msi126}{Wagner, 2005}),
      \code{B20Cost()}, \code{Y20Cost()}, and \code{H11Cost()} (biosynthetic
      costs from \href{https://doi.org/10.1038/s41467-018-06461-1}{Zhang et
        al., 2018}).

      \item Differential expression datasets (including old to new ID mappings
      in \file{uniprot_updates.csv}) and chemical analysis vignettes have been
      moved from \pkg{canprot} to to \pkg{JMDplots}.

      \item \code{read_fasta()}, \code{count_aa()}, and \code{sum_aa()} have
      been renamed and moved from \pkg{CHNOSZ} to \pkg{canprot}.

      \item \code{calc_metrics()} and \code{add_hull()} have been moved from
      \pkg{chem16S} to \pkg{canprot}.

      \item \code{protcomp()} has been renamed \code{human_aa()}.

      \item Add \code{add_cld()} (add compact letter display to boxplots).

      \item Add demos \file{thermophiles.R}, \file{locations.R}, and
      \file{redoxins.R}.

      \item Add vignettes \emph{Demos for canprot} and \emph{More about
        metrics}.

    }

}
