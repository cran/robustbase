% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/robustbase/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{robustbase}}
\encoding{UTF-8}

\section{CHANGES in robustbase VERSION 0.92-4 (2015-06-04)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Started this \file{NEWS.Rd} file, to eventually replace the \file{ChangeLog}
      \item \code{plot.lmrob()} also identifies largest residuals as
      \code{plot.lm()}.  Also gets new argument \code{panel}, and
      \code{add.smooth=TRUE} behavior.
      \item adapt to the fact that R 3.3.0 will have its own
      \code{sigma()} S3 generic.
      \item setup for having message translations (volunteers sought!).
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item more careful in \file{../src/mc.c} (valgrind, thanks to Brian)
      \item add missing documentation, better examples for \code{predict.lmrob}
      \item \code{warn.limit.*} checks in \code{lmrob*()}
    }
  }
}

\section{CHANGES in robustbase VERSION 0.92-3 (2015-01-14, svn r679)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The \sQuote{Co-Median} \code{covComed()} from Maria Anna,
      tweaked by Valentin and modified considerably by Martin.
      \item Also document (and export) \code{r6pack()} utility.
      \item New \code{smoothWgt()} function --- \dQuote{Biweight on a
	Stick} --- to be used as \code{wgtFUN} option for \code{covMcd()}
      or \code{covComed()}.
      \item New utility \code{colMedians()} and \code{rowMedians}, as we
      use columnwise medians in so many places.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Tweaks to \code{medcouple()}, after detecting flaws -- which
      may be inherent and have \emph{not} been removed.
      \item Improved its documentation and the \code{adjOutlyingness()}
      one, notably its \dQuote{central} case.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.92-2 (2014-11-22, svn r660)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{covMcd()} with new options (\code{kmini},
      \code{nmini}) now ok (sometimes wrong in 0.92-1).
    }
  }
}

\section{CHANGES in robustbase VERSION 0.92-1 (2014-11-18)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The deterministic MCD, via \code{covMcd(..., nsamp="deterministic")}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{adjOutlyingness()}: reverse the defaults of
      \code{clower} and \code{cupper} and fix an \dQuote{eternal}
      erronous \eqn{\pm}{+/-} swap; see new note in \file{man/adjOutlyingness.Rd}.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.92-0 (2014-11-18)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{nlrob()} now works with indexed vector parameters.
      \item new \code{outlierStats()} (Manuel).
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item got rid of Fortran compiler warnings about \bold{ancient} style.
      \item \code{nlrob(*, weigths)}, fixing R-forge bug #5988.
      \item \code{covMcd()} fix for \dQuote{MAD = 0} case (new
      \code{exactfit} code 3).
    }
  }
}

\section{CHANGES in robustbase VERSION 0.91-1 (2014-05-01)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}

\section{CHANGES in robustbase VERSION 0.91-0 (2014-04-24)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}

\section{CHANGES in robustbase VERSION 0.90-1 (2014-01-30)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}

\section{CHANGES in robustbase VERSION 0.9-8 (2013-06-14)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}