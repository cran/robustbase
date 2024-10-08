% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/robustbase/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{robustbase}}
\encoding{UTF-8}

  % --- 1.0-0  once we fix  adjOutlyingness() properly !
  %     =====

\section{CHANGES in robustbase VERSION 0.99-4-1 (2024-09-25, svn 1006)}{
  \subsection{Misc}{
    \itemize{
      \item small changes in \file{tests/lmrob-psifns.R} and \file{*.Rout.save}, notably
      avoiding small diffs on "M1mac" and "noLD" (\command{--disable-long-double}).
    }
  }
}

\section{CHANGES in robustbase VERSION 0.99-4 (2024-08-19, svn 1002)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{lmrob.fit()} gets new option \code{bare.only} for some
      speed up, in case a relatively \dQuote{bare} result is sufficient.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item fix another LaTeX no more wanting math in section
      titles (same problem as 0.99-3), now for \file{lmrob_simulation.Rnw}.
      Now also use \verb{\texorpdfstring} inside \verb{\title{}} of
      \file{psi_functions.Rnw} \emph{and} use \verb{newtheorem{..}}
      \emph{after} \verb{\hypersetup{...}} (thanks to Kurt).

      \item avoid some \pkg{ggplot2} warnings (about omitting obs.)

      \item do not warn about skipping equilibration for p=1.

      \item define \code{STRICT_R_HEADERS} now (src/*.[ch]).
    }
  }
}

\section{CHANGES in robustbase VERSION 0.99-3 (2024-07-01, svn r992)}{
  \subsection{Misc}{
    \itemize{
      \item Again adapt \file{../vignettes/lmrob_simulation.Rnw} to
      \pkg{ggplot2} changes (\code{ s/ fun.y / fun / }).

      \item \file{tests/subsample.R} now shows "Matrix" info (as that gave
      confusing error during Matrix / R-devel changes).

      \item Latest LaTeX does not like math in \verb{\section{.}}; hence
      replacing \verb{\psi} by \verb{\texorpdfstring{$\psi$}{psi}}.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.99-2 (2024-01-27, svn r989)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item fix thinko bug introduced exact fit checking in 0.99-0; thanks
      to the report by Thomas Mang.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.99-1 (2023-11-28, svn r988)}{
  \subsection{Misc}{
    \itemize{
      \item document the \code{...} argument(s) of \code{nlrob.control()}.
      \item \file{DESCRIPTION}: more URLSs.
      \item \file{*.Rd} tweaks; avoiding other \verb{NOTE}s.
      \item \file{src/lmrob.c,mc.c,monitor.c}: format \%lld etc
    }
  }
}

\section{CHANGES in robustbase VERSION 0.99-0 (2023-06-12, svn r928)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item In case it finds an \emph{exact fit}, i.e., \eqn{k > n/2}
      residuals are (practically) zero, \code{lmrob.S()} now returns that,
      \code{scale = 0}, and robustness weights \code{weights(.,	"robustness")}
      which are either 0 or 1 quickly and more consistently.

      This is \emph{not yet} satisfactorily working for the \dQuote{large n} case.

      **** \bold{TODO}: Optionally (but off by default?? for back compatibility??),
      \code{lmrob()} initialized by such an exact-fit estimator may use a
      (user-specified or automatic) scale \eqn{\hat{s} > 0}{s^ > 0} to
      continue providing efficient SM / MM estimates.

      \item \code{lmrob.control()} gets class \code{"lmrobCtrl"} which can be
      \code{print()}ed and \code{update()}d.

      \item New tuning const \code{zero.scale = 1e-10} instead of hardwired
      \code{EPS_ZERO} in \code{lmrob.S()}'s C code; additionally the exact
      fit / zero-residuals checking is finally y-scale equivariant, using
      \eqn{1/n ||y||_1 = } \code{mean(abs(y))}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed \code{covMcd()}: corrected the consistency correction
      factor for the reweighted estimate: was MCDCons(p, sum(w)/n),
      now is MCDCons(p, 0.975), see Croux and Haesbroeck (1999).

      Consequently, often \code{covMcd()} and estimates based on it, such
      as \code{BYlogreg()} or \code{ltsReg()} are slightly changed. % 3 of our ex. checks broke

      \item registered \code{residuals.lmrob.S} \emph{and} changed
      \code{formals} compatible with generic to \code{(object, ...)}.

      \item \code{lmrob.S()} in C level \code{refine_fast_s()} no longer
      modifies \code{beta_cand[]} which notably fixes a bug found in
      exact-fit examples previously resulting in \code{coef() == 0}
      \dQuote{results}, also for \code{lmrob()} which calls
      \code{lmrob.S()} for its default initial estimate.

      \item \code{print(<lmrob.S>)} finally works again as intended
      (\emph{not} printing a largish list); mostly by fixing internal
      \code{lmrob.control.minimal()}.

      \item \code{summary.lmrob(obj)} now also gets a \code{weights}
      component from \code{obj} if that is non-trivial; this may also
      correct \code{summary()} output in weighted cases.

      \item Minor fixes in \code{plot.lts()} related to par mfrow.

      \item ltsReg() Help file changed - the note strongly advising NOT
      to use LTS regression removed.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item argument \code{mf} e.g., in \code{lmrob.fit()} has been deprecated
      since 2017-12-09 and is defunct now.

      \item \code{lmrob()} tweaks, notably saying more for \code{trace.lev >= 3}.

      \item replaceed (almost) all \code{Calloc()} by \code{R_alloc()} in C code;
      e.g., prevents leakage from user interrupt.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.95-1 (2023-03-28, svn r903)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed \code{covMcd()}: \code{raw.only=TRUE} does not work in the case \code{nsamp="deterministic"}
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item New \code{lmc()} and \code{rmc()} for robust tail weight
      (\dQuote{robust kurtosis}). % ../man/lmc-rmc.Rd
    }
  }
  \subsection{Misc}{
    \itemize{
      \item Our \file{tests/*.R} no longer rely on the \CRANpkg{Matrix}
      \file{test-tools.R} collection but rather use our much smaller own
      one in \file{xtraR/test-tools.R}.

      \item fixed \code{\\} escape in \file{man/aircraft.Rd}.

      \item initialize \var{equed} in \file{src/lmrob.c} to avoid
      false positive compiler warning.

      \item document "internal" \code{.regularize.Mpsi()}.

      \item replace \code{()} by \code{(void)} arguments for \code{F77_SUB(println)}.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.95-0 (2022-04-02, svn r894)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{huberM(x, k)} now also works in small sample cases with
      \code{+/-Inf}, e.g. for \code{x = c(-Inf, 1)}.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item New \code{huberize()} utility, notably by default used in \code{mc()}.

      \item \code{mc()} now by default uses \code{doScale=FALSE} in order
      to guarantee convergence also in extreme cases when used \emph{together}
      with new \code{c.huberize = 1e11}.  Both are \emph{not} back
      compatible but clearly better in extreme cases, as assessed by
      Lukas Graz' BSc thesis.
      Note that the behavior of the new default \code{doScale=FALSE} has
      been amended in the C code to be more compatible in its \code{eps}
      checks with (the previous default) \code{doScale=TRUE}.

      Lastly, a \code{message()} is printed currently, at most once per \R session,
      whenever \code{mc()} is called without an explicit specification of \code{doScale=*},
      alerting users to the changed default.  This message can be turned
      off by setting \code{options(mc_doScale_quiet = TRUE)}.

      \item New dataset \code{data(x30o50)} as an example of potentially
      problematic data, not just for \code{mc()}.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item rather "truncating" \eqn{+/-}\code{Inf} to <large>, use
      \code{NA.OK = TRUE} in \code{Qn()} and \code{Sn()}'s \code{.C()} call.

      \item Our C code no longer uses the \code{DOUBLE_*} constants from S,
      but rather standard \code{<float.h>}.
    }
  }
}


\section{CHANGES in robustbase VERSION 0.93-9 (2021-09-27, svn r888)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{Qn(x)} and \code{Sn(x)} gain an \code{na.rm = FALSE}
      option and now work when \code{x} has \code{NA}s and contains
      \eqn{\pm}{+/-} \code{Inf}.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item No longer include \file{Rdefines.h} as it is somewhat deprecated.

      \item In \file{src/lmrob.c} add \code{USE_FC_LEN_T} and use \code{FCONE}
      rather than the older FCLEN approach.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.93-8 (2021-06-01, svn r879)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{scaleTau2()} gets new optional \code{iter = 1} and
      \code{tol.iter} arguments; mostly experimentally to see if or when
      iteration makes sense.

      \item \code{Qn(x, *)} gets new optional \code{k = .} to indicate the
      \dQuote{quantile} i.e., order statistic to be computed (with default as
      previously hard-coded).\cr
      Experimentally to try for cases where more than n/2 observations
      coincide (with the median), i.e., \code{x[i] == x0 == median(x[])},
      and hence \code{Qn(x)} and \code{mad(x)} are zero.

      \item \code{adjOutlyingness()} gets new option \code{IQRtype = 7}.
    }
  }
  \subsection{Tweaks}{
    \itemize{
      \item For tests: *again* differences found in the non-sensical
      \code{adjOutlyingness()} example (with large p/n, hence many "random" values
      in the order of 1e15).
      Disable the test for now (and record the result in *.Rout).
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item The \code{test()} utility in \file{tests/lmrob-ex12.R} no
      longer calls \code{matrix(x, n,4)} where the length of x does not
      match \code{4n}.   Similar change in \file{tests/mc-strict.R}
    }
  }
}

\section{CHANGES in robustbase VERSION 0.93-7 (2021-01-03, svn r865)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Use \command{\\CRANpkg\{.\}} in most places, providing web
      links to the respective CRAN package page.

      \item \code{adjOutlyingness()} now gains optional parameters to be
      passed to \code{mc()}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item update the internal man page, so new \code{checkRdContents()}
      is happy.
      \item fix several \samp{\\url{.}}'s that now are diagnosed as
      \sQuote{moved}.

      \item \code{adjOutlyingness()} finally works with \code{p.samp > p}.

      \item \code{scaleTau2()} now works with \code{Inf} and very large
      values, and obeys new \code{na.rm = FALSE} argument.

      \item add \code{check.environment=FALSE} to some of the
      \code{all.equal()} calls (for 'R-devel', i.e., future R 4.1.x).

      \item \code{wgt.himedian(numeric())} now returns \code{NA} instead of
      occasionally seg.faulting or inf.looping.  Ditto for a case when
      called from \code{Qn()}.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.93-6 (2020-03-20, svn r854)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{splitFrame()} now treats \code{character} columns also as
      categorical (the same as \code{factor}s).
    }
  }
  \subsection{Tweaks}{
    \itemize{
      \item Small updates, also in checks for newer compiler settings,
      e.g., \code{FCLEN} macro; also F77_*() etc, in order to fix 'LTO'
      issues.

      \item More careful or \emph{less} calling \code{intpr()}: correct "Rank" of array
      (for gfortran/gcc 10, when \code{-fallow-argument-mismatch} is not set).
    }
  }
}

\section{CHANGES in robustbase VERSION 0.93-5 (2019-05-09, svn r840)}{
  \subsection{Tweaks}{
    \itemize{
      \item For tests: try harder in \code{moreSessionInfo()} to detect
      openBLAS (should work for KH's debian openBLAS setup).
    }
  }
}

\section{CHANGES in robustbase VERSION 0.93-4 (2019-03-18, svn r839)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \env{R_CHECK_LENGTH_1*} fix in default value for \code{ask} in
      plot methods for \code{mcd} and \code{lts}

      \item \file{src/rowMedians.c} gets better PROTECT()ion thanks to
      Tomas Kalibera's "rchk"s.
    }
  }
  \subsection{Tweaks}{
    \itemize{
      \item losen the strict tests more, convergence issues on all
      non-standard BLAS platforms.  Now modularized the partly platform
      dependent testing into \code{moreSessionInfo()} in
      \file{./xtraR/platform-sessionInfo.R}.  Further, for now use
      \code{RNGversion("3.5.0")} e.g. in \file{../tests/mc-strict.R}.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.93-3 (2018-09-06, svn r827)}{
  \subsection{Tweaks}{
    \itemize{
      \item losen the strict tests, notably in
      \file{tests/nlregrob-tst.R} to work around \code{nls()}
      convergence issues on all non-standard BLAS platforms.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.93-2 (2018-07-26, svn r825)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Update \code{lmrob(<empty multivariate>)} to the change in R
      (incl 3.5.1 patched).

      \item \code{mc()} gets new optional \code{doScale} argument, and
      \emph{increased} defaults for the tolerances \preformatted{    eps1 = 1e-14, eps2 = 1e-15}
      such that it should converge by default in more cases.

      \item A \code{na.action} is now kept in \code{summary(lmrob(*))},
      and when \code{print()}ing it, a note about omitted observations,
      e.g., because of \code{NA}'s, is made as for \code{lm()}.

      \item Internal \code{lmrob.weights()}: more "resistant" in case
      scale=0, using na.rm=TRUE  (report only, no reprex).

      \item \code{lmrob(*, trace.lev >= 2)} now shows some information
      about the number of \code{find_scale()} iterations used (as these
      are now stored C internally).% TODO: return the maximum used to R?
    }
  }
}

\section{CHANGES in robustbase VERSION 0.93-1-1 (2018-07-18, svn r818)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \file{src/robustbase.h}: \code{is_redescender} now is
      \code{static inline}, needed for some compilers, e.g., on
      ubuntu 18.04.
      Fixing R-forge bug(s) 6588 (and 6590, 6593),
      \url{https://r-forge.r-project.org/tracker/index.php?func=detail&aid=6588&group_id=59&atid=302}
    }
  }
}

\section{CHANGES in robustbase VERSION 0.93-1 (2018-06-20, svn r815)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The \sQuote{Usage:}s in the data set help pages now say
      \code{data(<..>, package="robustbase")}.

      \item The \file{lmrob_simulation} vignette now should continue
      to work with upcoming package \CRANpkg{ggplot2}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{Mpsi(x, c, psi="huber", deriv=-1)} now gives rho(x)
      instead of mostly \code{Inf}.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.93-0 (2018-04-21, svn r805)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{.psi.const(*, "lqq")} now also gives a
      \code{"constants"} attribute.

      \item more examples and help on Mpsi() etc functions and tuning constants.

      \item The S estimator \code{lmrob.S()} and M-S estimator now both
      make use of the new \code{lmrob.control()} argument
      \code{scale.tol} which defaults to \code{1e-10}, its formerly
      hardwired value.

      \item \code{lmrob.S()} further gets a new option \code{only.scale = FALSE},
      which when true allows to only compute the S scale estimate.  In
      that case, but also generally, \code{trace.lev = 3} or larger also
      produces output showing the C level \code{find_scale()} iterations.

      \item (By Manuel Koller) There's now a small C API to call our
      \code{Mpsi()} etc from C code in other packages, as
      \code{C_psi()}, etc; using new \file{../inst/include/robustbase.h}.

      \item \code{nlrob()$call$algorithm} now always contains the
      algorithm used as a \code{\link{character}} string, compatibly
      with \code{\link{nls}()}.

      \item new data set \code{steamUse}.
      %% which also depicts the \dQuote{non-robust} behavior of \code{lmrob(setting="KS2014")}

      \item Vignette \file{lmrob_simulation.Rnw}: fixed the wrong
      \dQuote{emprical power} plots; with faster \CRANpkg{ggplot2}, remove
      all \code{eval=FALSE} for plots and longer store the
      \file{*.pdf}s.

      \code{nlrob()} gets \code{model} option to ask for the
      \code{model.frame} to be returned.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{lmrob(..., method = "S")} no longer necessarily
      produces a warning in \code{.vcov.w()}.

      \item \code{nlrob()} returns a correct \code{dataClasses} component.
    }
  }
  \subsection{Tests Tweaks}{
    \itemize{
      \item For use in non-R-internal BLAS/Lapack libraries, several
      \file{tests/*.R} examples have been tweaked.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.92-8 (2017-10-30, svn r778)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{fullRank()} utility for \code{adjOutlyingness}:

      \item \code{adjOutlyingness()}, gets new options \code{p.samp} and
      \code{trace.lev}, and when it fails to find enough good
      directions, now checks the rank and mentions \code{fullRank()} in
      case the matrix is (QR-)rank deficient.

      \item The \code{"lmrob"} method for \code{vcov()} gets optional
      argument \code{complete = TRUE}, where only \emph{complete =
	FALSE} is back compatible.

      \item improved (error) messages in singular case in
      \code{.vcov.avar1()}.

      \item \code{.psi.const()} is exported as well, and help is
      improved about using and setting non-default psi tuning constants.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item loosened some regression test tolerances (for alternatives
      to BLAS) in \file{tests/(mc-strict|poisson-ex)}.

      \item \code{scaleTau2(x, *)} gains an optional \code{sigma0}
      argument and now returns zero instead of \code{NaN} when
      \code{sigma0 == 0}.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.92-7 (2016-11-23, svn r742)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{vignette(psi_functions)}

      \item \code{adjOutlyingness(..., maxit.mult = max(100, p))} new option,
      allowing more iterations for finding \dQuote{good} projection directions.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{summary(nlrob(*, method = "CM"))} now works.
      \item \code{lmrob..MM..fit()} now works again when \code{x} and
      \code{y} are specified, but neither \code{method} nor \code{obj} is.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.92-6 (2016-05-28, svn r717)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Now provide an \code{"lmrob"} method for the standard \R
      generic function \code{hatvalues()}, and also export its lower
      level workhorse \code{.lmrob.hat()} (formerly hidden
      \code{lmrob.leverages()}), which now by default has
      \code{names(.)}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{.lmrob.hat()} (formerly \code{lmrob.leverages()}) has
      been corrected for the rank-deficient case.
      \item \code{classPC(m)} now also works for a 1-column matrix.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.92-5 (2015-07-21, svn r698)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \emph{Hidden} \code{print()} methods \code{print.summary.lmrob()}
      and \code{print.lmrob.S()} get a \code{showAlgo = TRUE} argument
      which the user can set to FALSE in order to suppress printing of
      the \dQuote{Algorithmic parameters}.
      \item import (remaining parts) from "base" packages.
      \item \code{summary(<nlrob>)} now also prints a summary on the residuals.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{summary(lmrob(.))}'s variance-covariance matrix is now
      called \code{cov} instead of \code{cov.unscaled} (because it
      \emph{is} scaled).  Code which has been using \code{vcov(<lmrob>)}
      or \code{<lmrob> $ cov}, or even \code{<summary.lmrob> $ cov} is not affected.
    }
  }
}

\section{CHANGES in robustbase VERSION 0.92-4 (2015-06-04, svn r696)}{
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
