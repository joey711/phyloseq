# Run from the repository root:
# Rscript --vanilla .github/diagnostics/benchmark-composition.R [baseline.R] [result.json]
# The optional baseline is a snapshot of R/composition-methods.R from before
# an optimization. No baseline/current code is installed or written to Git.
args <- commandArgs(trailingOnly=TRUE)
baseline <- if (length(args)) args[1] else NULL
output <- if (length(args) >= 2L) args[2] else NULL
pkgload::load_all(".", quiet=TRUE, export_all=FALSE)

load_implementation <- function(file) {
  env <- new.env(parent=asNamespace("phyloseq"))
  sys.source(file, envir=env)
  env
}
implementations <- list(after=load_implementation("R/composition-methods.R"))
if (!is.null(baseline)) {
  implementations <- c(list(before=load_implementation(baseline)), implementations)
}

allocated_bytes <- function(fun) {
  if (!isTRUE(capabilities("profmem"))) return(NA_real_)
  path <- tempfile("phyloseq-allocation-")
  utils::Rprofmem(path)
  tryCatch(invisible(fun()), finally=utils::Rprofmem(NULL))
  lines <- readLines(path, warn=FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  sum(bytes, na.rm=TRUE)
}

same_result <- function(left, right) {
  if (inherits(left, "SummarizedExperiment")) {
    stopifnot(identical(class(left), class(right)), identical(dimnames(left), dimnames(right)),
              identical(SummarizedExperiment::assayNames(left), SummarizedExperiment::assayNames(right)))
    for (assay in SummarizedExperiment::assayNames(left)) {
      stopifnot(isTRUE(all.equal(SummarizedExperiment::assay(left, assay),
                                SummarizedExperiment::assay(right, assay), tolerance=1e-12)))
    }
    # Hash representations may be versioned; scientific parameters must agree.
    a <- S4Vectors::metadata(left)$phyloseqV2$transformations
    b <- S4Vectors::metadata(right)$phyloseqV2$transformations
    stopifnot(identical(names(a), names(b)))
    for (name in names(a)) {
      for (field in c("source_hashes", "result_hashes", "fingerprint_version")) {
        a[[name]][[field]] <- b[[name]][[field]] <- NULL
      }
    }
    stopifnot(identical(a, b))
  } else {
    stopifnot(isTRUE(all.equal(left, right, tolerance=1e-12)))
  }
}

results <- list()
for (shape in list(c(20000L, 100L), c(4000L, 500L))) {
  values <- matrix(as.double(rep(1:17, length.out=prod(shape))), shape[1], shape[2],
                    dimnames=list(paste0("t", seq_len(shape[1])), paste0("s", seq_len(shape[2]))))
  input <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=values))
  original <- input
  prepared <- lapply(implementations, function(env) {
    list(relative=env$transform_composition(input, "relative"),
         clr=env$transform_composition(input, "clr"))
  })
  cases <- list(
    check_counts=function(env, state) env$check_composition(input, scale="counts"),
    relative=function(env, state) env$transform_composition(input, "relative"),
    hellinger=function(env, state) env$transform_composition(input, "hellinger"),
    clr=function(env, state) env$transform_composition(input, "clr"),
    check_clr=function(env, state) env$check_composition(state$clr, assay="clr", scale="clr"),
    clr_from_relative=function(env, state) env$transform_composition(
      state$relative, "clr", assay="relative", input_scale="proportion"))
  for (case in names(cases)) {
    functions <- lapply(names(implementations), function(name) {
      env <- implementations[[name]]
      state <- prepared[[name]]
      action <- cases[[case]]
      function() action(env, state)
    })
    names(functions) <- names(implementations)
    reference <- lapply(functions, function(fun) fun())
    if (length(reference) == 2L) same_result(reference[[1]], reference[[2]])
    rm(reference)
    times <- stats::setNames(lapply(functions, function(fun) numeric(5)), names(functions))
    for (iteration in 1:5) {
      order <- if (iteration %% 2L) names(functions) else rev(names(functions))
      for (name in order) {
        gc()
        times[[name]][iteration] <- system.time(invisible(functions[[name]]()), gcFirst=FALSE)[["elapsed"]]
      }
    }
    for (name in names(functions)) {
      row <- list(taxa=shape[1], samples=shape[2], operation=case, implementation=name,
                   median_seconds=stats::median(times[[name]]), elapsed_seconds=times[[name]],
                   allocated_bytes=allocated_bytes(functions[[name]]))
      results[[length(results) + 1L]] <- row
      cat(sprintf("%d x %d | %-18s | %-6s | %.3f s | %.1f MiB allocated\n",
                   shape[1], shape[2], case, name, row$median_seconds,
                   row$allocated_bytes / 1024^2))
    }
    stopifnot(identical(input, original))
  }
}
report <- list(
  benchmark="composition optimization", date=as.character(Sys.Date()),
  R=R.version.string, platform=R.version$platform,
  baseline_md5=if (!is.null(baseline)) unname(tools::md5sum(baseline)) else NULL,
  current_md5=unname(tools::md5sum("R/composition-methods.R")),
  repetitions=5L, allocation="Total bytes allocated per call; not peak resident memory",
  correctness="Assays, IDs, diagnostic reports and scientific parameters compared before timing",
  packages=as.list(vapply(c("SummarizedExperiment", "digest"), function(x) {
    as.character(utils::packageVersion(x))
  }, character(1))), results=results)
if (!is.null(output)) {
  dir.create(dirname(output), recursive=TRUE, showWarnings=FALSE)
  writeLines(jsonlite::toJSON(report, auto_unbox=TRUE, pretty=TRUE), output)
  cat("Results:", normalizePath(output), "\n")
}
