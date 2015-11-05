passesBlacklist <- function(prot.desc, bl, perl = TRUE, ignore.case = TRUE) {
  all(as.logical(lapply(bl, function(regex) !grepl(regex, prot.desc, perl = perl, 
    ignore.case = ignore.case))))
}

filterProtDesc <- function(prot.desc, filtrs, perl = TRUE, ignore.case = TRUE) {
  res <- gsub("\\s{2,}", " ", Reduce(function(processed.prot.desc, curr.regex) {
    gsub(curr.regex, " ", processed.prot.desc, perl = perl, ignore.case = ignore.case)
  }, filtrs, init = prot.desc), perl = TRUE)
  sub("^\\s+", "", sub("\\s+$", "", res))
}

loadRegexList <- function(path.to.list, dump.regex = "(?i)", ...) {
  as.character(lapply(readLines(path.to.list), function(x) gsub(dump.regex, "", 
    x, ...)))
}

extractProtDescription <- function(fasta.line, regexs = getOption("ahrd.prot.desc.regex", 
  c("^\\S+\\s+", "\\s+[0-9A-Z]+\\s+OS=.+$")), perl = TRUE, ignore.case = TRUE) {
  filterProtDesc(fasta.line, regexs, perl = perl, ignore.case = ignore.case)
}

extractBlacklistAndFilter <- function(fasta.lines, bl, filtrs, perl = TRUE, ignore.case = TRUE, 
  lapply.funk = mclapply) {
  prot.descs <- as.character(lapply.funk(fasta.lines, extractProtDescription, 
    perl = perl, ignore.case = ignore.case))
  bl.res <- prot.descs[as.logical(lapply.funk(prot.descs, passesBlacklist, bl = bl, 
    perl = perl, ignore.case = ignore.case))]
  as.character(lapply.funk(bl.res, filterProtDesc, filtrs = filtrs, perl = perl, 
    ignore.case = ignore.case))
}

descriptionFrequencies <- function(prot.descs, lapply.funk = mclapply, set.name = "") {
  upd <- unique(prot.descs)
  freq.df <- data.frame(stringsAsFactors = FALSE, description = upd, frequencies = as.numeric(lapply.funk(upd, 
    function(x) length(prot.descs[which(prot.descs == x)]))))
  freq.df$relative.freqs <- freq.df$frequencies/length(prot.descs)
  freq.df$frequency.of.relative.freqs <- 0
  unq.rel.freqs <- sort(unique(freq.df$relative.freqs))
  lapply(unq.rel.freqs, function(x) {
    freq.df[which(freq.df$relative.freqs == x), "frequency.of.relative.freqs"] <<- nrow(freq.df[which(freq.df$relative.freqs == 
      x), ])/nrow(freq.df)
  })
  cbind(freq.df, set.name = set.name)
}

toWords <- function(prot.desc, split.regex = "-|/|;|\\\\|,|:|\"|'|\\.|\\s+|\\||\\(|\\)", 
  token.blacklist = NULL, ...) {
  if (!is.null(token.blacklist) && length(token.blacklist) > 0) {
    prot.desc <- filterProtDesc(prot.desc, token.blacklist, ...)
  }
  wrds <- strsplit(tolower(prot.desc), split = "\\s+|-|\\(|\\)", perl = TRUE)[[1]]
  wrds[wrds != ""]
}

fMeasure <- function(ref.wrds, asgnd.wrds, f.beta = getOption("ahrd.f.beta", 2)) {
  true.pos <- length(intersect(ref.wrds, asgnd.wrds))
  precision <- true.pos/length(asgnd.wrds)
  recall <- true.pos/length(ref.wrds)
  # Compute the F-Measure, weighted with f.beta
  f.scr <- (1 + f.beta^2) * (precision * recall)/((f.beta^2 * precision) + recall)
  if (is.na(f.scr)) {
    0
  } else {
    f.scr
  }
}

referenceSetSimilarities <- function(ref.set.a, ref.set.b, ref.set.a.freqs) {
  as.numeric(unlist(mclapply(1:length(ref.set.a), function(i) {
    prot.desc <- ref.set.a[[i]]
    curr.best.f.scr <- 0
    for (desc.2.comp in ref.set.b) {
      iter.f.scr <- fMeasure(prot.desc, desc.2.comp)
      if (iter.f.scr > curr.best.f.scr) curr.best.f.scr <- iter.f.scr
      if (curr.best.f.scr == 1) break
    }
    rep(curr.best.f.scr, ref.set.a.freqs[[i]])
  })))
}

overlapScore <- function(queryStart,queryEnd,
  queryLength,subjectStart,subjectEnd, subjectLength ) {
    ((queryEnd - queryStart + 1.0) + (subjectEnd - subjectStart + 1.0))
				/ (queryLength + subjectLength)
}
