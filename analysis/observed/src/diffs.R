## calculate the difference betweent the early and late stages of
## assembly

getDiff <- function(dats, metrics){
  early <- dats[dats$assem == "early",]
  late <- dats[dats$assem == "late",]
  not.shared <- c(early$GenusSpecies[!early$GenusSpecies %in%
                                      late$GenusSpecies],
                  late$GenusSpecies[!late$GenusSpecies %in%
                                     early$GenusSpecies])
  early <- early[!early$GenusSpecies %in% not.shared,]
  late <- late[!late$GenusSpecies %in% not.shared,]
  late <- late[match(early$GenusSpecies, late$GenusSpecies),]
  diffs <- late[, metrics] - early[, metrics]
  colnames(diffs) <- paste("diff", metrics, sep=".")
  diffs$GenusSpecies <- early$GenusSpecies
  out <- merge(early, diffs,
               by = "GenusSpecies", all.x=TRUE)
  out <- merge(out, late[,c("GenusSpecies", metrics)],
               by = "GenusSpecies", all.x=TRUE)
  return(out)
}


## specialization models
specMods <- function(metric, diff.dats, type){

  y <- diff.dats[, paste(metric, "y", sep=".")]
  x <- diff.dats[, paste(metric, "x", sep=".")]

  out.mod <- lmer(y ~ x +
                    (1|Site) +
                    (1|GenusSpecies),
                  data=diff.dats)

  out.mod2 <- lmer(y ~ 1 + offset(x) +
                     (1|Site) +
                     (1|GenusSpecies),
                   data=diff.dats)
  out.test <- anova(out.mod, out.mod2)
  save(out.mod, file=sprintf("saved/mods/%s%s.Rdata", type, metric))
  return(out.test)
}
