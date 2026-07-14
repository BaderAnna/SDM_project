# =============================================================================
# 1 - TSS
# =============================================================================

# Install mecofun from Gitlab using the devtools package:
library(devtools)
devtools::install_git("https://gitup.uni-potsdam.de/macroecology/mecofun.git")

# citation
paste("Zurell, D. (2024).",
      "mecofun: useful functions for macroecology and species distribution modelling",
      "version 0.7.1. University of Potsdam, Potsdam.",
      "https://gitup.uni-potsdam.de/macroecology/mecofun")


library(mecofun)


# threshold
optimal.thresholds(DATA = NULL, threshold = 101, which.model = 1:(ncol(DATA)-2), 
                   model.names = NULL, na.rm = FALSE, opt.methods = NULL, req.sens, req.spec, 
                   obs.prev = ObsPrev, smoothing = 1, FPC, FNC)

data(Anguilla_train)
m1 <- glm(Angaus ~ poly(SegSumT,2), data=Anguilla_train, family='binomial')
preds_cv <- crossvalSDM(m1, kfold=5, traindat=Anguilla_train, colname_species = 'Angaus', colname_pred = 'SegSumT')
evalSDM(Anguilla_train$Angaus, preds_cv)



