test_that("slr works", {
  set.seed(1)
  HIV <- load_data() # Load HIV data
  X <- HIV[,1:60]
  X.adjusted <- sweep(X+1,rowSums(X+1),MARGIN = 1, FUN='/')
  y <- ifelse(HIV[,62] == "Pos", 1, 0)

  bp <- matrix(c(1, -1, -1, 1))
  rownames(bp) <- c('g_Bacteroides','g_RC9_gut_group', 'f_vadinBB60_g_unclassified','f_Erysipelotrichaceae_g_unclassified')
  colnames(bp) <- 'z1'

  expect_equal(slr(X.adjusted, y, method ='PC', family = 'binomial', threshold = 0.998)$sbp, bp)
})
