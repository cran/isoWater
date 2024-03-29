#Prep MWL
O = runif(10, -15, -2)
H = O * 8 + 10 + rnorm(10, 0, 6)
HO = data.frame(H, O)

MWL = mwl(HO, plot = FALSE)
test_that("mwl works", {
   expect_length(MWL, 6)
})

d = dex(HO, form = "both")
test_that("dex works", {
  expect_s3_class(d, "data.frame")
  expect_equal(d[,3], H - 8 * O)
  expect_error(dex(H))
  expect_error(dex(HO, "lcex", 12))
  expect_error(dex(HO, "jim"))
})

obs = iso(-60, -6, 0.5, 0.1, 0)
slope = c(5, 0.3)

test_that("mwlSource works", {
  expect_length(mwlSource(obs, MWL, slope, ngens = 1e2), 2)
  expect_length(mwlSource(obs, MWL, slope, ngens = 1e2, ncores = 2), 2)
})

obs2 = obs
class(obs2) = "data.frame"
test_that("mwlSource warnings work", {
  expect_warning(mwlSource(obs, slope = slope, stype = 2, ngens = 1e2))
  expect_warning(mwlSource(obs2, slope = slope, ngens = 1e2))
})

test_that("mwlSource errors work", {
  expect_error(mwlSource(obs, slope = slope, stype = 4))
})
