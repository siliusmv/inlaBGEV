test_that("standardise_vec() works on a simple vector", {
  x = rnorm(10)
  y = (x - mean(x)) / sd(x)
  expect_equal(standardise_vec(x), y)
})

test_that("standardise_vec() returns error for input with no length", {
  expect_error(standardise_vec(numeric(0)))
})

test_that("NA and NaN does not trouble standardise_vec()", {
  x = c(NA, rnorm(10), NaN)
  y = (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  expect_equal(standardise_vec(x), y)
})

test_that("standardise_vec() gives warning when x has sd of 0", {
  expect_warning(standardise_vec(rep(3, 3)))
})

test_that("standardise_vec() gives warning on logical variable", {
  expect_warning(standardise_vec(c(T, F, T)))
})

test_that("standardise() works on data.frames()", {
  x = rnorm(10)
  df1 = data.frame(x)
  df2 = df1 %>% dplyr::mutate(x = (x - mean(x)) / sd(x))
  expect_equal(standardise(df1), df2)
})

test_that("standardise() behaves nicely on logical df", {
  df = data.frame(x = c(T, F, T))
  expect_identical(df, standardise(df))
})

test_that("standardise_vec() works with external stats", {
  x = rnorm(10)
  y = rnorm(10)
  z = (x - mean(y)) / sd(y)
  expect_equal(z, standardise_vec(x, c(mean(y), sd(y))))
})

test_that("standardise_vec() gives error if stats don't have lengt 2", {
  expect_error(standardise_vec(rnorm(10), 3))
  expect_error(standardise_vec(rnorm(10), rep(1, 3)))
})

test_that("standardise_vec() gives warning if external stats have σ = 0", {
  expect_warning(standardise_vec(rnorm(10), c(1, 0)))
})

test_that("standardise_vec() gives error if external stats are not numerical", {
  expect_error(standardise_vec(rnorm(10), c(T, F)))
})

test_that("standardise_vec() gives error if external stats are NA/NaN", {
  expect_error(standardise_vec(rnorm(10), c(NA, NA)))
  expect_error(standardise_vec(rnorm(10), c(NaN, NaN)))
})

test_that("get_mean_and_sd_stats() compute correct stats", {
  df = data.frame(x = rnorm(10), y = rnorm(10))
  mat = apply(df, 2, function(x) c(mean(x), sd(x)))
  rownames(mat) = c("mean", "sd")
  expect_equal(mat, get_mean_and_sd_stats(df))
})

test_that("standardise() works with get_mean_and_sd_stats()", {
  df = data.frame(x = rnorm(10), y = rnorm(10), z = rnorm(10))
  stats1 = get_mean_and_sd_stats(df)
  df1 = standardise(df, stats1)
  stats2 = get_mean_and_sd_stats(df, c("x", "z"))
  df2 = standardise(df, stats2)
  df3 = standardise(df)
  expect_equal(standardise(df), df1)
  expect_false(identical(df1, df2))
})
