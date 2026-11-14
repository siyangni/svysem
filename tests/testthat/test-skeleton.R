test_that("the package namespace loads", {
  expect_true(is.function(svysem))
  expect_true(is.function(svycfa))
  expect_true(is.function(svygrowth))
})
