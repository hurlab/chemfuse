test_that("cf_to_csv validates input", {
  expect_error(cf_to_csv(NULL, "output.csv"))
  expect_error(cf_to_csv(data.frame(x = 1), NULL))
  expect_error(cf_to_csv(data.frame(x = 1), ""))
  expect_error(cf_to_csv("not_a_df", "output.csv"))
})

test_that("cf_to_csv creates a file", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  df <- data.frame(name = c("aspirin", "caffeine"), mw = c(180.2, 194.2))
  result <- cf_to_csv(df, tmp)

  expect_true(file.exists(tmp))
  expect_equal(result, tmp)

  # Verify file content
  written <- utils::read.csv(tmp)
  expect_equal(nrow(written), 2L)
  expect_true("name" %in% names(written))
})

test_that("cf_to_csv returns path invisibly", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  df <- data.frame(x = 1:3)
  expect_invisible(cf_to_csv(df, tmp))
})

test_that("cf_to_excel validates input", {
  expect_error(cf_to_excel(NULL, "output.xlsx"))
  expect_error(cf_to_excel(data.frame(x = 1), NULL))
  expect_error(cf_to_excel("not_a_df", "output.xlsx"))
})

test_that("cf_to_excel creates a file when openxlsx is available", {
  skip_if_not_installed("openxlsx")

  tmp <- tempfile(fileext = ".xlsx")
  on.exit(unlink(tmp))

  df <- data.frame(name = c("aspirin", "caffeine"), mw = c(180.2, 194.2))
  result <- cf_to_excel(df, tmp)

  expect_true(file.exists(tmp))
  expect_equal(result, tmp)
})

test_that("cf_to_excel errors without openxlsx", {
  skip_if(requireNamespace("openxlsx", quietly = TRUE), "openxlsx is installed")

  df <- data.frame(x = 1:3)
  expect_error(cf_to_excel(df, "output.xlsx"), regexp = "openxlsx")
})
