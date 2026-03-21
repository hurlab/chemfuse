test_that("cf_screen validates input", {
  expect_error(cf_screen(NULL))
  expect_error(cf_screen("nonexistent_file.csv"))
  expect_error(cf_screen(data.frame(x = 1:3)))  # no 'smiles' or 'name' column
})

test_that("cf_screen accepts a data frame with smiles column", {
  skip_if_no_python()

  df <- data.frame(smiles = c("CCO", "CCC"), stringsAsFactors = FALSE)
  result <- cf_screen(df)
  expect_s3_class(result, "tbl_df")
})

test_that("cf_screen accepts a CSV file path", {
  skip_if_no_python()

  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  df <- data.frame(smiles = c("CCO", "CCC"), stringsAsFactors = FALSE)
  write.csv(df, tmp, row.names = FALSE)

  result <- cf_screen(tmp)
  expect_s3_class(result, "tbl_df")
})

test_that("cf_screen handles data frame with name column", {
  skip_if_no_python()

  df <- data.frame(name = c("aspirin", "caffeine"), stringsAsFactors = FALSE)
  result <- cf_screen(df)
  expect_s3_class(result, "tbl_df")
})
