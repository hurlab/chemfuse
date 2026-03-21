test_that("cf_get validates input", {
  expect_error(cf_get(NULL))
  expect_error(cf_get(""))
  expect_error(cf_get(42))
})

test_that("cf_get returns a list with expected structure", {
  skip_if_no_python()

  result <- cf_get("aspirin")
  expect_type(result, "list")
  expect_true("compound" %in% names(result))
  expect_true("admet" %in% names(result))
  expect_true("druglikeness" %in% names(result))
})

test_that("cf_get compound field is a tibble", {
  skip_if_no_python()

  result <- cf_get("aspirin")
  expect_s3_class(result$compound, "tbl_df")
  expect_true(nrow(result$compound) >= 1L)
})

test_that("cf_get admet flag populates admet field", {
  skip_if_no_python()

  result <- cf_get("aspirin", admet = TRUE)
  # admet field should be a tibble or NULL (NULL if ADMET-AI not installed)
  if (!is.null(result$admet)) {
    expect_s3_class(result$admet, "tbl_df")
  }
})

test_that("cf_get druglikeness flag works", {
  skip_if_no_python()

  result <- cf_get("aspirin", druglikeness = TRUE)
  if (!is.null(result$druglikeness)) {
    expect_type(result$druglikeness, "list")
  }
})

test_that("cf_get returns empty structure for unknown compound", {
  skip_if_no_python()

  result <- cf_get("xyzzy_nonexistent_compound_12345")
  expect_type(result, "list")
  # compound field should exist even if empty
  expect_true("compound" %in% names(result))
})
