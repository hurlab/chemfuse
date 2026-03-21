test_that("cf_descriptors validates input", {
  expect_error(cf_descriptors(NULL))
  expect_error(cf_descriptors(character(0)))
  expect_error(cf_descriptors(123))
})

test_that("cf_descriptors returns a tibble with numeric columns", {
  skip_if_no_python()

  result <- cf_descriptors("CCO")
  expect_s3_class(result, "tbl_df")
  expect_true("smiles" %in% names(result))
  # All non-smiles/error columns should be numeric
  numeric_cols <- names(result)[!names(result) %in% c("smiles", "error")]
  for (col in numeric_cols) {
    expect_true(is.numeric(result[[col]]), info = paste("Column", col, "should be numeric"))
  }
})

test_that("cf_descriptors handles multiple SMILES", {
  skip_if_no_python()

  result <- cf_descriptors(c("CCO", "CCC", "CC(=O)O"))
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3L)
})

test_that("cf_admet validates input", {
  expect_error(cf_admet(NULL))
  expect_error(cf_admet(character(0)))
})

test_that("cf_admet returns a tibble", {
  skip_if_no_python()

  result <- cf_admet("CCO")
  expect_s3_class(result, "tbl_df")
})

test_that("cf_druglikeness validates input", {
  expect_error(cf_druglikeness(NULL))
  expect_error(cf_druglikeness(""))
  expect_error(cf_druglikeness(c("CCO", "CCC")))  # only one SMILES allowed
})

test_that("cf_druglikeness returns a list", {
  skip_if_no_python()

  result <- cf_druglikeness("CC(=O)Oc1ccccc1C(=O)O")
  expect_type(result, "list")
})
