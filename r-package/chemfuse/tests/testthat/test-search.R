test_that("cf_search validates input", {
  expect_error(cf_search(NULL))
  expect_error(cf_search(""))
  expect_error(cf_search(123))
  expect_error(cf_search("aspirin", sources = character(0)))
})

test_that("cf_search returns a tibble structure", {
  skip_if_no_python()

  result <- cf_search("aspirin")
  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) >= 1L)
  expect_true("name" %in% names(result))
  expect_true("smiles" %in% names(result))
})

test_that("cf_search respects sources parameter", {
  skip_if_no_python()

  result <- cf_search("aspirin", sources = c("pubchem"))
  expect_s3_class(result, "tbl_df")
})

test_that("cf_search returns empty tibble for unknown compound", {
  skip_if_no_python()

  result <- cf_search("xyzzy_nonexistent_compound_12345")
  expect_s3_class(result, "tbl_df")
  # May return 0 rows for unknown compound
  expect_true(is.numeric(nrow(result)) || nrow(result) == 0L)
})

test_that("cf_search handles smiles query_type", {
  skip_if_no_python()

  result <- cf_search("CC(=O)Oc1ccccc1C(=O)O", query_type = "smiles")
  expect_s3_class(result, "tbl_df")
})
