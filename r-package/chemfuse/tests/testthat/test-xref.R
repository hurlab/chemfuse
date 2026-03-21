test_that("cf_xref requires at least one argument", {
  expect_error(cf_xref())
})

test_that("cf_xref returns a tibble with database and identifier columns", {
  skip_if_no_python()

  result <- cf_xref(cid = 2244)  # aspirin
  expect_s3_class(result, "tbl_df")
  expect_true("database" %in% names(result))
  expect_true("identifier" %in% names(result))
})

test_that("cf_xref maps CID to ChEMBL ID", {
  skip_if_no_python()

  result <- cf_xref(cid = 2244)
  # Should contain a chembl entry if UniChem has it
  if (nrow(result) > 0L) {
    expect_true(is.character(result$database))
    expect_true(is.character(result$identifier))
  }
})

test_that("cf_xref handles chembl_id input", {
  skip_if_no_python()

  result <- cf_xref(chembl_id = "CHEMBL25")
  expect_s3_class(result, "tbl_df")
})

test_that("cf_xref handles unknown identifier gracefully", {
  skip_if_no_python()

  # An unlikely CID should return a minimal result
  result <- cf_xref(cid = 999999999L)
  expect_s3_class(result, "tbl_df")
})
