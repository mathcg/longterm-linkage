test_that("paper simulation generators return documented data structures", {
  time_dependent <- simulate_linkage_data(n = 100, setting = 1, seed = 11)
  expect_s3_class(time_dependent$df, "data.frame")
  expect_s3_class(time_dependent$tdata, "data.frame")
  expect_equal(nrow(time_dependent$df), 100)

  time_independent <- simulate_linkage_data(n = 100, setting = 6, seed = 11)
  expect_s3_class(time_independent$df, "data.frame")
  expect_equal(nrow(time_independent$df), 100)
})

test_that("synthetic PCPT data are bundled and can be prepared", {
  pcpt <- read_pcpt_synthetic()
  expect_equal(nrow(pcpt), 16518)
  expect_true(all(c("id", "txg", "linkage", "Y", "Delta") %in% names(pcpt)))

  prepared <- prepare_pcpt_linkage_data(pcpt[1:200, ])
  expect_s3_class(prepared$df, "data.frame")
  expect_s3_class(prepared$tdata, "data.frame")
  expect_true(all(c("T1", "T2", "tstart", "tstop", "death") %in% names(prepared$tdata)))
})

test_that("input validation rejects unsupported settings", {
  expect_error(simulate_linkage_data(100, setting = 8), "one integer")
  expect_error(
    run_linkage_simulation(100, setting = 6, time_dependent = TRUE),
    "time-independent"
  )
  expect_error(
    run_linkage_simulation(100, method = "aiplw", setting = 1),
    "settings 6 and 7"
  )
})
