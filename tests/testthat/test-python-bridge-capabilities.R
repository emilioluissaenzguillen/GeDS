test_that("Python bridge capabilities identify the tested R behavior", {
  capabilities <- getFromNamespace(".GeDS_python_bridge_capabilities", "GeDS")
  expect_setequal(capabilities, c("offset_glm_v1", "terms_matrix_v1"))
})
