# Plan and fetch in one step, approving the plan without asking.
plan_and_fetch <- function(
        ..., return_type = "matrix", coordinate_columns = NULL) {
    plan <- plan_fetch(..., confirm = FALSE)
    fetch_data(
        plan,
        return_type = return_type,
        coordinate_columns = coordinate_columns
    )
}
