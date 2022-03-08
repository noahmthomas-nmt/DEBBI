#' hello
#'
#' @param name A name
#'
#' @return A message
#' @export
#'
#' @examples hello("Alice")

hello <- function(name){
  message <- paste0("Hello ", name, "!")

  return(message)
}
