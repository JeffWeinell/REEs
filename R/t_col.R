#' Make Transparent Color
#' 
#' Makes a transparent version of a named color. Function is from "Transparent colors. Mark Gardener 2015. www.dataanalytics.org.uk"
#' 
#' @param color color name
#' @param percent % transparency (default is 50)
#' @param name Optional name for the color (default is NULL)
#' @return Returns transparent version of a named color.
#' @export
t_col <- function(color, percent = 50, name = NULL) {
	rgb.val <- col2rgb(color)
	t.col   <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],max = 255,alpha = (100 - percent) * 255 / 100,names = name)
	invisible(t.col)
}
