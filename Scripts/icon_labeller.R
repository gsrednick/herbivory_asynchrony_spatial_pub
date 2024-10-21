library(magick)
library(grid)

# Function to read image and convert to rasterGrob
read_and_convert_image <- function(file_path) {
  image <- magick::image_read(file_path) %>%
    magick::image_background("none") %>%
    grid::rasterGrob()
  return(image)
}

# List of file paths for the icons
icon_paths <- list(
  FH_NH = "./Figures/MS_figs/inset_icons/300ppi/FH_NH.png",
  FL_NH = "./Figures/MS_figs/inset_icons/300ppi/FL_NH.png",
  FH_NL = "./Figures/MS_figs/inset_icons/300ppi/FH_NL.png",
  FL_NL = "./Figures/MS_figs/inset_icons/300ppi/FL_NL.png"
)

# Read images and convert to rasterGrob
icons <- lapply(icon_paths, read_and_convert_image)

# List of labels associated with each icon
labels <- c("high herbivore_high nutrient" = "<img src=./Figures/MS_figs/inset_icons/300ppi/FH_NH.png",
            "low herbivore_high nutrient" = "<img src=./Figures/MS_figs/inset_icons/300ppi/FL_NH.png",
            "high herbivore_low nutrient" = "<img src=./Figures/MS_figs/inset_icons/300ppi/FH_NL.png",
            "low herbivore_low nutrient" = "<img src=./Figures/MS_figs/inset_icons/300ppi/FL_NL.png")





# Combine icons and labels into a named list
icon_list <- setNames(icons, labels)

# Define a custom labeller function
icon_labeller <- function(variable, value) {
  return(icon_list[value])
}

# Example usage of the labeller function
as_labeller(icon_labeller)




#### NEW

icon_markdown <- paste0("<img src=", list.files("Figures/MS_figs/inset_icons/300ppi/", full.names = TRUE))





theme_icon <- function(base_size = 20,
                       title_size = 20,
                       ...){
  # CUSTOM THEME:
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      # title
      plot.title = element_text(size = title_size),
      plot.title.position = "plot",

      #strip
      strip.text.x = element_markdown(size = 30),
      ...
    )
}

# other

icon_labeller <- function(variable, value) {
  # Split the facet level combination into herbivore and nutrient treatments
  treatments <- unlist(strsplit(value, "_"))
  herb_treat_sp <- treatments[1]
  nutrient_treat_sp <- treatments[2]

  # Construct the label for the icon
  icon_label <- paste0(herb_treat_sp, " herbivore_", nutrient_treat_sp, " nutrient")

  # Return the corresponding icon
  return(icon_markdown[icon_label])
}



# Modify ggplot code to use the custom labeller function and set strip width


image_to_grob <- function(image) {
  raster_grob <- grid::rasterGrob(image = image[[1]], interpolate = FALSE)
  return(raster_grob)
}

# Convert images to grobs
FH_NH_grob <- image_to_grob(FH_NH_image)
FL_NH_grob <- image_to_grob(FL_NH_image)
FH_NL_grob <- image_to_grob(FH_NL_image)
FL_NL_grob <- image_to_grob(FL_NL_image)

# Create a list of grobs and corresponding labels
facet_labels <- list(
  "FH_NH" = FH_NH_grob,
  "FL_NH" = FL_NH_grob,
  "FH_NL" = FH_NL_grob,
  "FL_NL" = FL_NL_grob
)

# Create a function to return grob based on label
custom_labeller <- function(variable, value) {
  return(facet_labels[value])
}



