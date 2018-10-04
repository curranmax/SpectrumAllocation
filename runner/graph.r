
library(ggplot2)
library(reshape)

theme_Publication <- function(use_major_gridlines = TRUE, use_minor_gridlines = FALSE, base_size=14, base_family="helvetica") {
      library(grid)
      library(ggthemes)

      major_elem <- element_blank()
      if (use_major_gridlines) {
      	major_elem <- element_line(colour="#f0f0f0")
      }

      minor_elem <- element_blank()
      if (use_minor_gridlines) {
      	minor_elem <- element_line(colour="#f0f0f0")
      }

      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(2)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(size = rel(1.5)), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = major_elem,
               panel.grid.minor = minor_elem,
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.spacing = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               legend.margin = margin(t = 0, unit = "cm"),
               plot.margin = unit(c(5,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

# Generic line plot function
make_line_plot <- function(data, labels, legend_position = "top", legend_direction = "horizontal", width = 6, height = 8, base_font_size = 14) {
	# Takes the multi-column data into the correct format
	data <- melt(data, id = c("x"))

	# Plots the data as lines
	plt <- ggplot(data = data)
	plt <- plt + geom_line(aes(x = data$x, y = value, color = variable), alpha = 0.75, size = 2.0)

	# Applies the theme and color scale
	plt <- plt + theme_Publication(base_size = base_font_size) + scale_colour_Publication(labels = labels)

	# Adjusts the legend
	plt <- plt + theme(legend.position = legend_position, legend.direction = legend_direction, legend.key.size = unit(1, "cm"),
						legend.title = element_text(face = "italic", size = rel(2)),
						legend.text = element_text(size = rel(2)),
						legend.background = element_rect(fill = alpha('white', 0.0)))

	# Overrides the sizes of the lines in the legend and the opacity
	plt <- plt + guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1.0)))	
}

# Makes plots for paper
make_time_npu_s <- function(data, f_name = "time_npu_s") {
	plt <- make_line_plot(data, c("2 SMs", "SM & KS", "SM & SU"))

	plt <- plt + labs(x = "Number of Selected PUs", y = "Time per SU Request (s)", colour = "")	
	ggsave(paste(c("graphs/oakland/", f_name, ".png"), collapse = ""), plot = plt, width = 8, height = 6)
}

make_time_nss_s <- function(data, f_name = "time_nss_s") {
	plt <- make_line_plot(data, c("2 SMs", "SM & KS", "SM & SU"))

	plt <- plt + labs(x = "Number of Selected SSs", y = "Time per SU Request (s)", colour = "")	
	ggsave(paste(c("graphs/oakland/", f_name, ".png"), collapse = ""), plot = plt, width = 8, height = 6)
}

make_time_bits <- function(data, f_name = "time_bits") {
	plt <- make_line_plot(data, c("2 SMs", "SM & KS", "SM & SU"))

	plt <- plt + labs(x = "Number of Bits per Value", y = "Time per SU Request (s)", colour = "")	
	ggsave(paste(c("graphs/oakland/", f_name, ".png"), collapse = ""), plot = plt, width = 8, height = 6)
}

make_secure_read <- function(data, f_name = "secure_read") {
	plt <- make_line_plot(data, c("2 SMs + Naive", "SM & KS + Naive", "2 SMs + Proposed", "SM & KS + Proposed"))

	plt <- plt + labs(x = "Size of Table (MB)", y = "Time per Secure Read(s)", colour = "")	
	ggsave(paste(c("graphs/oakland/", f_name, ".png"), collapse = ""), plot = plt, width = 8, height = 6)
}

make_secure_write <- function(data, f_name = "secure_write") {
	plt <- make_line_plot(data, c("2 SMs + Naive", "SM & KS + Naive", "2 SMs + Proposed", "SM & KS + Proposed"))

	plt <- plt + labs(x = "Selected PUs", y = "Time per Secure Write(s)", colour = "")	
	ggsave(paste(c("graphs/oakland/", f_name, ".png"), collapse = ""), plot = plt, width = 8, height = 6)
}

make_ple_cdf <- function(data, f_name = "ple_cdf") {
	plt <- make_line_plot(data, c("1 SS", "10 SS", "25 SS", "50 SS"), legend_position = "right", legend_direction = "vertical")

	plt <- plt + labs(x = "Error of Estimated Path Loss", y = "CDF", colour = "# Selected:")	
	ggsave(paste(c("graphs/oakland/", f_name, ".png"), collapse = ""), plot = plt, width = 8, height = 6)
}

# Reads in data and creates plots
time_npu_s_data = data.frame("x" = c(1,2,3,4), "a" = c(0.1, 0.2, 0.3, 0.4), "b" = c(0.2, 0.4, 0.6, 0.8), "c" = c(0.5, 1.0, 1.5, 2.0))
make_time_npu_s(time_npu_s_data)

# time_nss_s_data = data.frame("x" = c(1,2,3,4), "a" = c(0.1, 0.2, 0.3, 0.4), "b" = c(0.2, 0.4, 0.6, 0.8), "c" = c(0.5, 1.0, 1.5, 2.0))
# make_time_nss_s(time_nss_s_data)

# time_bits_data = data.frame("x" = c(1,2,3,4), "a" = c(0.1, 0.2, 0.3, 0.4), "b" = c(0.2, 0.4, 0.6, 0.8), "c" = c(0.5, 1.0, 1.5, 2.0))
# make_time_bits(time_bits_data)

# secure_read_data = data.frame("x" = c(1,2,3,4), "c" = c(0.5, 1.0, 1.5, 2.0), "d" = c(1.0, 2.0, 3.0, 4.0), "a" = c(0.1, 0.2, 0.3, 0.4), "b" = c(0.2, 0.4, 0.6, 0.8))
# make_secure_read(secure_read_data)

# secure_write_data = data.frame("x" = c(1,2,3,4), "c" = c(0.5, 1.0, 1.5, 2.0), "d" = c(1.0, 2.0, 3.0, 4.0), "a" = c(0.1, 0.2, 0.3, 0.4), "b" = c(0.2, 0.4, 0.6, 0.8))
# make_secure_write(secure_write_data)

# ple_cdf_data = data.frame("x" = c(1,2,3,4), "c" = c(0.5, 1.0, 1.5, 2.0), "d" = c(1.0, 2.0, 3.0, 4.0), "a" = c(0.1, 0.2, 0.3, 0.4), "b" = c(0.2, 0.4, 0.6, 0.8))
# make_ple_cdf(ple_cdf_data)
