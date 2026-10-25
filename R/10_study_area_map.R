# 10_study_area_map.R
# Study area map (Fig 3.1): sampling locations in Dhaka coloured by project,
# MRT Line 6 drawn as a line, with north arrow and scale bar.
#
# Coordinates are APPROXIMATE (WGS84, decimal degrees) and for illustration only.
# MRT Line 6 stations are placed at the published station sites; points marked
# `certain = FALSE` are rough locations and should be checked against the
# borehole logs before the map is finalised.
# Project assignment of the non-MRT sites is an assumption - check it.
#
# Basemap tiles (OpenStreetMap) are downloaded on the first
# run and cached in data/map_tiles/, so later runs work offline.

source("R/00_packages.R")
set.seed(123)

theme_set(theme_bw(base_size = 10, base_family = "Times New Roman"))

# ---- Locations -------------------------------------------------------------
sites <- tribble(
  ~site,               ~lat,    ~lon,    ~project,              ~certain,
  "Kawlar",            23.8570, 90.4180, "Aga Khan Academy",    FALSE,
  "Ashkona",           23.8625, 90.4235, "Aga Khan Academy",    FALSE,
  "Purbachal",         23.8400, 90.5000, "100 ft Khal project", FALSE,
  "300 ft road",       23.8320, 90.4650, "100 ft Khal project", FALSE,
  "Uttara",            23.8693, 90.3673, "MRT Line 6",          TRUE,
  "Pallabi",           23.8265, 90.3643, "MRT Line 6",          TRUE,
  "Mirpur 11",         23.8193, 90.3654, "MRT Line 6",          TRUE,
  "Mirpur 10",         23.8071, 90.3686, "MRT Line 6",          TRUE,
  "Kazipara",          23.7997, 90.3721, "MRT Line 6",          TRUE,
  "Shewrapara",        23.7905, 90.3757, "MRT Line 6",          TRUE,
  "Agargaon",          23.7782, 90.3801, "MRT Line 6",          TRUE,
  "Bijoy Sarani",      23.7665, 90.3833, "MRT Line 6",          TRUE,
  "Farmgate",          23.7579, 90.3897, "MRT Line 6",          TRUE,
  "Kawran Bazar",      23.7514, 90.3930, "MRT Line 6",          TRUE,
  "Shahbag",           23.7389, 90.3958, "MRT Line 6",          TRUE,
  "Dhaka University",  23.7329, 90.3976, "MRT Line 6",          TRUE,
  "Secretariat",       23.7298, 90.4074, "MRT Line 6",          TRUE,
  "Motijheel",         23.7280, 90.4176, "MRT Line 6",          TRUE
) %>%
  mutate(project = factor(project, levels = c("Aga Khan Academy",
                                              "100 ft Khal project",
                                              "MRT Line 6")))

# MRT Line 6 alignment: Uttara North to Motijheel, including the two Uttara
# stations that are not sampling sites (Uttara Center, Uttara South)
mrt_path <- tribble(
  ~lat,    ~lon,
  23.8693, 90.3673,   # Uttara North
  23.8597, 90.3655,   # Uttara Center
  23.8457, 90.3632,   # Uttara South
) %>%
  bind_rows(sites %>% filter(project == "MRT Line 6", site != "Uttara") %>%
              select(lat, lon))

sites_sf <- sf::st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)
mrt_sf <- sf::st_sfc(sf::st_linestring(as.matrix(mrt_path[, c("lon", "lat")])),
                     crs = 4326) %>%
  sf::st_sf(name = "MRT Line 6 alignment", geometry = .)

# ---- Map -------------------------------------------------------------------
tile_cache <- "data/map_tiles"
dir.create(tile_cache, recursive = TRUE, showWarnings = FALSE)
# Standard OpenStreetMap tiles (attribution required)
basemap_url <- "osm"

project_cols <- c("Aga Khan Academy" = "#D55E00",
                  "100 ft Khal project" = "#009E73",
                  "MRT Line 6" = "#0072B2")

# Labels: sites with uncertain coordinates get an asterisk
sites_lab <- sites %>% mutate(label = ifelse(certain, site, paste0(site, "*")))

fig3_1 <- ggplot() +
  ggspatial::annotation_map_tile(type = basemap_url, zoomin = 0,
                                 cachedir = tile_cache, progress = "none") +
  geom_sf(data = mrt_sf, aes(linetype = name), colour = "#0072B2", linewidth = 1.1) +
  geom_sf(data = sites_sf, aes(fill = project), shape = 21, colour = "black",
          size = 2.6, stroke = 0.4) +
  ggrepel::geom_text_repel(data = sites_lab, aes(x = lon, y = lat, label = label),
                           family = "Times New Roman", size = 2.8, seed = 123,
                           min.segment.length = 0.1, box.padding = 0.35,
                           bg.color = "white", bg.r = 0.12) +
  scale_fill_manual(values = project_cols, name = "Project") +
  scale_linetype_manual(values = "solid", name = NULL) +
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
                                    height = unit(1.1, "cm"), width = unit(0.9, "cm"),
                                    style = ggspatial::north_arrow_fancy_orienteering(
                                      text_family = "Times New Roman")) +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.25,
                              text_family = "Times New Roman") +
  coord_sf(xlim = c(90.345, 90.520), ylim = c(23.715, 23.885), crs = 4326,
           default_crs = 4326, expand = FALSE) +
  labs(x = "Longitude", y = "Latitude",
       caption = paste("* Approximate location, to be confirmed.",
                       "Basemap: © OpenStreetMap contributors.")) +
  guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.02),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = "grey60"),
        legend.spacing.y = unit(2, "pt"))

ggsave("outputs/figures/Fig3_1_study_area.png", fig3_1, width = 7, height = 7.2,
       units = "in", dpi = 300, device = ragg::agg_png)

# ---- Plain-English summary -------------------------------------------------
cat("\n--- Fig 3.1: study area ---\n")
print(as.data.frame(count(sites, project, name = "n_sites")), row.names = FALSE)
cat("\nSummary:\n")
cat("- ", nrow(sites), " locations mapped: ",
    sum(sites$project == "MRT Line 6"), " MRT Line 6 stations from Uttara to Motijheel, plus ",
    sum(sites$project != "MRT Line 6"), " sites in the east and north-east of Dhaka.\n", sep = "")
cat("- Approximate coordinates (please check):",
    paste(sites$site[!sites$certain], collapse = ", "), "\n")
cat("- Project assignment of the non-MRT sites is assumed; edit the `sites` table if it is wrong.\n")
cat("Saved outputs/figures/Fig3_1_study_area.png\n")
