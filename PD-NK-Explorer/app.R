# =============================================================
# app.R
# PD-NK Explorer
# An Interactive Platform for Exploring Natural Killer Cell
# Transcriptomic Biomarkers in Parkinson's Disease
#
# Run with: shiny::runApp()
# =============================================================

source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\global.R")

# ---------------------------------------------------------------
# UI
# ---------------------------------------------------------------
ui <- bs4DashPage(
  title = APP_TITLE,
  dark = NULL,
  help = NULL,
  header = bs4DashNavbar(
    title = bs4DashBrand(title = APP_TITLE, color = "primary"),
    status = "white"
  ),
  sidebar = bs4DashSidebar(
    status = "primary",
    bs4SidebarMenu(
      id = "sidebar_menu",
      bs4SidebarMenuItem("Home / Dashboard",        tabName = "dashboard",  icon = icon("house")),
      bs4SidebarMenuItem("Dataset Explorer",         tabName = "dataset",   icon = icon("database")),
      bs4SidebarMenuItem("Single-Cell Explorer",     tabName = "singlecell",icon = icon("circle-nodes")),
      bs4SidebarMenuItem("NK Cell Explorer",         tabName = "nkcell",    icon = icon("shield-virus")),
      bs4SidebarMenuItem("Differential Expression",  tabName = "dge",       icon = icon("chart-line")),
      bs4SidebarMenuItem("Pathway Analysis",         tabName = "pathway",   icon = icon("diagram-project")),
      bs4SidebarMenuItem("Cross-Dataset Validation", tabName = "validation",icon = icon("check-double")),
      bs4SidebarMenuItem("Biomarker Analysis",       tabName = "biomarker", icon = icon("dna")),
      bs4SidebarMenuItem("Download Results",         tabName = "download",  icon = icon("download")),
      bs4SidebarMenuItem("About",                    tabName = "about",     icon = icon("circle-info"))
    )
  ),
  body = bs4DashBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
      # hide help / dark-mode / controlbar buttons and footer; keep dropdown menus on top
      tags$style(HTML("
        [data-widget='control-sidebar'], .main-footer { display: none !important; }
        .dropdown-menu.show, .bootstrap-select .dropdown-menu { z-index: 99999 !important; }
      "))
    ),
    bs4TabItems(
      bs4TabItem(tabName = "dashboard",  mod_dashboard_ui("dashboard")),
      bs4TabItem(tabName = "dataset",    mod_dataset_ui("dataset")),
      bs4TabItem(tabName = "singlecell", mod_singlecell_ui("singlecell")),
      bs4TabItem(tabName = "nkcell",     mod_nkcell_ui("nkcell")),
      bs4TabItem(tabName = "dge",        mod_dge_ui("dge")),
      bs4TabItem(tabName = "pathway",    mod_pathway_ui("pathway")),
      bs4TabItem(tabName = "validation", mod_validation_ui("validation")),
      bs4TabItem(tabName = "biomarker",  mod_biomarker_ui("biomarker")),
      bs4TabItem(tabName = "download",   mod_download_ui("download")),
      bs4TabItem(tabName = "about",      mod_about_ui("about"))
    )
  )
)

# ---------------------------------------------------------------
# SERVER
# ---------------------------------------------------------------
server <- function(input, output, session) {
  mod_dashboard_server("dashboard")
  mod_dataset_server("dataset")
  mod_singlecell_server("singlecell")
  mod_nkcell_server("nkcell")
  mod_dge_server("dge")
  mod_pathway_server("pathway")
  mod_validation_server("validation")
  mod_biomarker_server("biomarker")
  mod_download_server("download")
  mod_about_server("about")
}

shinyApp(ui, server)