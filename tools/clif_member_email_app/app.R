library(shiny)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(DT)
library(bslib)

default_workbook_path <- "/Users/saborpete/Downloads/Untitled spreadsheet.xlsx"
email_pattern <- "[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}"

extract_first_email <- function(x) {
  value <- stringr::str_extract(x %||% "", regex(email_pattern, ignore_case = TRUE))
  ifelse(is.na(value) | value == "", NA_character_, stringr::str_to_lower(value))
}

clean_status <- function(x) {
  status <- trimws(x %||% "")
  ifelse(status == "", "Unknown", status)
}

derive_site <- function(affiliation, email, gmail) {
  source_text <- paste(affiliation %||% "", email %||% "", gmail %||% "", sep = " | ")

  dplyr::case_when(
    str_detect(source_text, regex("University of Chicago|UChicago|MacLean|Harris School|Pritzker|uchicago", ignore_case = TRUE)) ~ "University of Chicago",
    str_detect(source_text, regex("Northwestern|northwestern\\.edu", ignore_case = TRUE)) ~ "Northwestern",
    str_detect(source_text, regex("Johns Hopkins|jhmi\\.edu|jh\\.edu", ignore_case = TRUE)) ~ "Johns Hopkins",
    str_detect(source_text, regex("Brigham and Women|bwh\\.harvard\\.edu", ignore_case = TRUE)) ~ "Brigham and Women's",
    str_detect(source_text, regex("Intermountain", ignore_case = TRUE)) ~ "Intermountain",
    str_detect(source_text, regex("University of Michigan|med\\.umich\\.edu", ignore_case = TRUE)) ~ "University of Michigan",
    str_detect(source_text, regex("Rush|rush\\.edu", ignore_case = TRUE)) ~ "Rush",
    str_detect(source_text, regex("Emory|emory\\.edu", ignore_case = TRUE)) ~ "Emory",
    str_detect(source_text, regex("University of Pennsylvania|UPenn|Penn|upenn\\.edu|pennmedicine\\.upenn\\.edu", ignore_case = TRUE)) ~ "University of Pennsylvania",
    str_detect(source_text, regex("University of Minnesota|UMN|umn\\.edu", ignore_case = TRUE)) ~ "University of Minnesota",
    str_detect(source_text, regex("University of Virginia|UVA|uvahealth\\.org", ignore_case = TRUE)) ~ "University of Virginia",
    str_detect(source_text, regex("Weill Cornell|Cornell|med\\.cornell\\.edu", ignore_case = TRUE)) ~ "Weill Cornell",
    str_detect(source_text, regex("Oregon Health\\s*&\\s*Science University|OHSU|ohsu\\.edu", ignore_case = TRUE)) ~ "OHSU",
    str_detect(source_text, regex("Rutgers", ignore_case = TRUE)) ~ "Rutgers",
    str_detect(source_text, regex("Tufts|tuftsmedicine\\.org", ignore_case = TRUE)) ~ "Tufts",
    str_detect(source_text, regex("Sunnybrook|sunnybrook\\.ca", ignore_case = TRUE)) ~ "Sunnybrook",
    str_detect(source_text, regex("UCSF|ucsf\\.edu", ignore_case = TRUE)) ~ "UCSF",
    str_detect(source_text, regex("University of Colorado|cuanschutz\\.edu", ignore_case = TRUE)) ~ "University of Colorado",
    str_detect(source_text, regex("Yale|yale\\.edu", ignore_case = TRUE)) ~ "Yale",
    TRUE ~ "Other"
  )
}

read_member_workbook <- function(path, prefer_gmail = FALSE) {
  if (!nzchar(path) || !file.exists(path)) {
    stop(paste("Workbook not found:", path), call. = FALSE)
  }

  raw <- readxl::read_excel(path, sheet = 1) %>%
    mutate(
      Name = trimws(Name %||% ""),
      Affiliation = coalesce(Affiliation, ""),
      work_email = extract_first_email(Email),
      gmail_email = extract_first_email(Gmail),
      selected_email = ifelse(prefer_gmail & !is.na(gmail_email), gmail_email, work_email),
      selected_email = ifelse(is.na(selected_email) & !is.na(work_email), work_email, selected_email),
      selected_email = ifelse(is.na(selected_email) & !is.na(gmail_email), gmail_email, selected_email),
      member_status_clean = clean_status(`Member status`),
      site = derive_site(Affiliation, Email, Gmail)
    ) %>%
    filter(Name != "") %>%
    distinct(Name, selected_email, .keep_all = TRUE) %>%
    arrange(site, Name) %>%
    mutate(
      member_id = row_number(),
      member_label = ifelse(
        is.na(selected_email),
        paste0(Name, " (", site, ")"),
        paste0(Name, " (", site, " | ", selected_email, ")")
      )
    )

  raw
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

site_input_id <- function(site) {
  paste0("site_members__", gsub("[^A-Za-z0-9]+", "_", site))
}

site_select_all_id <- function(site) {
  paste0(site_input_id(site), "__select_all")
}

site_clear_id <- function(site) {
  paste0(site_input_id(site), "__clear")
}

ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  tags$head(
    tags$style(HTML("
      .app-header {
        background: linear-gradient(135deg, #0f3d3e 0%, #1f6f78 100%);
        color: white;
        padding: 20px 24px;
        border-radius: 18px;
        margin-bottom: 20px;
      }
      .metric-card {
        background: #f8fbfc;
        border: 1px solid #d9e7ea;
        border-radius: 14px;
        padding: 14px 16px;
        margin-bottom: 14px;
      }
      .metric-value {
        font-size: 1.5rem;
        font-weight: 700;
        color: #0f3d3e;
      }
      .small-note {
        color: #4f6b72;
        font-size: 0.95rem;
      }
      .panel-box {
        background: white;
        border: 1px solid #dde7ea;
        border-radius: 16px;
        padding: 16px;
        margin-bottom: 18px;
        box-shadow: 0 10px 25px rgba(14, 61, 62, 0.05);
      }
      textarea.form-control {
        font-family: Menlo, Consolas, monospace;
      }
      .site-roster {
        border: 1px solid #d9e7ea;
        border-radius: 14px;
        padding: 12px 14px;
        margin-bottom: 12px;
        background: #fbfdfd;
      }
      .site-roster h5 {
        margin-top: 0;
        margin-bottom: 8px;
        color: #0f3d3e;
      }
      .site-roster-header {
        display: flex;
        justify-content: space-between;
        align-items: flex-start;
        gap: 12px;
        margin-bottom: 8px;
      }
      .site-roster-actions {
        display: flex;
        gap: 8px;
        flex-wrap: wrap;
      }
      .site-roster-actions .btn {
        padding: 4px 10px;
      }
      .shiny-options-group {
        max-height: 220px;
        overflow-y: auto;
      }
    ")),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('copy-text', function(message) {
        function updateButton() {
          const button = document.getElementById(message.button_id);
          if (!button) return;
          const original = button.innerText;
          button.innerText = 'Copied';
          setTimeout(function() { button.innerText = original; }, 1200);
        }

        function fallbackCopy(text) {
          const helper = document.createElement('textarea');
          helper.value = text;
          helper.setAttribute('readonly', '');
          helper.style.position = 'fixed';
          helper.style.opacity = '0';
          document.body.appendChild(helper);
          helper.focus();
          helper.select();
          try {
            document.execCommand('copy');
          } catch (e) {}
          document.body.removeChild(helper);
          updateButton();
        }

        if (navigator.clipboard && window.isSecureContext) {
          navigator.clipboard.writeText(message.text)
            .then(updateButton)
            .catch(function() { fallbackCopy(message.text); });
        } else {
          fallbackCopy(message.text);
        }
      });
    "))
  ),
  titlePanel(NULL),
  div(
    class = "app-header",
    h2("CLIF Consortium Email Builder"),
    p(
      class = "mb-0",
      "Filter by site, narrow to specific people, and copy a ready-to-paste recipient list for abstracts and manuscripts."
    )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      div(
        class = "panel-box",
        textInput("workbook_path", "Workbook path", value = if (file.exists(default_workbook_path)) default_workbook_path else ""),
        fileInput("workbook_upload", "Or upload a workbook", accept = c(".xlsx", ".xls")),
        checkboxInput("prefer_gmail", "Use Gmail when available", value = FALSE),
        actionButton("reload_data", "Reload workbook", class = "btn-primary"),
        br(),
        br(),
        uiOutput("data_status")
      ),
      div(
        class = "panel-box",
        checkboxGroupInput("site_filter", "Sites", choices = NULL),
        selectizeInput("status_filter", "Member status", choices = NULL, multiple = TRUE),
        textInput("search_text", "Search people", placeholder = "Name, email, affiliation, site"),
        fluidRow(
          column(6, actionButton("select_all_filtered", "Select all visible")),
          column(6, actionButton("clear_selection", "Clear selection"))
        )
      ),
      div(
        class = "panel-box",
        h4("People by Selected Site"),
        uiOutput("site_member_checkboxes")
      )
    ),
    mainPanel(
      width = 8,
      fluidRow(
        column(4, div(class = "metric-card", div("Filtered members"), div(class = "metric-value", textOutput("filtered_n", inline = TRUE)))),
        column(4, div(class = "metric-card", div("Selected people"), div(class = "metric-value", textOutput("selected_n", inline = TRUE)))),
        column(4, div(class = "metric-card", div("Emails ready"), div(class = "metric-value", textOutput("email_n", inline = TRUE))))
      ),
      div(
        class = "panel-box",
        h4("Filtered Member List"),
        DTOutput("member_table")
      ),
      div(
        class = "panel-box",
        h4("Selected Recipients"),
        DTOutput("selected_table")
      ),
      div(
        class = "panel-box",
        h4("Recipient String"),
        textAreaInput("recipient_text", NULL, value = "", width = "100%", height = "130px", resize = "vertical"),
        fluidRow(
          column(4, actionButton("copy_semicolon", "Copy with semicolons", class = "btn-success")),
          column(4, actionButton("copy_comma", "Copy with commas")),
          column(4, downloadButton("download_selected", "Download CSV"))
        ),
        br(),
        uiOutput("mailto_ui"),
        div(class = "small-note", "Semicolons usually work best for Outlook. Commas are often better for Gmail and Apple Mail.")
      )
    )
  )
)

server <- function(input, output, session) {
  registered_site_observers <- reactiveVal(character(0))

  workbook_path <- reactive({
    uploaded <- input$workbook_upload
    if (!is.null(uploaded) && !is.null(uploaded$datapath) && nzchar(uploaded$datapath)) {
      uploaded$datapath
    } else {
      trimws(input$workbook_path %||% "")
    }
  })

  members <- eventReactive(
    list(input$reload_data, input$prefer_gmail, workbook_path()),
    {
      read_member_workbook(workbook_path(), prefer_gmail = isTRUE(input$prefer_gmail))
    },
    ignoreNULL = FALSE
  )

  observeEvent(members(), {
    df <- members()

    updateCheckboxGroupInput(
      session,
      "site_filter",
      choices = sort(unique(df$site)),
      selected = sort(unique(df$site))
    )

    updateSelectizeInput(
      session,
      "status_filter",
      choices = sort(unique(df$member_status_clean)),
      selected = sort(unique(df$member_status_clean)),
      server = TRUE
    )
  }, ignoreInit = FALSE)

  filtered_members <- reactive({
    req(members())
    df <- members()

    if (!is.null(input$site_filter) && length(input$site_filter) > 0) {
      df <- df %>% filter(site %in% input$site_filter)
    }

    if (!is.null(input$status_filter) && length(input$status_filter) > 0) {
      df <- df %>% filter(member_status_clean %in% input$status_filter)
    }

    query <- trimws(input$search_text %||% "")
    if (nzchar(query)) {
      df <- df %>%
        filter(
          str_detect(
            paste(Name, selected_email, site, member_status_clean, Affiliation, sep = " | "),
            regex(query, ignore_case = TRUE)
          )
        )
    }

    df
  })

  observe({
    req(members())
    current_sites <- input$site_filter %||% character(0)
    for (site in current_sites) {
      site_df <- filtered_members() %>%
        filter(site == !!site) %>%
        arrange(Name)

      input_id <- site_input_id(site)
      existing <- input[[input_id]] %||% character(0)
      valid_selected <- existing[existing %in% as.character(site_df$member_id)]

      updateCheckboxGroupInput(
        session,
        inputId = input_id,
        choices = setNames(as.character(site_df$member_id), site_df$member_label),
        selected = valid_selected
      )
    }
  })

  observeEvent(input$select_all_filtered, {
    req(members())
    for (site in input$site_filter %||% character(0)) {
      site_df <- filtered_members() %>%
        filter(site == !!site)
      updateCheckboxGroupInput(
        session,
        inputId = site_input_id(site),
        selected = as.character(site_df$member_id)
      )
    }
  })

  observeEvent(input$clear_selection, {
    req(members())
    for (site in unique(members()$site)) {
      updateCheckboxGroupInput(session, inputId = site_input_id(site), selected = character(0))
    }
  })

  observeEvent(members(), {
    known_sites <- registered_site_observers()
    new_sites <- setdiff(unique(members()$site), known_sites)

    lapply(new_sites, function(site) {
      local({
        current_site <- site
        observeEvent(input[[site_select_all_id(current_site)]], {
          site_df <- filtered_members() %>%
            filter(site == current_site)

          updateCheckboxGroupInput(
            session,
            inputId = site_input_id(current_site),
            selected = as.character(site_df$member_id)
          )
        }, ignoreInit = TRUE)

        observeEvent(input[[site_clear_id(current_site)]], {
          updateCheckboxGroupInput(
            session,
            inputId = site_input_id(current_site),
            selected = character(0)
          )
        }, ignoreInit = TRUE)
      })
    })

    registered_site_observers(sort(unique(c(known_sites, new_sites))))
  }, ignoreInit = TRUE)

  selected_members <- reactive({
    req(members())
    all_ids <- unlist(
      lapply(unique(members()$site), function(site) input[[site_input_id(site)]] %||% character(0)),
      use.names = FALSE
    )
    ids <- suppressWarnings(as.integer(all_ids))
    members() %>%
      filter(member_id %in% ids) %>%
      arrange(site, Name)
  })

  output$site_member_checkboxes <- renderUI({
    req(members())
    current_sites <- input$site_filter %||% character(0)

    if (length(current_sites) == 0) {
      return(div(class = "small-note", "Tick one or more sites above to see that site's roster."))
    }

    roster_panels <- lapply(current_sites, function(site) {
      site_df <- filtered_members() %>%
        filter(site == !!site) %>%
        arrange(Name)

      if (nrow(site_df) == 0) {
        return(
          div(
            class = "site-roster",
            div(
              class = "site-roster-header",
              div(
                h5(site),
                div(class = "small-note", "0 people shown")
              ),
              div(
                class = "site-roster-actions",
                actionButton(site_select_all_id(site), "Select all", class = "btn btn-outline-primary btn-sm"),
                actionButton(site_clear_id(site), "Clear", class = "btn btn-outline-secondary btn-sm")
              )
            ),
            div(class = "small-note", "No matching people for the current status/search filters.")
          )
        )
      }

      div(
        class = "site-roster",
        div(
          class = "site-roster-header",
          div(
            h5(site),
            div(class = "small-note", paste(nrow(site_df), "people shown"))
          ),
          div(
            class = "site-roster-actions",
            actionButton(site_select_all_id(site), "Select all", class = "btn btn-outline-primary btn-sm"),
            actionButton(site_clear_id(site), "Clear", class = "btn btn-outline-secondary btn-sm")
          )
        ),
        checkboxGroupInput(
          inputId = site_input_id(site),
          label = NULL,
          choices = setNames(as.character(site_df$member_id), site_df$member_label),
          selected = input[[site_input_id(site)]] %||% character(0)
        )
      )
    })

    tagList(roster_panels)
  })

  recipient_string <- reactive({
    emails <- selected_members() %>%
      filter(!is.na(selected_email), selected_email != "") %>%
      pull(selected_email) %>%
      unique()

    paste(emails, collapse = "; ")
  })

  output$data_status <- renderUI({
    path <- workbook_path()
    if (!nzchar(path)) {
      return(div(class = "small-note", "Enter a workbook path or upload the spreadsheet."))
    }

    div(
      class = "small-note",
      strong("Current source:"),
      tags$br(),
      path
    )
  })

  output$filtered_n <- renderText(nrow(filtered_members()))
  output$selected_n <- renderText(nrow(selected_members()))
  output$email_n <- renderText(sum(!is.na(selected_members()$selected_email) & selected_members()$selected_email != ""))

  output$member_table <- renderDT({
    df <- filtered_members() %>%
      transmute(
        Name,
        Site = site,
        Status = member_status_clean,
        Email = selected_email,
        `Work email` = work_email,
        Gmail = gmail_email,
        Affiliation
      )

    datatable(
      df,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 12, scrollX = TRUE)
    )
  })

  output$selected_table <- renderDT({
    df <- selected_members() %>%
      transmute(
        Name,
        Site = site,
        Status = member_status_clean,
        Email = selected_email,
        `Work email` = work_email,
        Gmail = gmail_email
      )

    datatable(
      df,
      rownames = FALSE,
      options = list(dom = "tip", pageLength = 10, scrollX = TRUE)
    )
  })

  observe({
    updateTextAreaInput(session, "recipient_text", value = recipient_string())
  })

  observeEvent(input$copy_semicolon, {
    session$sendCustomMessage(
      "copy-text",
      list(
        text = recipient_string(),
        button_id = "copy_semicolon"
      )
    )
  })

  observeEvent(input$copy_comma, {
    text <- selected_members() %>%
      filter(!is.na(selected_email), selected_email != "") %>%
      pull(selected_email) %>%
      unique() %>%
      paste(collapse = ", ")

    session$sendCustomMessage(
      "copy-text",
      list(
        text = text,
        button_id = "copy_comma"
      )
    )
  })

  output$download_selected <- downloadHandler(
    filename = function() {
      paste0("clif_email_selection_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      selected_members() %>%
        transmute(
          Name,
          Site = site,
          Status = member_status_clean,
          Email = selected_email,
          `Work email` = work_email,
          Gmail = gmail_email,
          Affiliation
        ) %>%
        write.csv(file, row.names = FALSE, na = "")
    }
  )

  output$mailto_ui <- renderUI({
    emails <- selected_members() %>%
      filter(!is.na(selected_email), selected_email != "") %>%
      pull(selected_email) %>%
      unique()

    if (length(emails) == 0) {
      return(div(class = "small-note", "Select at least one member with an email address to generate a mailto link."))
    }

    tags$a(
      href = paste0("mailto:", utils::URLencode(paste(emails, collapse = ","), reserved = TRUE)),
      target = "_blank",
      class = "btn btn-outline-secondary",
      "Open draft email"
    )
  })
}

shinyApp(ui, server)
