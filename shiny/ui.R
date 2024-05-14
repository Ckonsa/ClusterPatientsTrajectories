library(shinyWidgets)
library(shinycssloaders)

background_color <- "#e9f4fb"
spinner_color <- "#53a9e0"

ui <- fluidPage(
  titlePanel("Ravimitrajektooride klasterdamine"),
  setBackgroundColor(color = background_color),

  sidebarLayout(
    sidebarPanel(
      fileInput("data_file", "Vali andmefail", multiple = FALSE,
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      tags$hr(style="border-color: #8c8c8c;"),
      sliderInput("records_range", "Kirjete vahemik", min = 1, max = 100, value = c(45, 55)),
      sliderInput("trajectory_length", "Trajektoori pikkus", min = 1, max = 50, value = 25),
      sliderInput("clusters_amount", "Klastrite arv", min = 1, max = 10, value = 4),
      sliderInput("popular_drugs_amount", "Vaadeldavate sagedaste ravimite arv", min = 1, max = 20, value = 5),
      materialSwitch("use_dtw_algorithm", label = h5("Kasuta DTW algoritmi"), value = FALSE, status = "primary"),
      fluidRow(
        column(h5("ATC koodide maksumused tasemete kaupa"), width = 12),
        column(numericInput("level_one", label = "1", min = 1, max = 20, value = 14, step = 0.5), width = 4),
        column(numericInput("level_two", label = "2", min = 1, max = 20, value = 11, step = 0.5), width = 4),
        column(numericInput("level_three", label = "3", min = 1, max = 20, value = 7, step = 0.5), width = 4),
        column(numericInput("level_four", label = "4", min = 1, max = 20, value = 4, step = 0.5), width = 6),
        column(numericInput("level_five", label = "5", min = 1, max = 20, value = 2, step = 0.5), width = 6)
      ),
      fluidRow(
        column(actionButton("submit", "Salvesta"), width = 6),
        column(downloadButton("download_parameters", "Parameetrid"), width = 6, align = "right")
      ),
      tags$hr(style="border-color: #8c8c8c;"),
      withSpinner(plotOutput("elbow_plot", width = "100%"), type = 7, color = spinner_color)
    ),

    mainPanel(
      tabsetPanel(
        tabPanel(
          "Trajektooride klastrid",
          h3("Ravitrajektooride klasterdamine"),
          fluidRow(column(downloadButton("download_results", "Tulemused"), width = 12, align = "right")),
          br(),
          withSpinner(plotOutput("clustered_trajectories", width = "100%"), type = 7, color = spinner_color),
          br()
        ),
        tabPanel(
          "Ravimid klastrites",
          h3("Sagedaste ravimite kogus klastrites"),
          materialSwitch("normalize", label = h5("Normaliseeri"), value = FALSE, status = "primary"),
          withSpinner(plotOutput("clusters_plots", width = "100%", height = "100vh"), type = 7, color = spinner_color), br(),
          withSpinner(tableOutput("clusters_data"), type = 7, color = spinner_color)
        ),
        tabPanel(
          "Abi",
          h3("Vahekaardid"),

          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; text-align: justify;'>",
                      "<b> Vahekaart \"Trajektooride klastrid\".</b> Kuvatakse diagramm, kuidas toimub patsientide klasterdamine nende ravimitrajektooride
                      alusel. Diagrammi värvus annab aimu kui sarnased või erinevad on klasterdatud ravimitrajektoorid. Suurem väärtus tähendab suuremat
                      erinevust kahe trajektoori vahel. Juhul kui klasterdatavaid patsiente on rohkem kui 600, siis valitakse klasterdamiseks juhuslikult
                      600 patsienti. See on vajalik, et diagrammid oleksid loetavad ja rakendus ei võtaks tulemuste kuvamiseks liiga pikalt aega. Tulemuste
                      salvestamisel nupu “Tulemused” abil klasterdatakse ravimitrajektoorid kõigi soovitud patsientidega. Tulemused salvestatakse RDS-faili.
                      Lisaks klastritele on failis kasutatud parameetrid ja vahelehtedel “Trajektooride klastrid” ning “Ravimid klastrites” olevad diagrammid. ", "</div>")),

          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; text-align: justify;'>",
                      "<b> Vahekaart \"Ravimid klastrites\".</b> Kuvatakse lisainformatsiooni loodud klastrite kohta. Iga klastri puhul luuakse diagramm,
                      millelt saab lugeda, mitu korda on patsientidele välja kirjutatud andmestiku kõige sagedasemaid ravimeid. Diagrammi väärtusi on võimalik
                      normaliseerida diagrammi kohal asuva lüliti “Normaliseeri” abil. Normaliseerimiseks kasutatakse funktsiooni log10. Normaliseerimine
                      võimaldab diagramme omavahel hõlpsamini võrrelda. Diagrammide all kuvatakse iga klastri kohta selles olevate patsientide arv, sagedaseim ravim
                      ja keskmine trajektoori pikkus.", "</div>")),

          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; margin-bottom: 5%; text-align: justify;'>",
                      "<b> Vahekaart \"Abi\".</b> Kuvatakse informatsiooni vahekaartide ja külgriba osade kohta.", "</div>")),

          h3("Külgriba"),
          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; text-align: justify;'>",
                      "<b> Andmefail.</b> Võimaldab lisada andmefaili, milles olevaid patsiente soovitakse klasterdada ravimitrajektooride alusel.
                      Vaikeväärtusena kasutatakse näidisfail, milles pole tegelike patsientide andmeid.", "</div>")),
          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; text-align: justify;'>",
                      "<b> Kirjete vahemik.</b> Võimaldab filtreerida klasterdatavaid patsiente. Kasutajal on võimalik valida minimaalne ja
                      maksimaalne kirjete arv, mis peab igal patsiendil andmestikus olema.", "</div>")),
          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; text-align: justify;'>",
                      "<b> Trajektoori pikkus.</b> Võimaldab reguleerida, mitut ravimit patsiendi kohta vaadeldaks ehk kui pikk on iga patsiendi ravimitrajektoor.", "</div>")),
          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; text-align: justify;'>",
                      "<b> Klastrite arv.</b> Võimaldab reguleerida, mitmeks klastriks patsiendid jaotatakse.", "</div>")),
          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; text-align: justify;'>",
                      "<b> Vaadeldavate sagedaste ravimite arv.</b> Võimaldab reguleerida, mitut andmestikus enimesinenud ravimi kogust iga patsiendi puhul vaadatakse.", "</div>")),
          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; text-align: justify;'>",
                      "<b> DTW algoritm.</b> Võimaldab kasutada ravimitrajektooride võrdlemisel dünaamilise ajadeformatsiooni algoritmi. Vaikeväärtusena kasutatakse eukleidilise kauguse algoritmi.", "</div>")),
          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; text-align: justify;'>",
                      "<b> ATC koodide maksumused.</b> Võimaldab määrata maksumuse igale ATC koodi tasemele. Maksumused peavad olema mittenegatiivsed.", "</div>")),
          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; text-align: justify;'>",
                      "<b> Parameetrid.</b> Võimaldab alla laadida viimasena kasutatud parameetrid CSV-failina.", "</div>")),
          HTML(paste0("<div style='font-size: 15px; width: 95%; margin: 1%; text-align: justify;'>",
                      "<b> Diagramm \"Keskmise standardhälbe sõltuvus klastrite arvust\".</b> Kuvab, kuidas klastrite keskmine standardhälve sõltub klastrite arvust. See võimaldab kasutajal valida optimaalse klastrite arvu.", "</div>"))
        )
      )
    ),
    position = "right"
  )
)
