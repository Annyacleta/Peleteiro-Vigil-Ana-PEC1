
# Creación de estructura de datos básica

library(SummarizedExperiment)
library(tidyverse)

# Leo archivo completo como texto
file <- "ST003680_AN006041_metadata_data.txt"
raw_lines <- readLines(file)

# Localizo las secciones de datos
start_data <- grep("MS_METABOLITE_DATA_START", raw_lines) + 2
end_data   <- grep("MS_METABOLITE_DATA_END", raw_lines) - 1
start_meta_rt <- grep("METABOLITES_START", raw_lines) + 1
end_meta_rt   <- grep("METABOLITES_END", raw_lines) - 1

# Leer matriz de expresión (metabolitos × muestras)
expr_raw <- read.delim(file, skip = start_data - 1, nrows = 175, check.names = FALSE) # 175 metabolitos = filas
expr_raw_HEADER <- read.delim(file, skip = start_data - 2, nrows = 177, check.names = FALSE)

# Convertir a numérico y trasponer (muestras × metabolitos -> metabolitos x muestras)
expr_data <- apply(expr_raw[, -1], 2, as.numeric) %>% 
  t() %>%
  as.data.frame()

# Asignar nombres correctos
rownames(expr_data) <- colnames(expr_raw_HEADER)[-1]  # Nombres de muestras sin término "Samples" 
colnames(expr_data) <-   expr_raw[[1]] # Nombres de metabolitos sin término "Factor"
sample_names <- rownames(expr_data)

# COLDATA -- Metadatos de muestras
condition <- c(rep("KC", 8), rep("SC", 8)) #condiciones experimentales: keto y estándar

col_data <- data.frame(
  condition = condition,
  row.names = sample_names,
  stringsAsFactors = FALSE
)


# ROWDATA -- metadatos de metabolitos, como nombre y tiempo de retención
metabolite_meta <- read.delim(file,
                              skip = start_meta_rt - 1,
                              nrows = end_meta_rt - start_meta_rt,
                              stringsAsFactors = FALSE)

# Asigno nombres de columnas directamente
colnames(metabolite_meta) <- c("metabolite_name", "retention_time")

row_data <- metabolite_meta


# Creo la clase SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(intensity = t(expr_data)),  # Matriz numérica con todo el dataset
  rowData = DataFrame(row_data),   # Metabolitos en filas
  colData = DataFrame(col_data)    # Muestras en columnas
)



# Creación de lista de metadatos relativos al proyecto y condiciones del experimento

metadata <- list(
  source = "METABOLOMICS WORKBENCH",
  datatrack_id = "5543",
  study_id = "ST003680",
  analysis_id = "AN006041",
  project_id = "PR002283",
  version = "1",
  created_on = "January 21, 2025, 1:16 pm",
  
  project = list(
    title = "Ketogenic diet suppresses colorectal cancer through the gut microbiome long chain fatty acid stearate",
    summary = "Our manuscript entitled \"Ketogenic diet suppresses colorectal cancer through the gut microbiome long chain fatty acid stearate\" describes a reduced colonic tumor burden upon ketogenic diet (KD) consumption in a CRC mouse model with a humanized microbiome. Importantly, we demonstrate a causal relationship through microbiome transplantation into germ-free mice, whereby alterations in the gut microbiota were maintained in the absence of continued selective pressure from the KD. Specifically, we identify a shift toward bacterial species that produce stearic acid in ketogenic conditions, whereas consumers were depleted, resulting in elevated levels of free stearate in the gut lumen. This microbial product demonstrates tumor-suppressing properties by inducing apoptosis in cancer cells and decreasing colonic Th17 immune cell populations. As part of this study, we used different metabolomics workflows to study metabolites in mouse fecal and plasma samples as well as the used rodent diet.",
    institute = "University of Luxembourg",
    researcher = list(
      last_name = "Letellier",
      first_name = "Elisabeth"
    ),
    address = "6, avenue du Swing, Belval, Esch, 4367, Luxembourg",
    email = "madita.brauer@uni.lu",
    phone = "(+352) 46 66 44 6954"
  ),
  
  study = list(
    title = "Ketogenic diet suppresses colorectal cancer through the gut microbiome long chain fatty acid stearate - untargeted LCMS data from CMT experiment",
    summary = "Germ-free mice received cecal content from mice having received a ketogenic diet (KC) or from mice having received a standard diet (SC) and were subjected to AOM/DSS treatment as described in Tsenkova et al. (2025). Mouse fecal pellets were analyzed by a LC-MS approach to characterize the fecal metabolome of mice from different treatment groups. Metabolites were extracted from mouse fecal samples as described in the method file. Resulting samples were used for LC-MS untargeted metabolomics screen using a Vanquish UHPLC (ThermoFisher Scientific), coupled to a Q Exactive HF mass spectrometer (ThermoFisher Scientific).",
    institute = "University of Luxembourg",
    researcher = list(
      last_name = "Letellier",
      first_name = "Elisabeth"
    ),
    address = "6, avenue du Swing, Belval, Esch, 4367, Luxembourg",
    email = "madita.brauer@uni.lu",
    phone = "(+352) 46 66 44 6954"
  ),
  
  subject = list(
    subject_type = "Mammal",
    subject_species = "Mus musculus",
    taxonomy_id = "10090"
  ),
  
  subject_sample_factors = list(
    description = "SUBJECT(optional)[tab]SAMPLE[tab]FACTORS(NAME:VALUE pairs separated by |)[tab]Raw file names and additional sample data",
    samples = list(
      KC_16 = list(source = "mouse feces", condition = "KC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_KC_16_1.mzML"),
      KC_18 = list(source = "mouse feces", condition = "KC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_KC_18_1.mzML"),
      KC_23 = list(source = "mouse feces", condition = "KC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_KC_23_1.mzML"),
      KC_28 = list(source = "mouse feces", condition = "KC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_KC_28_1.mzML"),
      KC_34 = list(source = "mouse feces", condition = "KC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_KC_34_1.mzML"),
      KC_38 = list(source = "mouse feces", condition = "KC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_KC_38_1.mzML"),
      KC_42 = list(source = "mouse feces", condition = "KC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_KC_42_1.mzML"),
      KC_47 = list(source = "mouse feces", condition = "KC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_KC_47_1.mzML"),
      SC_11 = list(source = "mouse feces", condition = "SC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_SC_11_1.mzML"),
      SC_14 = list(source = "mouse feces", condition = "SC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_SC_14_1.mzML"),
      SC_20 = list(source = "mouse feces", condition = "SC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_SC_20_1.mzML"),
      SC_21 = list(source = "mouse feces", condition = "SC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_SC_21_1.mzML"),
      SC_30 = list(source = "mouse feces", condition = "SC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_SC_30_1.mzML"),
      SC_32 = list(source = "mouse feces", condition = "SC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_SC_32_1.mzML"),
      SC_35 = list(source = "mouse feces", condition = "SC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_SC_35_1.mzML"),
      SC_41 = list(source = "mouse feces", condition = "SC", raw_file = "MT011_LC_GlasgowMA_Mouse_feces_SC_41_1.mzML")
    )
  ),
  
  sampleprep = list(
    summary = "Mouse fecal pellets were placed in 0.5mL Precellys® tubes (VWR) containing five ceramic beads each and MilliQ® water was added to each sample at a 1:16 dry weight to water ratio. The samples were homogenized at 6000rpm, for two 30-second-long cycles at 4°C, in a Precellys®24 Homogenizer (Bertin Corp.). Samples were then incubated at 4°C for ten minutes, then centrifuged at maximum speed for 10 minutes at 4°C. Samples were maintained on ice and in the dark. The supernatant was used for further downstream processing. An ISM was prepared (2μg/mL of ribitol, pentanedioic-d6 acid and d-mannose and 10μg/mL tridecanoic-d25 acid, 6-chloropurine riboside, 4-chloro-DL-phenylalanine, Nε-trifluoroacetyl-L-lysine and thionicotinamide adenine dinucleotide (Sigma-Aldrich) in MilliQ® water). 40μL of the ISM was added to 100μL of the fecal supernatant fluid. 80μL of this mixture was added to 320μL of methanol, vortexed thoroughly, incubated for five minutes at 4°C at maximum speed in an Eppendorf ThermoMixer, and then centrifugated for five minutes at 4°C at maximum speed. 350μL of supernatant were added to 280μL of chloroform. 180μL of MilliQ® water were added, the samples were vortexed thoroughly, incubated for ten minutes at 4°C at maximum speed in an ThermoMixer (Eppendorf), then centrifugated for five minutes at 4°C at maximum speed. The extract was then split – 200μL of the upper (polar) phase were aliquoted into a GC vial with a micro-insert and the rest of the upper (polar) phase was filtered through a PHENEX-RC syringe filter (Phenomenex), and 200μL were retained in Eppendorf tubes for further processing. Then, the temperature of the SpeedVac® was increased to 25°C for 25 minutes to avoid water condensation on the surface of the glass vial. Samples were protected from light exposure and stored at -80°C. Samples for LC-MS analysis (polar) were reconstituted in 80μL 50% ACN in water, transferred into LC vials for LC-MS analysis and analyzed on a Vanquish UHPLC (ThermoFisher Scientific), coupled to a Q Exactive HF mass spectrometer (ThermoFisher Scientific)."
  ),
  
  chromatography = list(
    type = "HILIC",
    instrument_name = "Thermo Vanquish",
    column_name = "SeQuant ZIC- pHILIC (150 x 2.1mm,5um)",
    solvent_a = "100% water, 0.01% formic acid",
    solvent_b = "90% acetonitrile, 0.1% formic acid, 20 mM ammonium acetate",
    flow_gradient = "90 % B for 1.5 min, followed by a decrease to 20% B within 15 min, then 20% B for 2 min, increase to 90% B within 2 min, then 90% B for 13 min",
    flow_rate = "200 µl/ml",
    column_temperature = "50°C"
  ),
  
  analysis = list(
    type = "MS"
  ),
  
  ms = list(
    instrument_name = "Thermo Q Exactive HF-X Orbitrap",
    instrument_type = "Orbitrap",
    ms_type = "ESI",
    ion_mode = "UNSPECIFIED",
    comments = "Resolution - 120,000, m/z range - 60-900, AGC target - 1e6, maximum injection time – 70 ms. ddMS2 was applied using the following settings: Resolution - 30,000, AGC target - 5e5, maximum injection time – 70, topN - 5, Normalized collision energy – 20. Samples were acquired in positive and negative ionization mode simultaneously (polarity switching)."
  )
)

metadata(se) <- metadata



# Guardo objeto SummarizeExperiment en formato binario .Rda
save(se, file = "SummarizedExperiment_metabolomics.Rda")

