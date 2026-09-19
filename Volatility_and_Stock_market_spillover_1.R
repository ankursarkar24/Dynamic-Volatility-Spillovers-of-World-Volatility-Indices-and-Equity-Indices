library(quantmod)

# 1. Define Tickers
group1_tickers <- c("VIX" = "^VIX", "OVX" = "^OVX", "GVZ" = "^GVZ")
group2_tickers <- c("SP500"  = "^GSPC",
                    "FTSE"   = "^FTSE",
                    "DAX"    = "^GDAXI",
                    "NIKKEI" = "^N225",
                    "NIFTY"  = "^NSEI",
                    "HSI"    = "^HSI",
                    "TSX"    = "^GSPTSE")

all_tickers <- c(group1_tickers, group2_tickers)

# 2. Historical Window
start_date <- "2016-01-01"
end_date   <- "2026-06-30"

price_list <- list()

for (sym_name in names(all_tickers)) {
  ticker_sym <- all_tickers[sym_name]
  
  # Download series
  raw_obj <- getSymbols(ticker_sym, src = "yahoo", 
                        from = start_date, to = end_date, auto.assign = FALSE)
  
  # Extract Adjusted Close
  adj_close <- Ad(raw_obj)
  
  # Fill internal single-day missing values if any exist within the ticker's own series
  adj_close <- na.locf(adj_close, na.rm = TRUE)
  
  price_list[[sym_name]] <- adj_close
}

# 3. Merge into a common matrix
merged_raw <- do.call(merge, price_list)

# 4. Handle cross-exchange holidays:
# Forward-fill remaining holiday gaps across mismatched national calendars, 
# then drop leading NAs where series didn't start on the exact same calendar day
aligned_prices <- na.omit(na.locf(merged_raw))
colnames(aligned_prices) <- names(all_tickers)

# 5. Volatility Construction: VOL_{i,t} = (100 * ln(P_t / P_{t-1}))^2
returns_pct <- 100 * diff(log(aligned_prices))[-1, ]
vol_matrix  <- returns_pct^2
vol_df      <- as.data.frame(vol_matrix)

# 6. Verification check (Must output 0)
cat("Any remaining missing values:", sum(is.na(vol_df)), "\n")
cat("Cleaned sample dimensions:", nrow(vol_df), "rows x", ncol(vol_df), "columns\n")





# ==============================================================================
# Table 2: Descriptive Statistics (Mean, Var, Skew, Kurt, JB, Q(10))
# ==============================================================================
table2 <- data.frame(
  Mean        = colMeans(vol_df),
  Variance    = apply(vol_df, 2, var),
  Skewness    = apply(vol_df, 2, moments::skewness),
  Ex_Kurtosis = apply(vol_df, 2, moments::kurtosis) - 3,
  JB_Stat     = apply(vol_df, 2, function(x) jarque.bera.test(x)$statistic),
  JB_p        = apply(vol_df, 2, function(x) jarque.bera.test(x)$p.value),
  Q10_Stat    = apply(vol_df, 2, function(x) Box.test(x, lag = 10, type = "Ljung-Box")$statistic),
  Q10_p       = apply(vol_df, 2, function(x) Box.test(x, lag = 10, type = "Ljung-Box")$p.value)
)
print("=== Table 2: Descriptive Statistics ===")
print(round(table2, 3))

# ==============================================================================
# Table 3: Bivariate Pearson Correlation Matrix
# ==============================================================================
table3 <- cor(vol_df)
print("=== Table 3: Correlation Matrix ===")
print(round(table3, 3))

# ==============================================================================
# Table 4: Unit Root Tests (ADF & ERS Point Optimal Tests)
# ==============================================================================
library(urca)

# 1. ADF Test with constant/drift
adf_vals <- sapply(vol_df, function(x) {
  test <- ur.df(x, type = "drift", selectlags = "BIC")
  test@teststat[1] # Extracts tau2 statistic for H0: gamma = 0
})

# 2. ERS Point Optimal Test (model = "constant", type = "P-test")
ers_vals <- sapply(vol_df, function(x) {
  test <- ur.ers(x, type = "P-test", model = "constant")
  test@teststat[1] # Extracts P-statistic
})

# 3. Assemble Table 4
table4 <- data.frame(
  ADF_Statistic = adf_vals,
  ADF_Crit_5pct = -2.86,
  ERS_P_Stat    = ers_vals,
  ERS_Crit_5pct = 3.26
)

cat("=== Table 4: Unit Root Tests ===\n")
print(round(table4, 3))


library(urca)

# 1. Run ADF test (drift / intercept specification, lag selection via BIC/SC)
adf_models <- lapply(vol_df, function(x) {
  ur.df(x, type = "drift", selectlags = "BIC")
})

# Extract tau2 test statistic (drift specification) and its exact 5% critical value
adf_stat <- sapply(adf_models, function(m) m@teststat[1, "tau2"])
adf_crit_5pct <- sapply(adf_models, function(m) m@cval["tau2", "5pct"])

# 2. Run ERS Point Optimal test (constant / intercept specification)
ers_models <- lapply(vol_df, function(x) {
  ur.ers(x, type = "P-test", model = "constant")
})

# Extract ERS Point Optimal P-statistic and its exact 5% critical value
ers_stat <- sapply(ers_models, function(m) m@teststat[1])
ers_crit_5pct <- sapply(ers_models, function(m) m@cval[1, "5pct"])

# 3. Assemble dynamic Table 4 (matching Table 4 of Koç, 2026)
table4 <- data.frame(
  ADF_Test_Stat = adf_stat,
  ADF_5pct_Crit = adf_crit_5pct,
  ERS_P_Stat    = ers_stat,
  ERS_5pct_Crit = ers_crit_5pct
)

cat("=== Table 4: Dynamic Unit Root Test Results ===\n")
print(round(table4, 4))


library(ConnectednessApproach)

# 1. Optimal specification matching Koç (2026) Section 3.2
dca <- ConnectednessApproach(
  vol_df,
  nlag   = 1,             # Lag length p = 1 per Schwarz Criterion (SC)
  nfore  = 10,            # 10-step ahead GFEVD horizon
  model  = "TVP-VAR",     # Passes model type directly
  Connectedness_Method = "time",
  corrected = FALSE       # Baseline TCI (set TRUE for corrected cTCI)
)

# 2. View Table 5 (Connectedness Matrix)
print("=== Table 5: Dynamic Connectedness Matrix ===")
print(dca$TABLE)





library(ConnectednessApproach)
library(zoo)

# 1. Convert volatility data to a zoo object using index dates
vol_zoo <- as.zoo(vol_df)

# 2. TVP-VAR Estimation matching Koç (2026) Section 3.2
dca <- ConnectednessApproach(
  vol_zoo,
  nlag          = 1,                  # p = 1 per Schwarz Criterion (SC)
  nfore         = 10,                 # Forecast horizon H = 10 days
  model         = "TVP-VAR",
  connectedness = "Time",             # Correct argument name (capital 'T')
  corrected     = FALSE,              # Baseline TCI (use TRUE for cTCI)
  VAR_config    = list(
    TVPVAR = list(
      kappa1 = 0.99,                  # Coefficient forgetting factor
      kappa2 = 0.99,                  # Covariance forgetting factor
      prior  = "BayesPrior"           # Prior specification
    )
  )
)

# 3. Display Table 5: Dynamic Connectedness Matrix
cat("\n=== Table 5: Dynamic Connectedness Matrix ===\n")
print(dca$TABLE)




# Define the two sub-systems
groups <- list(
  "Vol_Indices" = c(1, 2, 3), # VIX, OVX, GVZ
  "Stock_Markets" = c(4:10)   # SP500, FTSE, DAX, NIKKEI, NIFTY, HSI, TSX
)

# Decompose Spillovers
internal_dca <- InternalConnectedness(dca, groups = groups)
external_dca <- ExternalConnectedness(dca, groups = groups)

# Report Overall TCI Decomposition
cat("\n--- Connectedness Breakdown ---\n")
cat("Total System TCI       :", round(mean(dca$TCI), 2), "%\n")
cat("Internal TCI (Within)  :", round(mean(internal_dca$TCI), 2), "%\n")
cat("External TCI (Between) :", round(mean(external_dca$TCI), 2), "%\n")






# Figure 1: Dynamic Total vs. Group Connectedness
PlotTCI(dca)

# Figure 2: Dynamic External (Cross-Group) Connectedness
PlotTCI(external_dca)

# Figure 3: Dynamic Total Net Directional Spillovers
PlotNET(dca)
\
# Figure 7: Average Net Pairwise Directional Connectedness (NPDC) Network
PlotNetwork(dca)
















library(ConnectednessApproach)
library(ggplot2)
library(zoo)

# 1. Define groups matching the 10-variable matrix
groups <- list(
  "Volatility_Indices" = c(1, 2, 3),
  "Stock_Markets"      = c(4:10)
)

# 2. Decompose into Internal and External Connectedness (Gabauer & Gupta, 2018)
dca_internal <- InternalConnectedness(dca, groups = groups)
dca_external <- ExternalConnectedness(dca, groups = groups)


names(dca_internal)
dim(dca_internal$TCI)
colnames(dca_internal$TCI)

library(ggplot2)
library(zoo)
library(ConnectednessApproach)

# 1. Recover true historical dates from your returns/aligned_prices object
real_dates <- as.Date(index(returns_pct))

# 2. Check if lengths match
if (length(real_dates) != nrow(vol_df)) {
  # If row count was reduced by nlag = 1 during estimation:
  real_dates <- tail(real_dates, length(dca$TCI))
}

# 3. Compute Group-Specific within-TCI (Group 1 vs Group 2)
# Group 1: VIX, OVX, GVZ (columns 1 to 3)
vol_zoo_g1 <- as.zoo(vol_df[, 1:3], order.by = real_dates)
dca_g1 <- ConnectednessApproach(
  vol_zoo_g1, nlag = 1, nfore = 10, model = "TVP-VAR", 
  connectedness = "Time", corrected = FALSE
)

# Group 2: 7 Stock Market Indices (columns 4 to 10)
vol_zoo_g2 <- as.zoo(vol_df[, 4:10], order.by = real_dates)
dca_g2 <- ConnectednessApproach(
  vol_zoo_g2, nlag = 1, nfore = 10, model = "TVP-VAR", 
  connectedness = "Time", corrected = FALSE
)

# 4. Construct df_fig1 with real calendar dates and zero NAs
df_fig1 <- data.frame(
  Date          = real_dates,
  Total_TCI     = as.numeric(dca$TCI),
  Internal_TCI  = as.numeric(dca_internal$TCI),
  Vol_Indices   = as.numeric(dca_g1$TCI),
  Stock_Markets = as.numeric(dca_g2$TCI)
)

cat("Successfully rebuilt df_fig1 with real calendar dates:\n")
head(df_fig1)

p_fig1 <- ggplot(df_fig1, aes(x = Date)) +
  geom_area(aes(y = Total_TCI, fill = "Total TCI"), alpha = 0.9) +
  geom_line(aes(y = Internal_TCI, color = "Internal TCI"), linewidth = 0.8) +
  geom_line(aes(y = Vol_Indices, color = "Volatility Indices"), linewidth = 0.8) +
  geom_line(aes(y = Stock_Markets, color = "Stock Markets"), linewidth = 0.8) +
  scale_fill_manual(name = "", values = c("Total TCI" = "black")) +
  scale_color_manual(name = "", values = c(
    "Internal TCI"       = "#E41A1C", 
    "Volatility Indices" = "#4DAF4A", 
    "Stock Markets"      = "#377EB8"
  )) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = "", y = "TCI (%)", title = "Figure 1. Dynamic Total, Internal, and Group-Specific Connectedness") +
  theme_classic(base_size = 12) +
  theme(
    legend.position   = c(0.22, 0.82),
    legend.background = element_rect(fill = "transparent"),
    axis.line         = element_line(color = "black")
  )

print(p_fig1)





# Extract external TCI
ext_tci_vlibrary(ggplot2)

# 1. Build df_fig2 using real_dates
df_fig2 <- data.frame(
  Date         = real_dates,
  External_TCI = as.numeric(dca_external$TCI)
)

# 2. Identify global peak
peak_val   <- max(df_fig2$External_TCI, na.rm = TRUE)
peak_date  <- df_fig2$Date[which.max(df_fig2$External_TCI)]
peak_label <- paste0("Peak: ", round(peak_val, 2), "% (", format(peak_date, "%d %b %Y"), ")")

# 3. Figure 2 Plot (Warning-Free, Publication Format)
p_fig2 <- ggplot(df_fig2, aes(x = Date, y = External_TCI)) +
  geom_line(color = "black", linewidth = 0.7) +
  # Event markers (vertical dashed lines matching paper's stress regimes)
  geom_vline(xintercept = as.Date("2020-03-11"), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = as.Date("2022-02-24"), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = as.Date("2022-03-16"), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = as.Date("2025-01-20"), linetype = "dashed", color = "grey40") +
  # Peak point and annotation via annotate() to avoid row-recycling warnings
  annotate("point", x = peak_date, y = peak_val, color = "black", size = 2) +
  annotate("text", x = peak_date, y = peak_val + (max(df_fig2$External_TCI) * 0.05), 
           label = peak_label, size = 3.5, fontface = "bold") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = "", y = "External TCI (%)", title = "Figure 2. Dynamic Cross-Group Connectedness") +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold")
  )

print(p_fig2)








calc_regime_stats <- function(df, start_d, end_d, regime_name) {
  sub_df <- subset(df, Date >= as.Date(start_d) & Date <= as.Date(end_d))
  if (nrow(sub_df) == 0) return(NULL)
  
  mean_val <- mean(sub_df$External_TCI, na.rm = TRUE)
  max_val  <- max(sub_df$External_TCI, na.rm = TRUE)
  peak_d   <- sub_df$Date[which.max(sub_df$External_TCI)]
  
  data.frame(
    Period            = regime_name,
    Mean_External_TCI = round(mean_val, 2),
    Max_External_TCI  = round(max_val, 2),
    Peak_Date         = format(peak_d, "%d %B %Y")
  )
}

table6 <- rbind(
  calc_regime_stats(df_fig2, "2016-01-01", "2019-12-31", "Pre-COVID"),
  calc_regime_stats(df_fig2, "2020-01-01", "2021-12-31", "COVID-19"),
  calc_regime_stats(df_fig2, "2022-01-01", "2022-12-31", "Russia-Ukraine & Tightening"),
  calc_regime_stats(df_fig2, "2023-01-01", "2024-12-31", "Post-Crisis Expansion"),
  calc_regime_stats(df_fig2, "2025-01-01", "2026-06-30", "Recent Macro & Tariff Uncertainty")
)

cat("\n=== Table 6: External Connectedness Across Selected Market Regimes ===\n")
print(table6, row.names = FALSE)


net_total    <- dca$NET
net_internal <- dca_internal$NET

var_names <- colnames(net_total)
fig3_list <- list()

for (v in var_names) {
  tmp <- data.frame(
    Date     = real_dates,
    Variable = v,
    Total    = as.numeric(net_total[, v]),
    Internal = as.numeric(net_internal[, v])
  )
  fig3_list[[v]] <- tmp
}

df_fig3 <- do.call(rbind, fig3_list)
df_fig3$Variable <- factor(df_fig3$Variable, levels = var_names)

# Figure 3 Plot
p_fig3 <- ggplot(df_fig3, aes(x = Date)) +
  geom_area(aes(y = Total), fill = "black", alpha = 0.85) +
  geom_line(aes(y = Internal), color = "red", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  facet_wrap(~ Variable, ncol = 3, scales = "free_y") +
  labs(x = "", y = "NET Connectedness", 
       title = "Figure 3. Dynamic Total and Within-Group Net Directional Connectedness") +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_fig3)






library(ggplot2)

# 1. Extract NET series from baseline and internal decomposition
net_total    <- dca$NET
net_internal <- dca_internal$NET

var_names <- colnames(net_total)
fig3_list <- list()

for (v in var_names) {
  tmp <- data.frame(
    Date     = real_dates,
    Variable = v,
    Total    = as.numeric(net_total[, v]),
    Internal = as.numeric(net_internal[, v])
  )
  fig3_list[[v]] <- tmp
}

df_fig3 <- do.call(rbind, fig3_list)
df_fig3$Variable <- factor(df_fig3$Variable, levels = var_names)

# 2. Determine common y-axis boundaries with slight padding
y_min <- min(c(df_fig3$Total, df_fig3$Internal), na.rm = TRUE)
y_max <- max(c(df_fig3$Total, df_fig3$Internal), na.rm = TRUE)
y_pad <- (y_max - y_min) * 0.05

# 3. Figure 3 Plot with Fixed/Common Y-Axis Scale
p_fig3 <- ggplot(df_fig3, aes(x = Date)) +
  geom_area(aes(y = Total), fill = "black", alpha = 0.85) +
  geom_line(aes(y = Internal), color = "red", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.4) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(limits = c(y_min - y_pad, y_max + y_pad)) +
  # Using scales = "fixed" enforces the exact same y-scale across all 10 panels
  facet_wrap(~ Variable, ncol = 3, scales = "fixed") +
  labs(
    x = "", 
    y = "NET Connectedness", 
    title = "Figure 3. Dynamic Total and Within-Group Net Directional Connectedness"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_fig3)















library(ggplot2)

plot_npdc_figure <- function(dca_obj, base_var, target_vars, fig_title) {
  # Extract NPDC 3D array: dim is [Time, Variables, Variables]
  npdc_arr <- dca_obj$NPDC
  
  plot_list <- list()
  for (tgt in target_vars) {
    # NPDC from base_var to tgt
    # In Gabauer's package: NPDC[t, i, j] = phi[t, j, i] - phi[t, i, j]
    val_series <- as.numeric(npdc_arr[, base_var, tgt])
    
    tmp <- data.frame(
      Date = real_dates,
      Pair = paste0(base_var, " - ", tgt),
      NPDC = val_series
    )
    plot_list[[tgt]] <- tmp
  }
  
  df_plot <- do.call(rbind, plot_list)
  # Preserve the exact target ordering
  df_plot$Pair <- factor(df_plot$Pair, levels = paste0(base_var, " - ", target_vars))
  
  p <- ggplot(df_plot, aes(x = Date, y = NPDC)) +
    geom_area(fill = "black", alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    facet_wrap(~ Pair, ncol = 2, scales = "free_y") +
    labs(x = "", y = "NPDC", title = fig_title) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# The 7 target stock markets
stock_targets <- c("SP500", "FTSE", "DAX", "NIKKEI", "NIFTY", "HSI", "TSX")

dim(dca$NPDC)
dimnames(dca$NPDC)


library(ggplot2)

plot_npdc_figure <- function(dca_obj, base_var, target_vars, fig_title) {
  npdc_arr <- dca_obj$NPDC
  all_names <- colnames(vol_df)
  
  # Find integer index of base_var
  base_idx <- which(all_names == base_var)
  if (length(base_idx) == 0) stop(paste("Base variable", base_var, "not found in dataset"))
  
  # Determine if Time is dimension 1 or dimension 3
  dims <- dim(npdc_arr)
  time_in_dim1 <- (dims[1] == length(real_dates))
  
  plot_list <- list()
  for (tgt in target_vars) {
    tgt_idx <- which(all_names == tgt)
    if (length(tgt_idx) == 0) stop(paste("Target variable", tgt, "not found in dataset"))
    
    # Extract pairwise series: Base -> Target
    if (time_in_dim1) {
      val_series <- as.numeric(npdc_arr[, base_idx, tgt_idx])
    } else {
      val_series <- as.numeric(npdc_arr[base_idx, tgt_idx, ])
    }
    
    tmp <- data.frame(
      Date = real_dates,
      Pair = paste0(base_var, " - ", tgt),
      NPDC = val_series
    )
    plot_list[[tgt]] <- tmp
  }
  
  df_plot <- do.call(rbind, plot_list)
  df_plot$Pair <- factor(df_plot$Pair, levels = paste0(base_var, " - ", target_vars))
  
  p <- ggplot(df_plot, aes(x = Date, y = NPDC)) +
    geom_area(fill = "black", alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    facet_wrap(~ Pair, ncol = 2, scales = "free_y") +
    labs(x = "", y = "NPDC", title = fig_title) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# The 7 target stock markets
stock_targets <- c("SP500", "FTSE", "DAX", "NIKKEI", "NIFTY", "HSI", "TSX")


# Figure 4: VIX to 7 Stock Markets
p_fig4 <- plot_npdc_figure(
  dca_obj     = dca,
  base_var    = "VIX",
  target_vars = stock_targets,
  fig_title   = "Figure 4. NPDC: VIX to International Stock Markets"
)
print(p_fig4)

# Figure 5: OVX to 7 Stock Markets
p_fig5 <- plot_npdc_figure(
  dca_obj     = dca,
  base_var    = "OVX",
  target_vars = stock_targets,
  fig_title   = "Figure 5. NPDC: OVX to International Stock Markets"
)
print(p_fig5)

# Figure 6: GVZ to 7 Stock Markets
p_fig6 <- plot_npdc_figure(
  dca_obj     = dca,
  base_var    = "GVZ",
  target_vars = stock_targets,
  fig_title   = "Figure 6. NPDC: GVZ to International Stock Markets"
)
print(p_fig6)









library(ggplot2)

plot_npdc_figure <- function(dca_obj, base_var, target_vars, fig_title) {
  npdc_arr <- dca_obj$NPDC
  all_names <- colnames(vol_df)
  
  # Find integer index of base_var
  base_idx <- which(all_names == base_var)
  if (length(base_idx) == 0) stop(paste("Base variable", base_var, "not found in dataset"))
  
  # Determine if Time is dimension 1 or dimension 3
  dims <- dim(npdc_arr)
  time_in_dim1 <- (dims[1] == length(real_dates))
  
  plot_list <- list()
  for (tgt in target_vars) {
    tgt_idx <- which(all_names == tgt)
    if (length(tgt_idx) == 0) stop(paste("Target variable", tgt, "not found in dataset"))
    
    # Extract pairwise series: Base -> Target
    if (time_in_dim1) {
      val_series <- as.numeric(npdc_arr[, base_idx, tgt_idx])
    } else {
      val_series <- as.numeric(npdc_arr[base_idx, tgt_idx, ])
    }
    
    tmp <- data.frame(
      Date = real_dates,
      Pair = paste0(base_var, " - ", tgt),
      NPDC = val_series
    )
    plot_list[[tgt]] <- tmp
  }
  
  df_plot <- do.call(rbind, plot_list)
  df_plot$Pair <- factor(df_plot$Pair, levels = paste0(base_var, " - ", target_vars))
  
  # Calculate 5-unit rounded limits across all panels
  raw_min <- min(df_plot$NPDC, na.rm = TRUE)
  raw_max <- max(df_plot$NPDC, na.rm = TRUE)
  
  y_min_5 <- floor(raw_min / 5) * 5
  y_max_5 <- ceiling(raw_max / 5) * 5
  
  p <- ggplot(df_plot, aes(x = Date, y = NPDC)) +
    geom_area(fill = "black", alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    # Fixed shared y-axis with tick breaks every 5 units
    scale_y_continuous(
      limits = c(y_min_5, y_max_5),
      breaks = seq(y_min_5, y_max_5, by = 5)
    ) +
    facet_wrap(~ Pair, ncol = 2, scales = "fixed") +
    labs(x = "", y = "NPDC", title = fig_title) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# The 7 target stock markets
stock_targets <- c("SP500", "FTSE", "DAX", "NIKKEI", "NIFTY", "HSI", "TSX")


# Figure 4: VIX to 7 Stock Markets
p_fig4 <- plot_npdc_figure(
  dca_obj     = dca,
  base_var    = "VIX",
  target_vars = stock_targets,
  fig_title   = "Figure 4. NPDC: VIX to International Stock Markets"
)
print(p_fig4)

# Figure 5: OVX to 7 Stock Markets
p_fig5 <- plot_npdc_figure(
  dca_obj     = dca,
  base_var    = "OVX",
  target_vars = stock_targets,
  fig_title   = "Figure 5. NPDC: OVX to International Stock Markets"
)
print(p_fig5)

# Figure 6: GVZ to 7 Stock Markets
p_fig6 <- plot_npdc_figure(
  dca_obj     = dca,
  base_var    = "GVZ",
  target_vars = stock_targets,
  fig_title   = "Figure 6. NPDC: GVZ to International Stock Markets"
)
print(p_fig6)







library(ConnectednessApproach)
library(igraph)

# Option A: Gabauer's built-in network visualizer
# Setting method = "NPDC" renders net directional pairwise arrows
PlotNetwork(dca, method = "NPDC")

# ------------------------------------------------------------------------------
# Option B: Fully Styled igraph Matching Figure 7 of the Paper
# ------------------------------------------------------------------------------
# 1. Extract average NPDC matrix and average NET / TO measures
# Average over the time dimension of NPDC array
dims <- dim(dca$NPDC)
time_in_dim1 <- (dims[1] == length(real_dates))

if (time_in_dim1) {
  avg_npdc <- apply(dca$NPDC, c(2, 3), mean, na.rm = TRUE)
} else {
  avg_npdc <- apply(dca$NPDC, c(1, 2), mean, na.rm = TRUE)
}

colnames(avg_npdc) <- colnames(vol_df)
rownames(avg_npdc) <- colnames(vol_df)

# Retain only positive net flows (from transmitter i -> receiver j)
net_adj <- avg_npdc
net_adj[net_adj < 0] <- 0

# 2. Build igraph object
g <- graph_from_adjacency_matrix(net_adj, mode = "directed", weighted = TRUE, diag = FALSE)

# 3. Calculate average TO and NET from baseline estimation
avg_net <- colMeans(dca$NET, na.rm = TRUE)
avg_to  <- colMeans(dca$TO, na.rm = TRUE)

# Color assignment: Blue = Net Transmitter (>0), Yellow = Net Receiver (<0)
V(g)$color <- ifelse(avg_net[V(g)$name] > 0, "#2B6CB0", "#E6AF2E")
V(g)$frame.color <- "grey40"
V(g)$label.color <- "black"
V(g)$label.font  <- 2
V(g)$label.cex   <- 0.85

# Node sizing proportional to TO spillovers
V(g)$size <- 12 + 18 * (avg_to[V(g)$name] - min(avg_to)) / (max(avg_to) - min(avg_to) + 1e-5)

# Edge styling (thickness proportional to average NPDC strength)
edge_w <- E(g)$weight
E(g)$width <- 0.5 + 4 * (edge_w - min(edge_w)) / (max(edge_w) - min(edge_w) + 1e-5)
E(g)$arrow.size <- 0.35
E(g)$arrow.width <- 0.8
E(g)$color <- rgb(0.3, 0.3, 0.3, 0.5)

# Plot circular layout matching Figure 7
par(mar = c(1, 1, 2, 1))
plot(g, layout = layout_in_circle(g), 
     main = "Figure 7. Average NPDC Network of Global Volatility Indices and International Stock Markets")
legend("bottomright", legend = c("Net Transmitter", "Net Receiver"), 
       fill = c("#2B6CB0", "#E6AF2E"), bty = "n", cex = 0.85)


library(ConnectednessApproach)

cat("\n--- Running Robustness Models (Table 7) ---\n")

# 1. Alternative Forecast Horizons: H = 5 and H = 20
dca_h5  <- ConnectednessApproach(vol_zoo, nlag = 1, nfore = 5,  model = "TVP-VAR", connectedness = "Time", corrected = FALSE)
dca_h20 <- ConnectednessApproach(vol_zoo, nlag = 1, nfore = 20, model = "TVP-VAR", connectedness = "Time", corrected = FALSE)

# 2. Alternative Lag Specification: Lag = 2
dca_p2  <- ConnectednessApproach(vol_zoo, nlag = 2, nfore = 10, model = "TVP-VAR", connectedness = "Time", corrected = FALSE)

# 3. Standard Diebold-Yilmaz (DY 2012) Rolling-Window VAR
dca_dy200 <- ConnectednessApproach(vol_zoo, nlag = 1, nfore = 10, window.size = 200, model = "VAR", connectedness = "Time", corrected = FALSE)
dca_dy250 <- ConnectednessApproach(vol_zoo, nlag = 1, nfore = 10, window.size = 250, model = "VAR", connectedness = "Time", corrected = FALSE)

# ------------------------------------------------------------------------------
# Extract Metrics to Assemble Table 7
# ------------------------------------------------------------------------------
extract_row_table7 <- function(model_obj, model_label, base_tci_vector) {
  tci_ser <- as.numeric(model_obj$TCI)
  
  # Dates of this model's TCI
  m_dates <- index(model_obj$TCI)
  if (is.null(m_dates) || is.numeric(m_dates)) {
    # If index lost, align from tail of real_dates
    m_dates <- tail(real_dates, length(tci_ser))
  }
  
  # Net measures across all variables
  net_means <- colMeans(model_obj$NET, na.rm = TRUE)
  transmitters <- names(sort(net_means[net_means > 0], decreasing = TRUE))
  receivers    <- names(sort(net_means[net_means < 0], decreasing = FALSE))
  
  # Correlation with baseline TCI (over overlapping periods)
  n_len <- length(tci_ser)
  base_sub <- tail(base_tci_vector, n_len)
  cor_val <- round(cor(tci_ser, base_sub, use = "complete.obs"), 4)
  
  # Max TCI and peak date
  max_tci_val <- max(tci_ser, na.rm = TRUE)
  peak_date   <- m_dates[which.max(tci_ser)]
  
  data.frame(
    Measure                 = model_label,
    Average_TCI             = round(mean(tci_ser, na.rm = TRUE), 2),
    Main_Net_Transmitters   = paste(transmitters, collapse = ", "),
    Main_Net_Receivers      = paste(receivers, collapse = ", "),
    TCI_Corr_With_Baseline  = cor_val,
    Maximum_TCI             = round(max_tci_val, 2),
    Peak_Date               = format(as.Date(peak_date), "%d %B %Y")
  )
}

base_tci_vec <- as.numeric(dca$TCI)

# Build all 6 columns/rows of Table 7
t7_baseline <- extract_row_table7(dca,       "Baseline TVP-VAR", base_tci_vec)
t7_h5       <- extract_row_table7(dca_h5,    "H = 5",            base_tci_vec)
t7_h20      <- extract_row_table7(dca_h20,   "H = 20",           base_tci_vec)
t7_p2       <- extract_row_table7(dca_p2,    "Lag = 2",          base_tci_vec)
t7_dy200    <- extract_row_table7(dca_dy200, "DY VAR (200-Day)", base_tci_vec)
t7_dy250    <- extract_row_table7(dca_dy250, "DY VAR (250-Day)", base_tci_vec)

table7_raw <- rbind(t7_baseline, t7_h5, t7_h20, t7_p2, t7_dy200, t7_dy250)

# Transpose to match the vertical column structure of Table 7 in Koç (2026)
table7_transposed <- as.data.frame(t(table7_raw[, -1]))
colnames(table7_transposed) <- table7_raw$Measure

cat("\n=== Table 7: Robustness Checks across Alternative Specifications ===\n")
print(table7_transposed)










# ==============================================================================
# Replicating Table 7 Formatting (Koç, 2026 Style)
# ==============================================================================

format_table7_paper_style <- function(dca, dca_h5, dca_h20, dca_p2, dca_dy200, dca_dy250, base_dates) {
  
  # List all 6 model specifications
  models <- list(
    "Baseline TVP-VAR"   = dca,
    "H = 5"              = dca_h5,
    "H = 20"             = dca_h20,
    "Lag = 2"            = dca_p2,
    "DY VAR (200-Day)"   = dca_dy200,
    "DY VAR (250-Day)"   = dca_dy250
  )
  
  col_names <- names(models)
  
  # 1. Baseline vectors for comparison
  base_tci <- as.numeric(dca$TCI)
  base_net <- colMeans(dca$NET, na.rm = TRUE)
  base_tx  <- names(sort(base_net[base_net > 0], decreasing = TRUE))
  base_rx  <- names(sort(base_net[base_net < 0], decreasing = FALSE))
  
  # Placeholders for each row
  row_tci          <- character(6)
  row_transmitters <- character(6)
  row_receivers    <- character(6)
  row_changes      <- character(6)
  row_corr         <- character(6)
  row_max_tci      <- character(6)
  row_peak_date    <- character(6)
  
  for (i in seq_along(models)) {
    m_name <- col_names[i]
    m_obj  <- models[[i]]
    
    tci_vec <- as.numeric(m_obj$TCI)
    m_dates <- index(m_obj$TCI)
    if (is.null(m_dates) || is.numeric(m_dates)) {
      m_dates <- tail(base_dates, length(tci_vec))
    }
    
    # Average TCI
    row_tci[i] <- sprintf("%.2f", mean(tci_vec, na.rm = TRUE))
    
    # Net Transmitters & Receivers
    m_net <- colMeans(m_obj$NET, na.rm = TRUE)
    cur_tx <- names(sort(m_net[m_net > 0], decreasing = TRUE))
    cur_rx <- names(sort(m_net[m_net < 0], decreasing = FALSE))
    
    if (i == 1) {
      row_transmitters[i] <- paste(cur_tx, collapse = ", ")
      row_receivers[i]    <- paste(cur_rx, collapse = ", ")
      row_changes[i]      <- "–"
      row_corr[i]         <- "1.000"
      row_max_tci[i]      <- sprintf("%.2f", max(tci_vec, na.rm = TRUE))
      row_peak_date[i]    <- format(as.Date(m_dates[which.max(tci_vec)]), "%d %B %Y")
    } else {
      # Shorthand conventions used in Koç (2026) Table 7:
      # Use "Same" if identical set
      row_transmitters[i] <- ifelse(setequal(cur_tx, base_tx), "Same", paste(cur_tx, collapse = ", "))
      row_receivers[i]    <- ifelse(setequal(cur_rx, base_rx), "Same", paste(cur_rx, collapse = ", "))
      
      # Role changes (number of assets whose net sign flipped vs baseline)
      sign_base <- sign(base_net)
      sign_cur  <- sign(m_net[names(base_net)])
      flips     <- names(which(sign_base != sign_cur))
      row_changes[i] <- if (length(flips) == 0) "0" else paste0(length(flips), " (", paste(flips, collapse = ", "), ")")
      
      # Table 7 leaves dashed entries for H = 5 and H = 20 on secondary metrics
      if (m_name %in% c("H = 5", "H = 20")) {
        row_corr[i]      <- "–"
        row_max_tci[i]   <- "–"
        row_peak_date[i] <- "–"
      } else {
        # Exact date-aligned correlation
        n_overlap   <- min(length(tci_vec), length(base_tci))
        sub_base    <- tail(base_tci, n_overlap)
        sub_cur     <- tail(tci_vec, n_overlap)
        row_corr[i] <- sprintf("%.4f", cor(sub_cur, sub_base, use = "complete.obs"))
        
        row_max_tci[i]   <- sprintf("%.2f", max(tci_vec, na.rm = TRUE))
        row_peak_date[i] <- format(as.Date(m_dates[which.max(tci_vec)]), "%d %B %Y")
      }
    }
  }
  
  # Assemble final publication matrix
  final_tab <- data.frame(
    "Baseline TVP-VAR" = character(7),
    "H = 5"            = character(7),
    "H = 20"           = character(7),
    "Lag = 2"          = character(7),
    "DY VAR (200-Day)" = character(7),
    "DY VAR (250-Day)" = character(7),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  final_tab[1, ] <- row_tci
  final_tab[2, ] <- row_transmitters
  final_tab[3, ] <- row_receivers
  final_tab[4, ] <- row_changes
  final_tab[5, ] <- row_corr
  final_tab[6, ] <- row_max_tci
  final_tab[7, ] <- row_peak_date
  
  rownames(final_tab) <- c(
    "Average TCI (%)",
    "Main net transmitters",
    "Main net receivers",
    "Transmitter/receiver role changes",
    "Dynamic TCI correlation with baseline",
    "Maximum TCI (%)",
    "Peak date"
  )
  
  return(final_tab)
}

# Run the function with your estimated models and real calendar dates
table7_publication <- format_table7_paper_style(
  dca         = dca,
  dca_h5      = dca_h5,
  dca_h20     = dca_h20,
  dca_p2      = dca_p2,
  dca_dy200   = dca_dy200,
  dca_dy250   = dca_dy250,
  base_dates  = real_dates
)

# Display Table 7 exactly as it appears in Koç (2026)
print(as.matrix(table7_publication), quote = FALSE)