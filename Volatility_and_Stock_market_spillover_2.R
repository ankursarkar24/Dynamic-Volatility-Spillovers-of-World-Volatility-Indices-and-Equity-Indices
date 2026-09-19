# ==============================================================================
# 0. Load Libraries
# ==============================================================================
suppressPackageStartupMessages({
  library(quantmod)
  library(urca)
  library(moments)
  library(tseries)
  library(ConnectednessApproach)
  library(zoo)
  library(ggplot2)
  library(igraph)
})

# ==============================================================================
# 1. Data Ingestion & Alignment
# ==============================================================================
group1_tickers <- c("VIX" = "^VIX", "OVX" = "^OVX", "GVZ" = "^GVZ")
group2_tickers <- c("SP500"  = "^GSPC",
                    "FTSE"   = "^FTSE",
                    "DAX"    = "^GDAXI",
                    "NIKKEI" = "^N225",
                    "NIFTY"  = "^NSEI",
                    "HSI"    = "^HSI",
                    "TSX"    = "^GSPTSE")

all_tickers <- c(group1_tickers, group2_tickers)

start_date <- "2016-01-01"
end_date   <- "2026-06-30"

price_list <- list()

for (sym_name in names(all_tickers)) {
  ticker_sym <- all_tickers[sym_name]
  raw_obj    <- getSymbols(ticker_sym, src = "yahoo", 
                           from = start_date, to = end_date, auto.assign = FALSE)
  adj_close  <- na.locf(Ad(raw_obj), na.rm = TRUE)
  price_list[[sym_name]] <- adj_close
}

merged_raw     <- do.call(merge, price_list)
aligned_prices <- na.omit(na.locf(merged_raw))
colnames(aligned_prices) <- names(all_tickers)

# Group 1 (Volatility Indices): Preserve index levels
vol_g1 <- log(aligned_prices[-1, names(group1_tickers)])

# Group 2 (Stock Indices): Daily realized variance proxy (returns squared)
prices_g2   <- aligned_prices[, names(group2_tickers)]
returns_pct <- 100 * diff(log(prices_g2))[-1, ]
vol_g2      <- returns_pct^2

# Combine keeping date index intact
vol_xts  <- merge(vol_g1, vol_g2)
vol_zoo  <- as.zoo(vol_xts)
vol_df   <- as.data.frame(vol_xts)
real_dates <- as.Date(index(vol_xts))

cat("Any remaining missing values:", sum(is.na(vol_df)), "\n")
cat("Cleaned sample dimensions:", nrow(vol_df), "rows x", ncol(vol_df), "columns\n")

# ==============================================================================
# 2. Table 2: Descriptive Statistics
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
cat("\n=== Table 2: Descriptive Statistics ===\n")
print(round(table2, 3))

# ==============================================================================
# 3. Table 3: Correlation Matrix
# ==============================================================================
table3 <- cor(vol_df)
cat("\n=== Table 3: Correlation Matrix ===\n")
print(round(table3, 3))

# ==============================================================================
# 4. Table 4: Unit Root Tests (ADF & ERS Point Optimal)
# ==============================================================================
adf_models <- lapply(vol_df, function(x) ur.df(x, type = "drift", selectlags = "BIC"))
ers_models <- lapply(vol_df, function(x) ur.ers(x, type = "P-test", model = "constant"))

table4 <- data.frame(
  ADF_Test_Stat = sapply(adf_models, function(m) m@teststat[1]),
  ADF_5pct_Crit = sapply(adf_models, function(m) m@cval[1, "5pct"]),
  ERS_P_Stat    = sapply(ers_models, function(m) m@teststat[1]),
  ERS_5pct_Crit = sapply(ers_models, function(m) m@cval[1, "5pct"])
)
cat("\n=== Table 4: Dynamic Unit Root Test Results ===\n")
print(round(table4, 4))

# ==============================================================================
# 5. Dynamic Connectedness (TVP-VAR Baseline)
# ==============================================================================
dca <- ConnectednessApproach(
  vol_zoo,
  nlag          = 1,
  nfore         = 10,
  model         = "TVP-VAR",
  connectedness = "Time",
  corrected     = FALSE,
  VAR_config    = list(
    TVPVAR = list(
      kappa1 = 0.99,
      kappa2 = 0.99,
      prior  = "BayesPrior"
    )
  )
)

cat("\n=== Table 5: Dynamic Connectedness Matrix ===\n")
print(dca$TABLE)

# Sub-system groups
groups <- list(
  "Volatility_Indices" = c(1, 2, 3),
  "Stock_Markets"      = c(4:10)
)

dca_internal <- InternalConnectedness(dca, groups = groups)
dca_external <- ExternalConnectedness(dca, groups = groups)

cat("\n--- Connectedness Breakdown ---\n")
cat("Total System TCI       :", round(mean(dca$TCI, na.rm = TRUE), 2), "%\n")
cat("Internal TCI (Within)  :", round(mean(dca_internal$TCI, na.rm = TRUE), 2), "%\n")
cat("External TCI (Between) :", round(mean(dca_external$TCI, na.rm = TRUE), 2), "%\n")

# ==============================================================================
# 6. Figure 1: Dynamic Total, Internal, and Group Connectedness
# ==============================================================================
vol_zoo_g1 <- vol_zoo[, 1:3]
dca_g1 <- ConnectednessApproach(vol_zoo_g1, nlag = 1, nfore = 10, model = "TVP-VAR", 
                                connectedness = "Time", corrected = FALSE)

vol_zoo_g2 <- vol_zoo[, 4:10]
dca_g2 <- ConnectednessApproach(vol_zoo_g2, nlag = 1, nfore = 10, model = "TVP-VAR", 
                                connectedness = "Time", corrected = FALSE)

# Match estimation length
plot_dates <- tail(real_dates, length(dca$TCI))

df_fig1 <- data.frame(
  Date          = plot_dates,
  Total_TCI     = as.numeric(dca$TCI),
  Internal_TCI  = as.numeric(dca_internal$TCI),
  Vol_Indices   = as.numeric(dca_g1$TCI),
  Stock_Markets = as.numeric(dca_g2$TCI)
)

p_fig1 <- ggplot(df_fig1, aes(x = Date)) +
  geom_area(aes(y = Total_TCI, fill = "Total TCI"), alpha = 0.85) +
  geom_line(aes(y = Internal_TCI, color = "Internal TCI"), linewidth = 0.8) +
  geom_line(aes(y = Vol_Indices, color = "Volatility Indices"), linewidth = 0.8) +
  geom_line(aes(y = Stock_Markets, color = "Stock Markets"), linewidth = 0.8) +
  scale_fill_manual(name = "", values = c("Total TCI" = "grey40")) +
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
    legend.position        = "inside",
    legend.position.inside = c(0.22, 0.82),
    legend.background      = element_rect(fill = "transparent"),
    axis.line              = element_line(color = "black")
  )
print(p_fig1)

# ==============================================================================
# 7. Figure 2 & Table 6: Cross-Group Connectedness & Market Regimes
# ==============================================================================
df_fig2 <- data.frame(
  Date         = plot_dates,
  External_TCI = as.numeric(dca_external$TCI)
)

peak_val   <- max(df_fig2$External_TCI, na.rm = TRUE)
peak_date  <- df_fig2$Date[which.max(df_fig2$External_TCI)]
peak_label <- paste0("Peak: ", round(peak_val, 2), "% (", format(peak_date, "%d %b %Y"), ")")

p_fig2 <- ggplot(df_fig2, aes(x = Date, y = External_TCI)) +
  geom_line(color = "black", linewidth = 0.7) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = as.Date("2022-02-24"), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = as.Date("2022-03-16"), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = as.Date("2025-01-20"), linetype = "dashed", color = "grey40") +
  annotate("point", x = peak_date, y = peak_val, color = "red", size = 2) +
  annotate("text", x = peak_date, y = peak_val + (peak_val * 0.05), 
           label = peak_label, size = 3.5, fontface = "bold") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = "", y = "External TCI (%)", title = "Figure 2. Dynamic Cross-Group Connectedness") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
print(p_fig2)

calc_regime_stats <- function(df, start_d, end_d, regime_name) {
  sub_df <- subset(df, Date >= as.Date(start_d) & Date <= as.Date(end_d))
  if (nrow(sub_df) == 0) return(NULL)
  data.frame(
    Period            = regime_name,
    Mean_External_TCI = round(mean(sub_df$External_TCI, na.rm = TRUE), 2),
    Max_External_TCI  = round(max(sub_df$External_TCI, na.rm = TRUE), 2),
    Peak_Date         = format(sub_df$Date[which.max(sub_df$External_TCI)], "%d %B %Y")
  )
}

table6 <- rbind(
  calc_regime_stats(df_fig2, "2016-01-01", "2019-12-31", "Pre-COVID"),
  calc_regime_stats(df_fig2, "2020-01-01", "2021-12-31", "COVID-19"),
  calc_regime_stats(df_fig2, "2022-01-01", "2022-12-31", "Russia-Ukraine & Tightening"),
  calc_regime_stats(df_fig2, "2023-01-01", "2024-12-31", "Post-Crisis Expansion"),
  calc_regime_stats(df_fig2, "2025-01-01", "2026-06-30", "Macro & Tariff Uncertainty")
)
cat("\n=== Table 6: External Connectedness Across Selected Market Regimes ===\n")
print(table6, row.names = FALSE)

# ==============================================================================
# 8. Figure 3: Dynamic Net Directional Connectedness
# ==============================================================================
net_total    <- dca$NET
net_internal <- dca_internal$NET
var_names    <- colnames(net_total)

fig3_list <- lapply(var_names, function(v) {
  data.frame(
    Date     = plot_dates,
    Variable = v,
    Total    = as.numeric(net_total[, v]),
    Internal = as.numeric(net_internal[, v])
  )
})
df_fig3 <- do.call(rbind, fig3_list)
df_fig3$Variable <- factor(df_fig3$Variable, levels = var_names)

y_min <- min(c(df_fig3$Total, df_fig3$Internal), na.rm = TRUE)
y_max <- max(c(df_fig3$Total, df_fig3$Internal), na.rm = TRUE)
y_pad <- (y_max - y_min) * 0.05

p_fig3 <- ggplot(df_fig3, aes(x = Date)) +
  geom_area(aes(y = Total), fill = "black", alpha = 0.85) +
  geom_line(aes(y = Internal), color = "red", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.4) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(limits = c(y_min - y_pad, y_max + y_pad)) +
  facet_wrap(~ Variable, ncol = 3, scales = "fixed") +
  labs(x = "", y = "NET Connectedness", 
       title = "Figure 3. Dynamic Total and Within-Group Net Directional Connectedness") +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )
print(p_fig3)

# ==============================================================================
# 9. Figures 4, 5, 6: Pairwise NPDC Plots
# ==============================================================================
plot_npdc_figure <- function(dca_obj, base_var, target_vars, fig_title) {
  npdc_arr  <- dca_obj$NPDC
  all_names <- colnames(vol_zoo)
  base_idx  <- which(all_names == base_var)
  
  dims <- dim(npdc_arr)
  time_in_dim1 <- (dims[1] == length(plot_dates))
  
  plot_list <- lapply(target_vars, function(tgt) {
    tgt_idx <- which(all_names == tgt)
    val_series <- if (time_in_dim1) {
      as.numeric(npdc_arr[, base_idx, tgt_idx])
    } else {
      as.numeric(npdc_arr[base_idx, tgt_idx, ])
    }
    data.frame(
      Date = plot_dates,
      Pair = paste0(base_var, " - ", tgt),
      NPDC = val_series
    )
  })
  
  df_plot <- do.call(rbind, plot_list)
  df_plot$Pair <- factor(df_plot$Pair, levels = paste0(base_var, " - ", target_vars))
  
  y_min_5 <- floor(min(df_plot$NPDC, na.rm = TRUE) / 5) * 5
  y_max_5 <- ceiling(max(df_plot$NPDC, na.rm = TRUE) / 5) * 5
  
  ggplot(df_plot, aes(x = Date, y = NPDC)) +
    geom_area(fill = "black", alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(limits = c(y_min_5, y_max_5), breaks = seq(y_min_5, y_max_5, by = 5)) +
    facet_wrap(~ Pair, ncol = 2, scales = "fixed") +
    labs(x = "", y = "NPDC", title = fig_title) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1)
    )
}

stock_targets <- names(group2_tickers)

print(plot_npdc_figure(dca, "VIX", stock_targets, "Figure 4. NPDC: VIX to International Stock Markets"))
print(plot_npdc_figure(dca, "OVX", stock_targets, "Figure 5. NPDC: OVX to International Stock Markets"))
print(plot_npdc_figure(dca, "GVZ", stock_targets, "Figure 6. NPDC: GVZ to International Stock Markets"))

# ==============================================================================
# 10. Figure 7: Average NPDC Network Diagram
# ==============================================================================
dims <- dim(dca$NPDC)
time_in_dim1 <- (dims[1] == length(plot_dates))

if (time_in_dim1) {
  avg_npdc <- apply(dca$NPDC, c(2, 3), mean, na.rm = TRUE)
} else {
  avg_npdc <- apply(dca$NPDC, c(1, 2), mean, na.rm = TRUE)
}

colnames(avg_npdc) <- colnames(vol_zoo)
rownames(avg_npdc) <- colnames(vol_zoo)

net_adj <- avg_npdc
net_adj[net_adj < 0] <- 0

g <- graph_from_adjacency_matrix(net_adj, mode = "directed", weighted = TRUE, diag = FALSE)

avg_net <- colMeans(dca$NET, na.rm = TRUE)
avg_to  <- colMeans(dca$TO, na.rm = TRUE)

V(g)$color       <- ifelse(avg_net[V(g)$name] > 0, "#2B6CB0", "#E6AF2E")
V(g)$frame.color <- "grey40"
V(g)$label.color <- "black"
V(g)$label.font  <- 2
V(g)$label.cex   <- 0.85
V(g)$size        <- 12 + 18 * (avg_to[V(g)$name] - min(avg_to)) / (max(avg_to) - min(avg_to) + 1e-5)

edge_w     <- E(g)$weight
E(g)$width <- 0.5 + 4 * (edge_w - min(edge_w)) / (max(edge_w) - min(edge_w) + 1e-5)
E(g)$arrow.size  <- 0.35
E(g)$arrow.width <- 0.8
E(g)$color       <- rgb(0.3, 0.3, 0.3, 0.5)

par(mar = c(1, 1, 2, 1))
plot(g, layout = layout_in_circle(g), 
     main = "Figure 7. Average NPDC Network of Global Volatility Indices and International Stock Markets")
legend("bottomright", legend = c("Net Transmitter", "Net Receiver"), 
       fill = c("#2B6CB0", "#E6AF2E"), bty = "n", cex = 0.85)

# ==============================================================================
# 11. Table 7: Robustness Checks
# ==============================================================================
cat("\n--- Running Robustness Models (Table 7) ---\n")

dca_h5    <- ConnectednessApproach(vol_zoo, nlag = 1, nfore = 5,  model = "TVP-VAR", connectedness = "Time", corrected = FALSE)
dca_h20   <- ConnectednessApproach(vol_zoo, nlag = 1, nfore = 20, model = "TVP-VAR", connectedness = "Time", corrected = FALSE)
dca_p2    <- ConnectednessApproach(vol_zoo, nlag = 2, nfore = 10, model = "TVP-VAR", connectedness = "Time", corrected = FALSE)
dca_dy200 <- ConnectednessApproach(vol_zoo, nlag = 1, nfore = 10, window.size = 200, model = "VAR", connectedness = "Time", corrected = FALSE)
dca_dy250 <- ConnectednessApproach(vol_zoo, nlag = 1, nfore = 10, window.size = 250, model = "VAR", connectedness = "Time", corrected = FALSE)

format_table7_paper_style <- function(models_list, base_model) {
  col_names <- names(models_list)
  base_net  <- colMeans(base_model$NET, na.rm = TRUE)
  base_tx   <- names(sort(base_net[base_net > 0], decreasing = TRUE))
  base_rx   <- names(sort(base_net[base_net < 0], decreasing = FALSE))
  base_tci_zoo <- base_model$TCI
  
  final_mat <- matrix("", nrow = 7, ncol = length(models_list),
                      dimnames = list(
                        c("Average TCI (%)",
                          "Main net transmitters",
                          "Main net receivers",
                          "Transmitter/receiver role changes",
                          "Dynamic TCI correlation with baseline",
                          "Maximum TCI (%)",
                          "Peak date"),
                        col_names
                      ))
  
  for (i in seq_along(models_list)) {
    m_name <- col_names[i]
    m_obj  <- models_list[[i]]
    tci_z  <- m_obj$TCI
    
    final_mat[1, i] <- sprintf("%.2f", mean(tci_z, na.rm = TRUE))
    
    m_net  <- colMeans(m_obj$NET, na.rm = TRUE)
    cur_tx <- names(sort(m_net[m_net > 0], decreasing = TRUE))
    cur_rx <- names(sort(m_net[m_net < 0], decreasing = FALSE))
    
    if (i == 1) {
      final_mat[2, i] <- paste(cur_tx, collapse = ", ")
      final_mat[3, i] <- paste(cur_rx, collapse = ", ")
      final_mat[4, i] <- "–"
      final_mat[5, i] <- "1.0000"
      final_mat[6, i] <- sprintf("%.2f", max(tci_z, na.rm = TRUE))
      final_mat[7, i] <- format(as.Date(index(tci_z)[which.max(tci_z)]), "%d %B %Y")
    } else {
      final_mat[2, i] <- ifelse(setequal(cur_tx, base_tx), "Same", paste(cur_tx, collapse = ", "))
      final_mat[3, i] <- ifelse(setequal(cur_rx, base_rx), "Same", paste(cur_rx, collapse = ", "))
      
      sign_base <- sign(base_net)
      sign_cur  <- sign(m_net[names(base_net)])
      flips     <- names(which(sign_base != sign_cur))
      final_mat[4, i] <- if (length(flips) == 0) "0" else paste0(length(flips), " (", paste(flips, collapse = ", "), ")")
      
      if (m_name %in% c("H = 5", "H = 20")) {
        final_mat[5, i] <- "–"
        final_mat[6, i] <- "–"
        final_mat[7, i] <- "–"
      } else {
        # Safe calendar timestamp alignment
        merged_tci <- na.omit(merge.zoo(tci_z, base_tci_zoo))
        final_mat[5, i] <- sprintf("%.4f", cor(merged_tci[, 1], merged_tci[, 2]))
        final_mat[6, i] <- sprintf("%.2f", max(tci_z, na.rm = TRUE))
        final_mat[7, i] <- format(as.Date(index(tci_z)[which.max(tci_z)]), "%d %B %Y")
      }
    }
  }
  return(as.data.frame(final_mat))
}

robustness_models <- list(
  "Baseline TVP-VAR"   = dca,
  "H = 5"              = dca_h5,
  "H = 20"             = dca_h20,
  "Lag = 2"            = dca_p2,
  "DY VAR (200-Day)"   = dca_dy200,
  "DY VAR (250-Day)"   = dca_dy250
)

table7_publication <- format_table7_paper_style(robustness_models, dca)
cat("\n=== Table 7: Robustness Checks across Alternative Specifications ===\n")
print(as.matrix(table7_publication), quote = FALSE)














format_table7_paper_style <- function(models_list, base_model, date_vector) {
  col_names <- names(models_list)
  base_net  <- colMeans(base_model$NET, na.rm = TRUE)
  base_tx   <- names(sort(base_net[base_net > 0], decreasing = TRUE))
  base_rx   <- names(sort(base_net[base_net < 0], decreasing = FALSE))
  
  # Build a formal zoo object for the baseline TCI
  base_tci_raw <- as.numeric(base_model$TCI)
  base_dates   <- tail(date_vector, length(base_tci_raw))
  base_tci_zoo <- zoo(base_tci_raw, order.by = base_dates)
  
  final_mat <- matrix("", nrow = 7, ncol = length(models_list),
                      dimnames = list(
                        c("Average TCI (%)",
                          "Main net transmitters",
                          "Main net receivers",
                          "Transmitter/receiver role changes",
                          "Dynamic TCI correlation with baseline",
                          "Maximum TCI (%)",
                          "Peak date"),
                        col_names
                      ))
  
  for (i in seq_along(models_list)) {
    m_name  <- col_names[i]
    m_obj   <- models_list[[i]]
    tci_raw <- as.numeric(m_obj$TCI)
    
    # Dynamically align calendar dates for this specific model (e.g. rolling windows)
    m_dates <- tail(date_vector, length(tci_raw))
    tci_z   <- zoo(tci_raw, order.by = m_dates)
    
    final_mat[1, i] <- sprintf("%.2f", mean(tci_raw, na.rm = TRUE))
    
    m_net  <- colMeans(m_obj$NET, na.rm = TRUE)
    cur_tx <- names(sort(m_net[m_net > 0], decreasing = TRUE))
    cur_rx <- names(sort(m_net[m_net < 0], decreasing = FALSE))
    
    if (i == 1) {
      final_mat[2, i] <- paste(cur_tx, collapse = ", ")
      final_mat[3, i] <- paste(cur_rx, collapse = ", ")
      final_mat[4, i] <- "–"
      final_mat[5, i] <- "1.0000"
      final_mat[6, i] <- sprintf("%.2f", max(tci_raw, na.rm = TRUE))
      final_mat[7, i] <- format(m_dates[which.max(tci_raw)], "%d %B %Y")
    } else {
      final_mat[2, i] <- ifelse(setequal(cur_tx, base_tx), "Same", paste(cur_tx, collapse = ", "))
      final_mat[3, i] <- ifelse(setequal(cur_rx, base_rx), "Same", paste(cur_rx, collapse = ", "))
      
      sign_base <- sign(base_net)
      sign_cur  <- sign(m_net[names(base_net)])
      flips     <- names(which(sign_base != sign_cur))
      final_mat[4, i] <- if (length(flips) == 0) "0" else paste0(length(flips), " (", paste(flips, collapse = ", "), ")")
      
      if (m_name %in% c("H = 5", "H = 20")) {
        final_mat[5, i] <- "–"
        final_mat[6, i] <- "–"
        final_mat[7, i] <- "–"
      } else {
        # Both objects are now validated zoo time series
        merged_tci      <- na.omit(merge.zoo(tci_z, base_tci_zoo, all = FALSE))
        final_mat[5, i] <- sprintf("%.4f", cor(as.numeric(merged_tci[, 1]), as.numeric(merged_tci[, 2])))
        final_mat[6, i] <- sprintf("%.2f", max(tci_raw, na.rm = TRUE))
        final_mat[7, i] <- format(m_dates[which.max(tci_raw)], "%d %B %Y")
      }
    }
  }
  return(as.data.frame(final_mat))
}

table7_publication <- format_table7_paper_style(robustness_models, dca, real_dates)
cat("\n=== Table 7: Robustness Checks across Alternative Specifications ===\n")
print(as.matrix(table7_publication), quote = FALSE)





















format_table7_paper_style <- function(models_list, base_model, date_vector) {
  col_names <- names(models_list)
  
  # 1. Baseline metrics
  base_net  <- colMeans(base_model$NET, na.rm = TRUE)
  base_tx   <- names(sort(base_net[base_net > 0], decreasing = TRUE))
  base_rx   <- names(sort(base_net[base_net < 0], decreasing = FALSE))
  
  # Build baseline lookup dataframe: Date + TCI
  base_tci_raw <- as.numeric(base_model$TCI)
  base_dates   <- as.Date(tail(date_vector, length(base_tci_raw)))
  df_base      <- data.frame(Date = base_dates, Base_TCI = base_tci_raw)
  
  # 2. Result Matrix Template
  final_mat <- matrix("", nrow = 7, ncol = length(models_list),
                      dimnames = list(
                        c("Average TCI (%)",
                          "Main net transmitters",
                          "Main net receivers",
                          "Transmitter/receiver role changes",
                          "Dynamic TCI correlation with baseline",
                          "Maximum TCI (%)",
                          "Peak date"),
                        col_names
                      ))
  
  # 3. Process each specification
  for (i in seq_along(models_list)) {
    m_name  <- col_names[i]
    m_obj   <- models_list[[i]]
    tci_raw <- as.numeric(m_obj$TCI)
    
    # Extract matching calendar dates for this specific model
    m_dates <- as.Date(tail(date_vector, length(tci_raw)))
    
    # Row 1: Average TCI
    final_mat[1, i] <- sprintf("%.2f", mean(tci_raw, na.rm = TRUE))
    
    # Net Directional Rankings
    m_net  <- colMeans(m_obj$NET, na.rm = TRUE)
    cur_tx <- names(sort(m_net[m_net > 0], decreasing = TRUE))
    cur_rx <- names(sort(m_net[m_net < 0], decreasing = FALSE))
    
    if (i == 1) {
      final_mat[2, i] <- paste(cur_tx, collapse = ", ")
      final_mat[3, i] <- paste(cur_rx, collapse = ", ")
      final_mat[4, i] <- "–"
      final_mat[5, i] <- "1.0000"
      final_mat[6, i] <- sprintf("%.2f", max(tci_raw, na.rm = TRUE))
      final_mat[7, i] <- format(m_dates[which.max(tci_raw)], "%d %B %Y")
    } else {
      # Rows 2 & 3: Transmitters & Receivers
      final_mat[2, i] <- ifelse(setequal(cur_tx, base_tx), "Same", paste(cur_tx, collapse = ", "))
      final_mat[3, i] <- ifelse(setequal(cur_rx, base_rx), "Same", paste(cur_rx, collapse = ", "))
      
      # Row 4: Role Changes (sign flips)
      sign_base <- sign(base_net)
      sign_cur  <- sign(m_net[names(base_net)])
      flips     <- names(which(sign_base != sign_cur))
      final_mat[4, i] <- if (length(flips) == 0) "0" else paste0(length(flips), " (", paste(flips, collapse = ", "), ")")
      
      # Rows 5 to 7: Correlation, Max, and Peak Date
      if (m_name %in% c("H = 5", "H = 20")) {
        final_mat[5, i] <- "–"
        final_mat[6, i] <- "–"
        final_mat[7, i] <- "–"
      } else {
        # Safe base R alignment via merge() on Date column
        df_cur <- data.frame(Date = m_dates, Model_TCI = tci_raw)
        merged_df <- merge(df_cur, df_base, by = "Date")
        
        cor_val <- cor(merged_df$Model_TCI, merged_df$Base_TCI, use = "complete.obs")
        final_mat[5, i] <- sprintf("%.4f", cor_val)
        final_mat[6, i] <- sprintf("%.2f", max(tci_raw, na.rm = TRUE))
        final_mat[7, i] <- format(m_dates[which.max(tci_raw)], "%d %B %Y")
      }
    }
  }
  
  return(as.data.frame(final_mat))
}

table7_publication <- format_table7_paper_style(robustness_models, dca, real_dates)

cat("\n=== Table 7: Robustness Checks across Alternative Specifications ===\n")
print(as.matrix(table7_publication), quote = FALSE)


























# ==============================================================================
# SECTION 12: PUBLICATION-GRADE VISUALIZATIONS
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(grid)
  library(igraph)
})

# ------------------------------------------------------------------------------
# Universal Publication Theme & Palette
# ------------------------------------------------------------------------------
theme_pub <- function(base_size = 11, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title         = element_text(face = "bold", size = rel(1.15), color = "#1A202C", margin = margin(b = 6)),
      plot.subtitle      = element_text(size = rel(0.92), color = "#4A5568", margin = margin(b = 14)),
      plot.caption       = element_text(size = rel(0.75), color = "#718096", hjust = 0, margin = margin(t = 10)),
      panel.grid.major.y = element_line(color = "#E2E8F0", linewidth = 0.4),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.line.x        = element_line(color = "#4A5568", linewidth = 0.5),
      axis.ticks.x       = element_line(color = "#4A5568", linewidth = 0.4),
      axis.ticks.length  = unit(3, "pt"),
      axis.text          = element_text(color = "#2D3748", size = rel(0.85)),
      axis.title         = element_text(face = "bold", color = "#1A202C", size = rel(0.9)),
      strip.background   = element_rect(fill = "#F7FAFC", color = "#CBD5E0", linewidth = 0.6),
      strip.text         = element_text(face = "bold", color = "#2D3748", size = rel(0.88)),
      plot.margin        = margin(t = 12, r = 14, b = 12, l = 12)
    )
}

event_markers <- list(
  data.frame(date = as.Date("2020-03-11"), label = "COVID-19"),
  data.frame(date = as.Date("2022-02-24"), label = "Ukraine"),
  data.frame(date = as.Date("2022-03-16"), label = "Fed Hike")
)

# ------------------------------------------------------------------------------
# Figure 1: Total, Internal, and Group Connectedness
# ------------------------------------------------------------------------------
p1_pub <- ggplot(df_fig1, aes(x = Date)) +
  geom_area(aes(y = Total_TCI, fill = "Total System TCI"), alpha = 0.22) +
  geom_line(aes(y = Total_TCI, color = "Total System TCI"), linewidth = 1.05) +
  geom_line(aes(y = Internal_TCI, color = "Internal (Within-Group)"), linewidth = 0.85, linetype = "solid") +
  geom_line(aes(y = Vol_Indices, color = "Volatility Indices"), linewidth = 0.75, linetype = "twodash") +
  geom_line(aes(y = Stock_Markets, color = "Stock Markets"), linewidth = 0.75, linetype = "longdash") +
  scale_fill_manual(name = "", values = c("Total System TCI" = "#1A365D")) +
  scale_color_manual(
    name = "",
    values = c(
      "Total System TCI"        = "#1A365D",
      "Internal (Within-Group)" = "#C53030",
      "Volatility Indices"      = "#2F855A",
      "Stock Markets"           = "#DD6B20"
    )
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
  labs(
    title    = "Figure 1. Dynamic Total, Internal, and Group-Specific Volatility Connectedness",
    subtitle = "Time-Varying Parameter VAR estimation across global volatility indices and international equities (2016–2026)",
    x        = NULL,
    y        = "TCI (%)",
    caption  = "Note: Total TCI represents generalized forecast error variance decomposition (H = 10, TVP-VAR kappa1 = kappa2 = 0.99)."
  ) +
  theme_pub() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.24, 0.82),
    legend.background      = element_rect(fill = alpha("white", 0.85), color = "#CBD5E0", linewidth = 0.4),
    legend.key.width       = unit(20, "pt"),
    legend.text            = element_text(size = rel(0.85), face = "bold")
  )

print(p1_pub)

# ------------------------------------------------------------------------------
# Figure 2: Cross-Group Connectedness with Event Overlays
# ------------------------------------------------------------------------------
p2_pub <- ggplot(df_fig2, aes(x = Date, y = External_TCI)) +
  geom_area(fill = "#2B6CB0", alpha = 0.12) +
  geom_line(color = "#2B6CB0", linewidth = 0.95) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype = "dashed", color = "#A0AEC0", linewidth = 0.5) +
  geom_vline(xintercept = as.Date("2022-02-24"), linetype = "dashed", color = "#A0AEC0", linewidth = 0.5) +
  geom_vline(xintercept = as.Date("2022-03-16"), linetype = "dashed", color = "#A0AEC0", linewidth = 0.5) +
  geom_vline(xintercept = as.Date("2025-01-20"), linetype = "dashed", color = "#A0AEC0", linewidth = 0.5) +
  annotate("text", x = as.Date("2020-03-11"), y = peak_val * 0.97, label = "COVID-19", 
           angle = 90, vjust = -0.6, hjust = 1, size = 2.9, color = "#4A5568", fontface = "italic") +
  annotate("text", x = as.Date("2022-02-24"), y = peak_val * 0.97, label = "Ukraine War", 
           angle = 90, vjust = -0.6, hjust = 1, size = 2.9, color = "#4A5568", fontface = "italic") +
  annotate("point", x = peak_date, y = peak_val, color = "#C53030", size = 2.5) +
  annotate("label", x = peak_date, y = peak_val + (peak_val * 0.04), 
           label = paste0("Global Peak: ", round(peak_val, 2), "%\n", format(peak_date, "%b %Y")),
           size = 3.1, fontface = "bold", color = "#742A2A", fill = "#FFF5F5", label.size = 0.3) +
  scale_y_continuous(limits = c(0, peak_val * 1.12), expand = c(0, 0)) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
  labs(
    title    = "Figure 2. Cross-System Information Flow: Volatility to Stock Markets",
    subtitle = "External Total Connectedness Index between implied volatility benchmarks and equity index variances",
    x        = NULL,
    y        = "External TCI (%)"
  ) +
  theme_pub()

print(p2_pub)

# ------------------------------------------------------------------------------
# Figure 3: Dynamic Total & Internal Net Directional Spillovers
# ------------------------------------------------------------------------------
p3_pub <- ggplot(df_fig3, aes(x = Date)) +
  geom_hline(yintercept = 0, color = "#4A5568", linewidth = 0.4) +
  geom_area(aes(y = Total, fill = Total > 0), alpha = 0.4) +
  geom_line(aes(y = Total, color = Total > 0), linewidth = 0.6) +
  geom_line(aes(y = Internal), color = "#C53030", linewidth = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = "#2B6CB0", "FALSE" = "#DD6B20"), guide = "none") +
  scale_color_manual(values = c("TRUE" = "#1A365D", "FALSE" = "#9C4221"), guide = "none") +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  scale_y_continuous(limits = c(y_min - y_pad, y_max + y_pad)) +
  facet_wrap(~ Variable, ncol = 3, scales = "fixed") +
  labs(
    title    = "Figure 3. Dynamic Total and Within-Group Net Spillovers",
    subtitle = "Blue = Net Volatility Transmitter (>0) | Orange = Net Receiver (<0) | Red Dashed = Internal Spillover",
    x        = NULL,
    y        = "Net Spillover (%)"
  ) +
  theme_pub(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p3_pub)

# ------------------------------------------------------------------------------
# Figures 4, 5, 6: Pairwise NPDC Plots
# ------------------------------------------------------------------------------
plot_npdc_pub <- function(dca_obj, base_var, target_vars, main_title) {
  npdc_arr     <- dca_obj$NPDC
  all_names    <- colnames(vol_zoo)
  base_idx     <- which(all_names == base_var)
  time_in_dim1 <- (dim(npdc_arr)[1] == length(plot_dates))
  
  plot_list <- lapply(target_vars, function(tgt) {
    tgt_idx <- which(all_names == tgt)
    val_series <- if (time_in_dim1) {
      as.numeric(npdc_arr[, base_idx, tgt_idx])
    } else {
      as.numeric(npdc_arr[base_idx, tgt_idx, ])
    }
    data.frame(
      Date = plot_dates,
      Pair = paste0(base_var, " \u2192 ", tgt),
      NPDC = val_series
    )
  })
  
  df_plot <- do.call(rbind, plot_list)
  df_plot$Pair <- factor(df_plot$Pair, levels = paste0(base_var, " \u2192 ", target_vars))
  
  y_min_val <- floor(min(df_plot$NPDC, na.rm = TRUE) / 5) * 5
  y_max_val <- ceiling(max(df_plot$NPDC, na.rm = TRUE) / 5) * 5
  
  ggplot(df_plot, aes(x = Date, y = NPDC)) +
    geom_hline(yintercept = 0, color = "#4A5568", linewidth = 0.4) +
    geom_area(aes(fill = NPDC > 0), alpha = 0.45) +
    geom_line(aes(color = NPDC > 0), linewidth = 0.55) +
    scale_fill_manual(values = c("TRUE" = "#2B6CB0", "FALSE" = "#DD6B20"), guide = "none") +
    scale_color_manual(values = c("TRUE" = "#1A365D", "FALSE" = "#9C4221"), guide = "none") +
    scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
    scale_y_continuous(limits = c(y_min_val, y_max_val), breaks = seq(y_min_val, y_max_val, by = 5)) +
    facet_wrap(~ Pair, ncol = 2, scales = "fixed") +
    labs(
      title    = main_title,
      subtitle = paste0("Positive values indicate ", base_var, " acts as a net transmitter of shocks to target equity"),
      x        = NULL,
      y        = "NPDC (%)"
    ) +
    theme_pub(base_size = 9.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

print(plot_npdc_pub(dca, "VIX", stock_targets, "Figure 4. Net Pairwise Spillovers: Equity Volatility (VIX) to World Markets"))
print(plot_npdc_pub(dca, "OVX", stock_targets, "Figure 5. Net Pairwise Spillovers: Crude Oil Volatility (OVX) to World Markets"))
print(plot_npdc_pub(dca, "GVZ", stock_targets, "Figure 6. Net Pairwise Spillovers: Gold Volatility (GVZ) to World Markets"))

# ------------------------------------------------------------------------------
# Figure 7: Publication Quality igraph Network Visualization
# ------------------------------------------------------------------------------
dims <- dim(dca$NPDC)
time_in_dim1 <- (dims[1] == length(plot_dates))

if (time_in_dim1) {
  avg_npdc <- apply(dca$NPDC, c(2, 3), mean, na.rm = TRUE)
} else {
  avg_npdc <- apply(dca$NPDC, c(1, 2), mean, na.rm = TRUE)
}

colnames(avg_npdc) <- colnames(vol_zoo)
rownames(avg_npdc) <- colnames(vol_zoo)

net_adj <- avg_npdc
net_adj[net_adj < 0] <- 0

g <- graph_from_adjacency_matrix(net_adj, mode = "directed", weighted = TRUE, diag = FALSE)

avg_net <- colMeans(dca$NET, na.rm = TRUE)
avg_to  <- colMeans(dca$TO, na.rm = TRUE)

# Vertex aesthetics
V(g)$is_tx       <- avg_net[V(g)$name] > 0
V(g)$color       <- ifelse(V(g)$is_tx, "#2B6CB0", "#E2E8F0")
V(g)$frame.color <- ifelse(V(g)$is_tx, "#1A365D", "#718096")
V(g)$frame.width <- 1.8
V(g)$label.color <- ifelse(V(g)$is_tx, "white", "#1A202C")
V(g)$label.font  <- 2
V(g)$label.cex   <- 0.85
V(g)$size        <- 16 + 18 * (avg_to[V(g)$name] - min(avg_to)) / (max(avg_to) - min(avg_to) + 1e-5)

# Edge aesthetics
edge_weights     <- E(g)$weight
E(g)$width       <- 0.6 + 4.5 * (edge_weights - min(edge_weights)) / (max(edge_weights) - min(edge_weights) + 1e-5)
E(g)$arrow.size  <- 0.35
E(g)$arrow.width <- 0.9
E(g)$color       <- rgb(0.2, 0.25, 0.3, 0.45)
E(g)$curved      <- 0.15

# Export rendering parameters
old_par <- par(no.readonly = TRUE)
par(mar = c(1, 1, 3, 1), bg = "white")

plot(
  g,
  layout = layout_in_circle(g),
  main   = "Figure 7. Average Net Directional Pairwise Network (2016–2026)"
)

legend(
  "bottomleft",
  legend = c("Net Transmitter", "Net Receiver"),
  pch    = 21,
  pt.bg  = c("#2B6CB0", "#E2E8F0"),
  col    = c("#1A365D", "#718096"),
  pt.cex = 1.8,
  bty    = "n",
  cex    = 0.9,
  text.font = 2
)
par(old_par)





