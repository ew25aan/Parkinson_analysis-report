
# Edit csv_path before running
csv_path <- "C:/Users/admin/Downloads/archive/Parkinsson disease.csv"

# load required package
if (!requireNamespace("car", quietly = TRUE)) install.packages("car", repos = "https://cloud.r-project.org")
library(car)

# read CSV
if (!file.exists(csv_path)) stop("CSV file not found: ", csv_path)
df <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE),
               error = function(e) stop("Error reading CSV: ", e$message))

cat("Loaded:", nrow(df), "rows x", ncol(df), "cols\n")

# get the names and locate columns
names(df) <- make.names(names(df))
mdvp_col <- grep("^MDVP.*Fo.*Hz", names(df), ignore.case = TRUE, value = TRUE)
if (length(mdvp_col) == 0) stop("Cannot find MDVP:Fo(Hz) column")
mdvp_col <- mdvp_col[1]

status_col <- if ("status" %in% names(df)) "status" else {
  poss <- grep("status|class|target", names(df), ignore.case = TRUE, value = TRUE)
  if (length(poss) == 0) stop("Cannot find a binary status column") else poss[1]
}

cat("Using columns:", mdvp_col, "and", status_col, "\n")

# prepare groups
df[[status_col]] <- as.factor(df[[status_col]])
if (all(levels(df[[status_col]]) %in% c("0", "1"))) {
  df[[status_col]] <- factor(df[[status_col]], levels = c("0", "1"), labels = c("Healthy", "Parkinson"))
} else if (length(levels(df[[status_col]])) != 2) {
  stop("Status column must be binary (two levels).")
}

g1 <- levels(df[[status_col]])[1]
g2 <- levels(df[[status_col]])[2]
x1 <- as.numeric(df[df[[status_col]] == g1, mdvp_col])
x2 <- as.numeric(df[df[[status_col]] == g2, mdvp_col])

# summary stats
summary_stats <- data.frame(
  Group = c(g1, g2),
  N = c(sum(!is.na(x1)), sum(!is.na(x2))),
  Mean = c(mean(x1, na.rm = TRUE), mean(x2, na.rm = TRUE)),
  SD = c(sd(x1, na.rm = TRUE), sd(x2, na.rm = TRUE)),
  Median = c(median(x1, na.rm = TRUE), median(x2, na.rm = TRUE)),
  Min = c(min(x1, na.rm = TRUE), min(x2, na.rm = TRUE)),
  Max = c(max(x1, na.rm = TRUE), max(x2, na.rm = TRUE))
)
print(summary_stats)
write.csv(summary_stats, "summary_stats_mdvp_fo.csv", row.names = FALSE)

# plots
png("boxplot_mdvp_fo.png", width = 1200, height = 800, res = 150)
boxplot(df[[mdvp_col]] ~ df[[status_col]], xlab = "Group", ylab = "MDVP:Fo(Hz)",
        main = "Boxplot of MDVP:Fo(Hz) by Group", col = c("lightblue", "salmon"), outline = TRUE)
dev.off()

png("hist_mdvp_fo_groups.png", width = 1200, height = 800, res = 150)
h1 <- hist(x1, breaks = 20, plot = FALSE)
h2 <- hist(x2, breaks = 20, plot = FALSE)
ymax <- max(h1$counts, h2$counts) * 1.2
plot(h1, col = rgb(0.9,0.6,0.0,0.5), xlim = range(c(x1, x2), na.rm = TRUE),
     xlab = "MDVP:Fo(Hz)", main = "Histogram of MDVP:Fo(Hz) by Group", ylim = c(0, ymax))
plot(h2, col = rgb(0,0.6,1,0.4), add = TRUE)
legend("topright", legend = c(g1, g2), fill = c(rgb(0.9,0.6,0.0,0.5), rgb(0,0.6,1,0.4)))
dev.off()

# assumption checks
sh_1 <- if (length(x1) >= 3 && length(x1) <= 5000) tryCatch(shapiro.test(x1), error = function(e) e) else NA
sh_2 <- if (length(x2) >= 3 && length(x2) <= 5000) tryCatch(shapiro.test(x2), error = function(e) e) else NA

levene_res <- tryCatch(leveneTest(df[[mdvp_col]] ~ df[[status_col]]), error = function(e) e)

# decide t-test type
use_welch <- TRUE
if (is.data.frame(levene_res)) {
  pval <- levene_res$`Pr(>F)`[1]
  if (!is.na(pval) && pval >= 0.05) use_welch <- FALSE
}

t_res <- if (use_welch) {
  t.test(df[[mdvp_col]] ~ df[[status_col]], var.equal = FALSE)
} else {
  t.test(df[[mdvp_col]] ~ df[[status_col]], var.equal = TRUE)
}

# print and save outputs
cat("\nShapiro-Wilk (group 1):\n"); print(sh_1)
cat("\nShapiro-Wilk (group 2):\n"); print(sh_2)
cat("\nLevene's test:\n"); print(levene_res)
cat("\nT-test:\n"); print(t_res)

sink("statistical_tests_output.txt")
cat("Summary statistics:\n"); print(summary_stats)
cat("\nShapiro-Wilk tests:\n"); if (!identical(sh_1, NA)) print(sh_1); if (!identical(sh_2, NA)) print(sh_2)
cat("\nLevene's test:\n"); print(levene_res)
cat("\nT-test:\n"); print(t_res)
sink()

cat("Outputs saved:\n - summary_stats_mdvp_fo.csv\n - boxplot_mdvp_fo.png\n - hist_mdvp_fo_groups.png\n - statistical_tests_output.txt\n")
