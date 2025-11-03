# ---- cor_plot.R --------------------------------------------------------------
# 依赖：tidyverse（ggplot2/dplyr/readr）、tools（判断扩展名）
# 若尚未安装：install.packages(c("tidyverse"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(tools)
})

# 读取任意常见表格
read_expression_table <- function(path) {
  ext <- tolower(file_ext(path))
  if (ext %in% c("csv")) {
    readr::read_csv(path, show_col_types = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("请先安装 readxl 包：install.packages('readxl')")
    }
    readxl::read_excel(path)
  } else {
    stop("不支持的文件格式：", ext)
  }
}

# 主函数：计算相关并绘图
# data_or_path: data.frame/tibble 或 文件路径（csv/tsv/xlsx）
# mir_col, mrna_col: 列名（miRNA 与 mRNA）
# method: "pearson" 或 "spearman"
# log_transform: 是否对数转换（log2(x+1)）
# out_png: 可选，若提供文件名则保存 PNG
cor_scatter_plot <- function(
    data_or_path,
    mir_col,
    mrna_col,
    method = c("pearson", "spearman"),
    log_transform = FALSE,
    out_png = NULL,
    point_alpha = 0.8,
    point_size = 2.2
) {
  method <- match.arg(method)
  
  # 读取数据
  df <- if (is.character(data_or_path)) read_expression_table(data_or_path) else as_tibble(data_or_path)
  
  # 基本检查
  if (!all(c(mir_col, mrna_col) %in% colnames(df))) {
    missing_cols <- setdiff(c(mir_col, mrna_col), colnames(df))
    stop("以下列在数据中未找到：", paste(missing_cols, collapse = ", "))
  }
  
  # 选取两列并清洗
  df2 <- df %>%
    dplyr::select(miRNA = all_of(mir_col), mRNA = all_of(mrna_col)) %>%
    mutate(
      miRNA = as.numeric(miRNA),
      mRNA  = as.numeric(mRNA)
    ) %>%
    filter(is.finite(miRNA), is.finite(mRNA))
  
  if (nrow(df2) < 3) stop("可用样本数 < 3，无法进行稳健的相关分析。")
  
  # 可选对数转换（常见于计数或强偏态表达数据）
  if (log_transform) {
    df2 <- df2 %>%
      mutate(
        miRNA = log2(miRNA + 1),
        mRNA  = log2(mRNA  + 1)
      )
  }
  
  # 相关与显著性检验
  ct <- suppressWarnings(cor.test(df2$miRNA, df2$mRNA, method = method, exact = FALSE))
  r  <- unname(ct$estimate)
  p  <- ct$p.value
  
  # 回归拟合（用于绘制线）
  fit <- lm(mRNA ~ miRNA, data = df2)
  
  # 构建标注文本（简洁显示）
  lab <- sprintf("%s r = %.3f\np = %s",
                 ifelse(method == "pearson", "Pearson", "Spearman"),
                 r,
                 ifelse(p < 2.2e-16, "<2.2e-16", formatC(p, format = "e", digits = 2)))
  
  # 绘图
  pplt <- ggplot(df2, aes(x = miRNA, y = mRNA)) +
    geom_point(alpha = point_alpha, size = point_size) +
    geom_smooth(method = "lm", se = TRUE) +
    annotate("text",
             x = Inf, y = -Inf, hjust = 1.05, vjust = -0.8,
             label = lab, size = 4) +
    labs(
      title = sprintf("Expression correlation: %s vs %s", mir_col, mrna_col),
      x = sprintf("%s expression%s", mir_col, ifelse(log_transform, " (log2+1)", "")),
      y = sprintf("%s expression%s", mrna_col, ifelse(log_transform, " (log2+1)", "")),
      caption = sprintf("n = %d | method = %s", nrow(df2), method)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(out_png)) {
    ggsave(filename = out_png, plot = pplt, width = 6, height = 5, dpi = 300)
  }
  
  list(
    r = r,
    p_value = p,
    method = method,
    n = nrow(df2),
    lm_coef = coef(fit),
    plot = pplt
  )
}

# -------------------------- 使用示例 ------------------------------------------
# 假设你的表达矩阵是 samples x genes 的宽表：
# | Sample | hsa-miR-21-5p | AKT1 | TP53 | ...
# 并保存在 "expr.csv"
#
res <- cor_scatter_plot(
  data_or_path = as_tibble(as.data.frame(t(sub_matrix))),
  mrna_col = "GGCT",
  mir_col = "hsa-mir-769",
  method = "pearson",       # "pearson" 或 "spearman"
  log_transform = FALSE,    # 若是count数据可设为 TRUE
  out_png = "miR769_GGCT_correlation_p.png"
)
print(res$plot)
res$r; res$p_value; res$lm_coef
# 用循环函数直接绘制多组----------------------------------------------------------------------------- 
# 'GGCT','hsa-mir-29a','hsa-mir-29b-1','hsa-mir-29b-2','hsa-mir-769'
for(i in c('hsa-mir-29a','hsa-mir-29b-1','hsa-mir-29b-2','hsa-mir-769')){
  res <- cor_scatter_plot(
    data_or_path = as_tibble(as.data.frame(t(sub_matrix))),
    mrna_col = "GGCT",
    mir_col = i,
    method = "pearson",       # "pearson" 或 "spearman"
    log_transform = FALSE,    # 若是count数据可设为 TRUE
    out_png = paste0(i,"_GGCT_correlation_p.png")
  )
}





















