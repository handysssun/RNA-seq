# 示例数据框
set.seed(123)
df <- data.frame(
  logfc = rnorm(1000, mean = 0, sd = 2),  # 随机生成 logfc 数据
  pvalue = runif(1000, min = 0, max = 1)  # 随机生成 pvalue 数据
)

# 定义颜色条件
df$color <- ifelse(df$logfc > 1 & df$pvalue < 0.05, "red",
                   ifelse(df$logfc < -1 & df$pvalue < 0.05, "blue", "gray"))

# 绘制直方图
library(ggplot2)
ggplot(df, aes(x = logfc, fill = color)) +
  geom_histogram(
    binwidth = 0.5,  # 设置柱形宽度为 0.5
    breaks = seq(floor(min(df$logfc)), ceiling(max(df$logfc)), by = 0.5),  # 设置分隔点
    color = "black", 
    alpha = 0.7
  ) +
  scale_fill_identity() +  # 使用定义的颜色
  labs(
    title = "Distribution of logfc (0.5 binwidth, centered at 0)",
    x = "logfc",
    y = "Count"
  ) +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")  # 添加 0 分隔线
