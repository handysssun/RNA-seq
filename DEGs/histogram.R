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
  geom_histogram(binwidth = 0.5, color = "black", alpha = 0.7) +
  scale_fill_identity() +  # 使用定义的颜色
  labs(
    title = "Distribution of logfc",
    x = "logfc",
    y = "Count"
  ) +
  theme_minimal()