library(iCAMP)
library(ape)
library(ggplot2)

# 2. 设置文件路径和参数
# 请替换为你的实际文件路径
wd <- "C:/Users/DELL/Desktop/iCAMP_test/"  # 输入文件的文件夹路径
save.wd <- "C:/Users/DELL/Desktop/iCAMP_test/"  # 输出文件的保存路径

# 文件名
com.file <- "relative_abundance.tsv"  # OTU 表（相对丰度数据）
tree.file <- "rep_seqs.nwk"  # 系统发育树文件

# 创建输出文件夹（如果不存在）
if (!dir.exists(save.wd)) dir.create(save.wd)

# 设置关键参数
prefix <- "iCAMP_analysis"  # 输出文件名前缀
rand.time <- 1000  # 随机化次数（测试时可设为 100 以节省时间）
nworker <- 4  # 并行计算线程数（根据 CPU 核心数调整，例如 4 或 8）
bin.size.limit <- 24  # 最小 bin 大小（根据 OTU 数量调整）

# 3. 读取数据
# 设置工作目录
setwd(wd)

# 读取 OTU 表（相对丰度数据）
comm <- read.table(com.file, header = TRUE, sep = "\t", row.names = 1, 
                   as.is = TRUE, stringsAsFactors = FALSE, comment.char = "")

# 转置 OTU 表（iCAMP 要求行为 OTU，列为样本）
comm <- t(comm)

# 检查 OTU 表是否为相对丰度（每列和应接近 1）
colSums(comm)  # 如果和接近 1，则数据已经标准化

# 读取系统发育树
tree <- read.tree(tree.file)

# 检查数据
head(comm)  # 查看 OTU 表前几行
plot(tree)  # 可视化系统发育树

# 4. 确保 OTU 名称一致
# 提取 OTU 表和树中的共有 OTU
common_otus <- intersect(colnames(comm), tree$tip.label)
if (length(common_otus) == 0) {
  stop("OTU 名称不匹配，请检查 comm 的列名和 tree$tip.label 是否一致")
}

# 筛选 OTU 表和树
comm <- comm[, common_otus]
tree <- drop.tip(tree, setdiff(tree$tip.label, common_otus))

# 5. 运行 iCAMP 分析
icamp.out <- icamp.big(
  comm = comm,  # OTU 表
  tree = tree,  # 系统发育树
  pd.wd = save.wd,  # 系统发育距离矩阵保存路径
  rand = rand.time,  # 随机化次数
  prefix = prefix,  # 输出文件名前缀
  nworker = nworker,  # 并行线程数
  bin.size.limit = bin.size.limit,  # 最小 bin 大小
  sp.check = TRUE,  # 检查物种名称一致性
  phylo.rand.scale = "within.bin",  # 系统发育随机化方法
  taxa.rand.scale = "across.all",  # 分类学随机化方法
  phylo.metric = "bMNTD",  # 系统发育指标（默认 bMNTD）
  sig.index = "SES.RC",  # 显著性检验方法
  detail.save = TRUE,  # 保存详细结果
  qp.save = TRUE,  # 保存量化过程结果
  ignore.zero = TRUE,  # 忽略零丰度
  unit.sum = NULL  # 数据为相对丰度，设置为 NULL
)

# 6. 保存结果
options(max.print = 999999)
sink(file.path(save.wd, "icamp_results.txt"))
print(icamp.out)
sink()

# 7. 可视化生态过程比例（可选）
# 提取生态过程结果
proc <- icamp.out$Pt  # 生态过程比例（Pt 包含样本对之间的过程比例）
proc_df <- as.data.frame(proc)

# 样本对名称可能需要手动提取，这里假设 proc_df 包含样本对和过程比例
# 如果 proc_df 的格式不正确，可检查 icamp.out$Pt 的结构
# 假设 proc_df 包含 "SamplePair"、"Process"、"Value" 列
# 如果需要调整，请根据实际输出修改
proc_df$SamplePair <- rownames(proc_df)  # 样本对名称
proc_df <- reshape2::melt(proc_df, id.vars = "SamplePair", 
                          variable.name = "Process", value.name = "Value")

# 绘制柱状图
ggplot(proc_df, aes(x = SamplePair, y = Value, fill = Process)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Sample Pair", y = "Relative Importance", fill = "Ecological Process") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Relative Importance of Ecological Processes")

# 保存图像
ggsave(file.path(save.wd, "ecological_processes.pdf"), width = 10, height = 6)

