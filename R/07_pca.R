# 07_pca.R — 기준 데이터 B의 PCA와 로딩
# 하는 일:
#   1. PLINK로 PC 10개와 SNP별 로딩을 계산한다
#   2. 각 PC의 설명 분산 비율을 계산한다
#   3. PC1-PC2, PC1-PC3, 설명 분산, 로딩 그림을 저장한다

source("R/00_shared.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)


# ── 1. PCA ───────────────────────────────────────────────────
# allele-wts: SNP별 로딩도 함께 저장한다
run_plink(c("--pfile", "work/B", "--pca", "allele-wts", "10", "--out", "work/B_pca"))

# 설명 분산 비율의 분모(전체 분산)를 구하기 위해 관계 행렬을 만든다
# PLINK의 --pca 결과에는 상위 10개 고유값만 있어서 전체 합을 따로 계산해야 한다
run_plink(c("--pfile", "work/B", "--make-rel", "square", "--out", "work/B_rel"))


# ── 2. 설명 분산 비율 ────────────────────────────────────────
# 전체 분산 = 관계 행렬의 대각선 합
rel <- as.matrix(data.table::fread("work/B_rel.rel", header = FALSE))
total_var <- sum(diag(rel))
rm(rel)

eigenval <- scan("work/B_pca.eigenval", quiet = TRUE)
pct <- 100 * eigenval / total_var

variance_table <- tibble(PC = seq_along(pct), eigenvalue = eigenval, pct_variance = pct)
write_csv(variance_table, "results/tables/table_pca_variance.csv")
print(variance_table)


# ── 3. 점수에 집단 정보 붙이기 ───────────────────────────────
pops <- read_tsv("work/pops.txt", col_types = cols(.default = "c")) |>
  rename(IID = `#IID`)

scores <- read_tsv("work/B_pca.eigenvec", show_col_types = FALSE) |>
  rename(IID = `#IID`) |>
  mutate(IID = as.character(IID)) |>
  inner_join(pops, by = "IID") |>
  mutate(
    SUPER = factor(SUPER, levels = super_levels),
    group_type = if_else(POP %in% admixed_pops, "최근 혼혈 집단", "그 외 집단")
  )

stopifnot(nrow(scores) == 2491)

# 대륙과 PEL의 평균 위치 (예측 확인용). PC의 부호는 임의이므로 값의 상대 위치만 본다
scores |>
  mutate(group = case_when(POP == "PEL" ~ "PEL", group_type == "최근 혼혈 집단" ~ POP, TRUE ~ as.character(SUPER))) |>
  group_by(group) |>
  summarise(n = n(), PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3), .groups = "drop") |>
  arrange(PC1) |>
  print(n = Inf)


# ── 4. 그림 3: PC 점수 ───────────────────────────────────────
ggsave("results/figures/fig3_pca_pc1_pc2.png",
       plot_pca(scores, "SUPER", pcs = c(1, 2), colours = super_colours,
                shape_by = "group_type", pct = pct, title = "PC1-PC2 (기준 B)"),
       width = 7.5, height = 6, dpi = 200, bg = "white")

ggsave("results/figures/fig3_pca_pc1_pc3.png",
       plot_pca(scores, "SUPER", pcs = c(1, 3), colours = super_colours,
                shape_by = "group_type", pct = pct, title = "PC1-PC3 (기준 B)"),
       width = 7.5, height = 6, dpi = 200, bg = "white")

ggsave("results/figures/fig3_variance.png", plot_variance(pct),
       width = 6, height = 3.5, dpi = 200, bg = "white")


# ── 5. 그림 3-1: 로딩 ────────────────────────────────────────
# 파일에는 SNP마다 대립유전자별로 두 줄이 있다 (부호만 반대). ALT 쪽 줄만 쓴다
allele_wts <- read_tsv("work/B_pca.eigenvec.allele",
                       col_types = cols(.default = "c", PC1 = "d", PC2 = "d", PC3 = "d", PC4 = "d")) |>
  filter(A1 == ALT) |>
  transmute(pos = as.numeric(sub(".*:", "", ID)), PC1, PC2, PC3, PC4)

stopifnot(nrow(allele_wts) == 8551)

loadings_long <- allele_wts |>
  pivot_longer(starts_with("PC"), names_to = "PC", values_to = "weight")

ggsave("results/figures/fig3b_loadings.png",
       plot_loadings(loadings_long, x_name = "22번 염색체 위치 (Mb)",
                     title = "PC1~PC4 로딩 (기준 B)"),
       width = 8, height = 7, dpi = 200, bg = "white")
