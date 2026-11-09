# 01_prep.R — 샘플 정리
# 하는 일:
#   1. 샘플별 집단 정보(panel)와 친척 목록을 읽는다
#   2. VCF에 실제로 들어 있는 샘플과 맞춰 본다
#   3. PLINK에 넘길 파일들을 work/에 만든다
#   4. 집단별 인원표(표 0)를 results/tables/에 저장한다

suppressPackageStartupMessages(library(tidyverse))

vcf_file     <- "data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
panel_file   <- "data/integrated_call_samples_v3.20130502.ALL.panel"
related_file <- "data/20140625_related_individuals.txt"

# 최근 혼혈 집단 (회귀분석에서 따로 다룰 집단)
# PEL은 혼혈이지만 원주민 조상 비율이 높아 따로 표시하고, 포함 여부는 Day 7에 판단한다
admixed_pops <- c("ASW", "ACB", "MXL", "PUR", "CLM")

dir.create("work", showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)


# ── 1. panel: 샘플별 집단·대륙·성별 ─────────────────────────────
# 머리글 줄 끝에만 빈 칸 2개가 더 있어서(6칸 vs 4칸) read_tsv는 경고를 낸다
# read.delim은 칸 수가 달라도 조용히 읽으므로 이걸로 읽고 4개 열만 고른다
panel <- read.delim(panel_file, colClasses = "character") |>
  as_tibble() |>
  select(sample, pop, super_pop, gender)


# ── 2. 친척 목록 ─────────────────────────────────────────────
# 확인 결과: 이 목록의 31명은 2,504명을 만들 때 "이미 제외된" 샘플이다.
#   예) HG00124(목록에 있음)는 VCF에 없고, 그 친척 HG00119는 VCF에 남아 있다.
# 2,504명 안에 남아 있는 친척은 03_kinship.R에서 KING 친족계수로 직접 찾아 work/related.txt를 만든다.
# 머리글이 "Sample "처럼 공백을 포함해서 열 이름을 직접 붙인다
related <- read_tsv(related_file, col_types = cols(.default = "c"),
                    col_names = c("sample", "pop", "gender", "reason"), skip = 1) |>
  mutate(across(everything(), str_trim))


# ── 3. VCF에 들어 있는 샘플 이름 ─────────────────────────────────
# VCF 앞부분의 설명 줄(##)을 지나 "#CHROM"으로 시작하는 줄에 샘플 이름이 있다
# 파일 전체(205MB)를 읽지 않고 그 줄까지만 읽는다
con <- gzfile(vcf_file, "r")
repeat {
  line <- readLines(con, n = 1)
  if (startsWith(line, "#CHROM")) break
}
close(con)
vcf_samples <- strsplit(line, "\t")[[1]][-(1:9)]  # 앞 9개 열은 위치·대립유전자 정보


# ── 4. 서로 맞는지 확인 ───────────────────────────────────────
cat("VCF 샘플 수:           ", length(vcf_samples), "\n")
cat("panel 샘플 수:         ", nrow(panel), "\n")
cat("panel에만 있는 샘플:   ", length(setdiff(panel$sample, vcf_samples)), "\n")
cat("VCF에만 있는 샘플:     ", length(setdiff(vcf_samples, panel$sample)), "\n")
cat("친척 목록 샘플 수:     ", nrow(related), "\n")
cat("그중 VCF에 있는 샘플:  ", sum(related$sample %in% vcf_samples), "\n")


# ── 5. 샘플 표 하나로 합치기 ──────────────────────────────────
samples <- panel |>
  filter(sample %in% vcf_samples) |>
  mutate(
    related = sample %in% related$sample,
    admixed = pop %in% admixed_pops
  )

write_tsv(samples, "work/samples.tsv")


# ── 6. PLINK용 파일 ──────────────────────────────────────────
# PLINK 2는 VCF를 읽을 때 샘플 이름을 IID로 쓴다. "#IID" 머리글을 붙이면 그 열로 인식한다

# 6-1. 집단 정보 (F_ST 계산 등에 사용)
samples |>
  transmute(`#IID` = sample, POP = pop, SUPER = super_pop) |>
  write_tsv("work/pops.txt")

# 6-2. 집단별 샘플 목록 (친척 제외) — 집단마다 따로 계산할 때 사용
for (p in unique(samples$pop)) {
  samples |>
    filter(pop == p, !related) |>
    transmute(`#IID` = sample) |>
    write_tsv(file.path("work", paste0("keep_", p, ".txt")))
}


# ── 7. 표 0: 집단별 인원 (친척 제외 전후) ─────────────────────────
table0 <- samples |>
  group_by(super_pop, pop) |>
  summarise(
    n_total    = n(),
    n_related  = sum(related),
    n_analysis = sum(!related),
    admixed    = first(admixed),
    .groups = "drop"
  ) |>
  arrange(factor(super_pop, levels = c("AFR", "EUR", "SAS", "EAS", "AMR")), pop)

write_csv(table0, "results/tables/table0_sample_counts.csv")

print(table0, n = Inf)
cat("\n대륙별 합계 (친척 제외 후):\n")
print(table0 |> group_by(super_pop) |> summarise(n = sum(n_analysis)))
