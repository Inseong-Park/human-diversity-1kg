# 04b_excess_het.R — 이형접합이 지나치게 많은 SNP를 찾아 데이터셋 V5 만들기
# 하는 일:
#   1. 기준 B에서 집단별로 SNP마다 관측(ho)과 기대(he) 이형접합도를 계산한다
#   2. 비혼혈 집단의 절반 이상에서 ho - he가 큰 SNP를 찾는다
#   3. 그 SNP만 뺀 데이터셋 V5를 만든다 (친척 제외, MAF, pruning은 B와 같다)
#   4. 빼기 전후의 집단별 He, Ho를 표로 저장한다
# 실행 순서: 02_qc.R → 04_diversity.R → 이 스크립트
# 이유는 notes.md Day 5의 "점검"과 비교 4 참고.

source("R/00_shared.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# 이상 SNP 기준. 관례 값이 아니라 이 분석에서 정한 값이다
excess_cutoff <- 0.1   # ho - he가 이 값보다 크면 "과잉"으로 본다

pops <- read_tsv("work/pops.txt", col_types = cols(.default = "c"))
all_pops <- sort(unique(pops$POP))
non_admixed <- setdiff(all_pops, admixed_pops)   # PEL을 포함한 21개 집단
min_pops <- ceiling(length(non_admixed) / 2)     # 21개 중 11개 이상


# ── 1. 집단별 SNP마다 ho, he 계산 ────────────────────────────
# ID("22:16052080")를 시각으로 잘못 읽지 않도록 개수 열만 숫자로 읽는다
calc_snp_het <- function(pop) {
  out <- paste0("work/gcount_B_", pop)
  run_plink(c("--pfile", "work/B",
              "--keep", paste0("work/keep_", pop, ".txt"),
              "--geno-counts", "--out", out))

  read_tsv(paste0(out, ".gcount"),
           col_types = cols(.default = "c", HOM_REF_CT = "d",
                            HET_REF_ALT_CTS = "d", TWO_ALT_GENO_CTS = "d")) |>
    transmute(
      ID,
      POP = pop,
      n = HOM_REF_CT + HET_REF_ALT_CTS + TWO_ALT_GENO_CTS,
      p = (2 * HOM_REF_CT + HET_REF_ALT_CTS) / (2 * n),
      he = 2 * p * (1 - p),
      ho = HET_REF_ALT_CTS / n
    )
}

snp_het <- map(all_pops, calc_snp_het) |> list_rbind()


# ── 2. 이상 SNP 찾기 ─────────────────────────────────────────
# 혼혈 집단은 빈도 구조가 달라서 기준을 정할 때 쓰지 않는다
excess <- snp_het |>
  filter(POP %in% non_admixed) |>
  group_by(ID) |>
  summarise(n_excess = sum(ho - he > excess_cutoff), .groups = "drop") |>
  filter(n_excess >= min_pops)

write_lines(excess$ID, "work/excess_het_snps.txt")
cat("이상 SNP:", nrow(excess), "개 →  work/excess_het_snps.txt\n")


# ── 3. 데이터셋 V5: B에서 이상 SNP만 뺀다 ─────────────────────
run_plink(c("--pfile", "work/B",
            "--exclude", "work/excess_het_snps.txt",
            "--make-pgen", "--out", "work/V5"))

count_rows <- function(path) sum(!startsWith(readLines(path, warn = FALSE), "#"))
cat("V5:", count_rows("work/V5.psam"), "명,", count_rows("work/V5.pvar"), "SNP\n")


# ── 4. 표 3b: 빼기 전후 비교 ─────────────────────────────────
table3b <- snp_het |>
  group_by(POP) |>
  summarise(
    He_B  = mean(he),
    Ho_B  = mean(ho),
    He_V5 = mean(he[!ID %in% excess$ID]),
    Ho_V5 = mean(ho[!ID %in% excess$ID]),
    .groups = "drop"
  ) |>
  mutate(gap_B = Ho_B - He_B, gap_V5 = Ho_V5 - He_V5) |>
  arrange(desc(He_B))

write_csv(table3b, "results/tables/table3b_excess_het.csv")
print(table3b, n = Inf)
