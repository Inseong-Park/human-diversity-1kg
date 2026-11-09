# 02_qc.R — PLINK 2 전처리
# 하는 일:
#   Part A (Day 3): VCF를 PLINK 2 형식(work/raw)으로 변환하고, 단계별로 변이가 몇 개 남는지 기록한다
#   Part B (Day 4 이후): 03_kinship.R로 친척 목록을 만든 뒤 분석 조건 B, V1~V4 데이터셋을 만든다

suppressPackageStartupMessages(library(tidyverse))

plink2   <- "bin/plink2.exe"
vcf_file <- "data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

dir.create("work", showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# PLINK를 실행하고, 실패하면 로그 마지막 부분을 보여주며 멈추는 함수
run_plink <- function(args) {
  out <- system2(plink2, args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cat(tail(out, 20), sep = "\n")
    stop("PLINK 실행 실패 (위 메시지 확인)")
  }
  invisible(out)
}


# ── Part A. VCF → PLINK 2 형식 ───────────────────────────────
# --snps-only just-acgt : 삽입·결실(indel)과 구조변이는 빼고, A/C/G/T 한 글자짜리 SNP만 남긴다
# --max-alleles 2       : 한 줄에 대립유전자가 3개 이상인 변이를 뺀다
# --set-all-var-ids @:# : 변이 ID를 "염색체:위치"로 통일한다 (원래 ID에는 rs번호와 '.'이 섞여 있음)
# --rm-dup exclude-all  : ID가 겹치는 변이, 즉 같은 위치에 여러 줄로 적힌 다중대립 변이를 모두 뺀다
# 변환은 몇 분 걸리므로, 이미 결과가 있으면 건너뛴다

if (!file.exists("work/raw.pgen")) {
  message("VCF 변환 중... (몇 분 걸립니다)")
  run_plink(c(
    "--vcf", vcf_file,
    "--snps-only", "just-acgt",
    "--max-alleles", "2",
    "--set-all-var-ids", "@:#",
    "--rm-dup", "exclude-all",
    "--make-pgen",
    "--out", "work/raw"
  ))
} else {
  message("work/raw.pgen이 이미 있어서 변환을 건너뜁니다")
}


# ── 단계별 변이 수 (표 1) ─────────────────────────────────────
# PLINK 로그에서 각 단계 숫자를 읽는다
# 로그 첫 부분에 한글 폴더 경로가 있어서, 글자 변환 경고가 나지 않게 useBytes = TRUE로 찾는다
log_text <- paste(readLines("work/raw.log", warn = FALSE), collapse = " ")
grab <- function(pattern) as.numeric(sub(pattern, "\\1", regmatches(log_text, regexpr(pattern, log_text, useBytes = TRUE)), useBytes = TRUE))

n_scanned  <- grab("--vcf: ([0-9]+) variants scanned")
n_snp_bi   <- grab("--snps-only, ([0-9]+) remaining")
n_final    <- grab("([0-9]+) variants remaining after main filters")

table1 <- tibble(
  step = c("1. VCF 전체 변이",
           "2. 이중대립 SNP만 남김 (indel 41,031 + 다중대립 6,280 + 구조변이 782 제거)",
           "3. 같은 위치 중복 제외 → work/raw"),
  n_variants = c(n_scanned, n_snp_bi, n_final)
) |>
  mutate(removed = lag(n_variants) - n_variants,
         pct_of_total = round(100 * n_variants / n_scanned, 2))

write_csv(table1, "results/tables/table1_variant_counts.csv")
print(table1)

# 샘플 수 확인
psam <- read_tsv("work/raw.psam", col_types = cols(.default = "c"))
cat("\n샘플 수:", nrow(psam), "\n")
# 참고: VCF에는 성별 정보가 없어서 PLINK가 모두 'ambiguous'로 읽는다.
#       22번 염색체(상염색체)만 쓰므로 분석에는 영향이 없다.


# ── Part B. 분석 조건별 데이터셋 ─────────────────────────────
# 순서: 02_qc.R Part A → 03_kinship.R → 02_qc.R Part B
if (!file.exists("work/related.txt")) {
  stop("work/related.txt가 없습니다. 먼저 R/03_kinship.R을 실행하세요")
}

# 조건 하나로 데이터셋 하나를 만드는 함수
# --remove           : 친척 목록의 샘플을 뺀다
# --maf              : 덜 흔한 대립유전자 빈도가 이 값 미만인 SNP를 뺀다
# --indep-pairwise   : 50개 창에서 5개씩 옮기며 r² > 0.2인 SNP 쌍의 한쪽을 뺀다
make_dataset <- function(name, maf, prune, remove_related) {
  base <- c("--pfile", "work/raw", "--maf", maf)
  if (remove_related) base <- c(base, "--remove", "work/related.txt")

  extract <- character(0)
  if (prune) {
    run_plink(c(base, "--indep-pairwise", "50", "5", "0.2",
                "--out", paste0("work/", name, "_prune")))
    extract <- c("--extract", paste0("work/", name, "_prune.prune.in"))
  }

  run_plink(c(base, extract, "--make-pgen", "--out", paste0("work/", name)))
}

conditions <- tribble(
  ~dataset, ~maf,   ~prune, ~remove_related, ~description,
  "B",      "0.05", TRUE,   TRUE,            "기준: 1촌 이상 제외, MAF 0.05, pruning",
  "V1",     "0.05", TRUE,   FALSE,           "친척 포함",
  "V2",     "0.01", TRUE,   TRUE,            "MAF 0.01",
  "V3",     "0.05", FALSE,  TRUE,            "pruning 생략",
  "V4",     "0.01", FALSE,  TRUE,            "MAF 0.01 + pruning 생략"
)

pwalk(select(conditions, name = dataset, maf, prune, remove_related), make_dataset)

# 만들어진 파일에서 샘플 수와 변이 수를 센다 (# 로 시작하는 줄은 머리글)
count_rows <- function(path) sum(!startsWith(readLines(path, warn = FALSE), "#"))

table1b <- conditions |>
  mutate(
    n_samples  = map_int(dataset, \(d) count_rows(paste0("work/", d, ".psam"))),
    n_variants = map_int(dataset, \(d) count_rows(paste0("work/", d, ".pvar")))
  ) |>
  select(dataset, description, n_samples, n_variants)

write_csv(table1b, "results/tables/table1b_dataset_summary.csv")
print(table1b)
