# 00_shared.R — 여러 스크립트에서 함께 쓰는 설정과 함수
# 사용법: 스크립트 맨 위에서 source("R/00_shared.R")

suppressPackageStartupMessages(library(tidyverse))

plink2 <- "bin/plink2.exe"

# 대륙을 항상 같은 순서로 그린다
super_levels  <- c("AFR", "EUR", "SAS", "EAS", "AMR")
super_colours <- c(AFR = "#2a78d6", EUR = "#eb6834", SAS = "#1baf7a",
                   EAS = "#eda100", AMR = "#e87ba4")

# 최근 혼혈 집단. PEL의 포함 여부는 Day 7에 판단한다
admixed_pops <- c("ASW", "ACB", "MXL", "PUR", "CLM")


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


# 그룹별 평균, 표준편차, 인원을 한 표로 만드는 함수
# 예: summarize_groups(het, Ho, c(SUPER, POP))
# 2단계 단일세포 분석에서도 그룹 이름만 바꿔 쓸 수 있다
summarize_groups <- function(data, value, group) {
  data |>
    group_by(across({{ group }})) |>
    summarise(
      mean = mean({{ value }}),
      sd   = sd({{ value }}),
      n    = n(),
      .groups = "drop"
    )
}
