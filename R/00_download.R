# 1000 Genomes Phase 3 데이터 다운로드
options(timeout = 3600)  # 큰 파일이라 시간 제한을 1시간으로 늘림

base <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
files <- c(
  # 서버의 현재 버전은 v5b (예전 v5a는 교체되어 없음)
  "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
  "integrated_call_samples_v3.20130502.ALL.panel",
  "20140625_related_individuals.txt"
)

dir.create("data", showWarnings = FALSE)
for (f in files) {
  dest <- file.path("data", f)
  if (!file.exists(dest)) {
    message("다운로드 중: ", f)
    ok <- tryCatch({
      download.file(paste0(base, f), dest, mode = "wb")
      TRUE
    }, error = function(e) { message("실패: ", conditionMessage(e)); FALSE })
    # 실패하면 덜 받아진 파일을 지워서 다음 실행 때 다시 받도록 함
    if (!ok && file.exists(dest)) file.remove(dest)
  }
}

# 받은 파일 크기 확인 (MB)
round(file.size(file.path("data", files)) / 1e6, 1)