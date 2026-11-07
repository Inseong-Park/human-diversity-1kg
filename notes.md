# 작업 노트

## 메모
- 함수 이름 충돌: 나중에 불러온 패키지가 우선. 헷갈리면 dplyr::filter()처럼 패키지 이름을 붙인다.
- .RData는 저장하지 않는다. 결과는 항상 스크립트로 다시 만든다.

## Day 1 — 환경 준비

### 설치한 것
- R 4.4.3 (원래 있음), RStudio 2026.10.1
- R 패키지: tidyverse(표 정리·그림), data.table(큰 표 읽기), geosphere(지구 위 거리), ape(계통수)
- PLINK 2 v2.0.0-b.1.3 → bin/plink2.exe
  - 1.9에는 없는 KING 친족계수, 모든 집단 쌍 F_ST 계산이 있고 큰 VCF를 빠르게 처리함
- PLINK 1.9 (수업 때 쓰던 것) → LD decay(--r2)에만 사용
- R에서 PLINK 부르기: system2("bin/plink2.exe", "--version")
  - 모든 PLINK 명령을 R 스크립트 안에서 실행해 과정 전체를 기록으로 남긴다

### 정한 규칙
- .RData 저장 끔 (Tools → Global Options → Never)
  - 결과는 항상 스크립트를 처음부터 다시 돌려서 만든다 → 재현 가능성
- 함수 이름 충돌: 나중에 불러온 패키지가 우선. 헷갈리면 dplyr::filter()처럼 패키지 이름을 붙인다
- 폴더: data(원본) / work(중간 결과) / R(스크립트) / results(그림·표) / bin(실행 파일)

### git / GitHub
- git = 내 컴퓨터의 변경 기록, GitHub = 그 기록을 올려두는 곳
- user.name은 커밋에 찍히는 작성자 이름(실명), user.email이 GitHub 계정과 연결하는 기준
- .gitignore로 data/, work/, bin/은 올리지 않음 (크고, 스크립트로 다시 만들 수 있음)
- 작업 단위마다: git add → git commit -m "..." → git push
  - add = 상자에 담기, commit = 봉하고 송장 붙이기(내 컴퓨터에만), push = 발송
- 저장소는 Private으로 시작. Day 14에 GitHub 프로필 정리 후 Public으로 전환
- 첫 커밋 메시지 오타(Initail)를 git commit --amend + git push --force-with-lease로 수정
  - 이미 push한 커밋은 원칙적으로 고치지 않는다. 혼자 쓰는 저장소의 첫 커밋이라 예외
- RStudio Terminal은 Ctrl+V가 잘 안 됨 → Shift+Insert 또는 오른쪽 클릭 Paste

## Day 2 — 데이터 받기와 샘플 정리

### 데이터 (1000 Genomes Phase 3, release/20130502)
- ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz (205.6MB)
- integrated_call_samples_v3.20130502.ALL.panel
- 20140625_related_individuals.txt

### 배운 것 1: 공개 데이터도 버전이 바뀐다
- 처음에 v5a로 받으려다 404 오류. 서버 목록을 확인하니 v5b로 교체되어 있었다
- README에 사용한 파일 버전을 정확히 적을 것
- 다운로드 실패 시 덜 받아진 파일을 지우도록 스크립트를 고침 (안 그러면 다음 실행 때 건너뜀)

### 배운 것 2: 파일 이름만 보고 용도를 짐작하면 안 된다
- 친척 목록 31명 중 VCF에 있는 샘플이 0명이었다
- HG00124(목록에 있음)는 VCF에 없고, 그 친척 HG00119만 남아 있음
  → 이 목록은 "앞으로 뺄 샘플"이 아니라 "2,504명을 만들 때 이미 뺀 샘플"의 기록
- 처음 계획한 비교 1(제공된 목록으로 친척 제외)은 아무도 빠지지 않아 의미가 없음
  → Day 4에 KING 친족계수로 2,504명 안에 남은 친척을 직접 찾는 방식으로 변경

### 결과
- VCF 2,504명 = panel 2,504명 (한쪽에만 있는 샘플 0)
- 26개 집단, 집단당 61명(ASW) ~ 113명(GWD)
- 대륙별: AFR 661, EUR 503, EAS 504, SAS 489, AMR 347
- 혼혈 집단 표시: ASW, ACB, MXL, PUR, CLM (PEL은 Day 7에 포함 여부 판단)
- 표 0: results/tables/table0_sample_counts.csv

### 확인할 질문
- 대륙별 인원이 비슷한데, 아프리카만 661명으로 많은 이유는? (혼혈 집단 ASW, ACB가 AFR로 분류됨)

## [비교 이름]
- 바꾼 조건:
- 예측 (실행 전에 작성):
- 실제 결과:
- 예측이 맞았나? 틀렸다면 왜?
- 수업과 연결되는 점: