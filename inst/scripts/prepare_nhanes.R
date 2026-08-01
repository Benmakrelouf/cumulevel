## ---------------------------------------------------------------------
##  Rebuild the shipped 'nhanes_smoking' data set from source.
##
##  The prepared data are already available as data(nhanes_smoking); this
##  script documents exactly how they were derived, so that the coding of
##  pack-years and of the mortality endpoint can be audited or altered.
##
##  Both source files are public domain. This script downloads the NHANES
##  components; the mortality file must be fetched separately from
##  https://www.cdc.gov/nchs/data-linkage/mortality-public.htm
##
##  Requires the 'haven' package to read the SAS transport files.
## ---------------------------------------------------------------------
library(cumulevel)
stopifnot(requireNamespace("haven", quietly = TRUE))

base <- "https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/"
get  <- function(f) haven::read_xpt(file.path(base, f))

demo <- get("DEMO_D.XPT")   # age, sex, income-to-poverty ratio
bmx  <- get("BMX_D.XPT")    # body mass index
smq  <- get("SMQ_D.XPT")    # smoking module

## Mortality file: fixed-width, see the NCHS layout document.
lmf <- read.fwf("NHANES_2005_2006_MORT_2019_PUBLIC.dat",
                widths = c(14, 1, 1, 3, 1, 1, 1, 4, 8, 8, 3, 3),
                col.names = c("SEQN", "eligstat", "mortstat", "ucod_leading",
                              "diabetes", "hyperten", "dodqtr", "dodyear",
                              "wgt_new", "sa_wgt_new", "permth_int",
                              "permth_exm"))

d <- merge(demo[, c("SEQN", "RIDAGEYR", "RIAGENDR", "INDFMPIR")],
           bmx[, c("SEQN", "BMXBMI")], by = "SEQN")
d <- merge(d, smq, by = "SEQN")
d <- merge(d, lmf[, c("SEQN", "eligstat", "mortstat", "permth_int")],
           by = "SEQN")

## Pack-years = (years smoked x cigarettes per day) / 20.
## Never-smokers (fewer than 100 lifetime cigarettes) are assigned zero.
never  <- d$SMQ020 == 2
start  <- d$SMD030                      # age of regular onset
stop_a <- ifelse(d$SMQ040 %in% c(1, 2), d$RIDAGEYR, d$SMD055)
cpd    <- ifelse(d$SMQ040 %in% c(1, 2), d$SMD070, d$SMD057)
years  <- pmax(0, stop_a - start)
py     <- ifelse(never, 0, years * cpd / 20)

## Five ordered levels: never, <=10, <=20, <=40, >40.
d$packyear_level <- cut(py, breaks = c(-Inf, 0, 10, 20, 40, Inf),
                        labels = FALSE)

## 10-year all-cause mortality. The 2005-2006 cycle is chosen because the
## 2019 linked file gives every participant at least ~156 months of
## follow-up, so no one is administratively censored before the horizon
## and the binary endpoint is legitimate.
d$died10 <- as.integer(d$mortstat == 1 & d$permth_int <= 120)

d <- subset(d, eligstat == 1 & RIDAGEYR >= 20)
keep <- c("died10", "packyear_level", "RIDAGEYR", "RIAGENDR", "BMXBMI",
          "INDFMPIR")
d <- na.omit(d[, keep])
names(d) <- c("died10", "packyear_level", "age", "sex", "bmi", "pir")

fit <- weighted_level(died10 ~ packyear_level | age + sex + bmi + pir,
                      data = d, J = 3, seed = 2005)
print(summary(fit))
plot(fit)

## NOTE: this analysis is unweighted for the NHANES complex survey design.
## A design-based extension is straightforward but is not implemented here.
