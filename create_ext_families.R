t1 <- Sys.time() # (For measuring time of script)
setwd("../")
library(data.table)

civil_status <- readRDS("01_data/civil_status_231025.RDS")
setkey(civil_status, pair_id)

rel_pairs <- readRDS("01_data/rel_pairs.RDS")
table(rel_pairs$type, rel_pairs$r)


# Remove UZ twins
rel_pairs <- rel_pairs[r!=.75]


## Coerce Opposite-Sex relatives to always have the female in the leftmost column
os_pairs <- rbind(rel_pairs[sex==0],rel_pairs[sex==0,.(i_id = j_id, j_id = i_id, r, type, sex)])
os_pairs[, i_female := i_id %in% civil_status$female_id]
os_pairs <- os_pairs[i_female == TRUE]
os_pairs <- os_pairs[,i_female := NULL]

rel_pairs <- rbind(rel_pairs[sex!=0],os_pairs)



########################################################
## FIND ALL EXT FAMILY UNITS
########################################################

RELmale <- rel_pairs[sex==1]
RELfemale <- rel_pairs[sex==2]
RELoppos <- rel_pairs[sex==0]

# Attach marriage id for all individuals
# Here, we allow multiple marriages per person, which we will later 
# sort by number of available relatives

RELmale <- merge(
  RELmale, civil_status[, .(male_id, c2id = pair_id)],
  by.x = "i_id", by.y = "male_id", all.x = TRUE, allow.cartesian = TRUE
)

RELmale <- merge(
  RELmale, civil_status[, .(male_id, c3id = pair_id)],
  by.x = "j_id", by.y = "male_id", all.x = TRUE, allow.cartesian = TRUE
)



RELfemale <- merge(
  RELfemale, civil_status[, .(female_id, c2id = pair_id)],
  by.x = "i_id", by.y = "female_id", all.x = TRUE, allow.cartesian = TRUE
)

RELfemale <- merge(
  RELfemale, civil_status[, .(female_id, c3id = pair_id)],
  by.x = "j_id", by.y = "female_id", all.x = TRUE, allow.cartesian = TRUE
)



RELoppos <- merge(
  RELoppos, civil_status[, .(female_id, c2id = pair_id)],
  by.x = "i_id", by.y = "female_id", all.x = TRUE, allow.cartesian = TRUE
)

RELoppos <- merge(
  RELoppos, civil_status[, .(male_id, c3id = pair_id)],
  by.x = "j_id", by.y = "male_id", all.x = TRUE, allow.cartesian = TRUE
)



########################################################
## Attach siblings-in-law 
# Only consider same-sex full siblings

#################################
# Female partners of male relatives
eligble_female_relatives <-  rel_pairs[type == "full sibs" & sex != 0]

# Force female_id to be "i_id"
eligble_female_relatives[,rel_id := paste0(i_id, j_id)]
eligble_female_relatives <- rbind(eligble_female_relatives,eligble_female_relatives[,.(i_id = j_id, j_id = i_id, r, type, sex, rel_id)])
eligble_female_relatives[civil_status, j_pair_id := i.pair_id, on = c("j_id" = "female_id")]
eligble_female_relatives <- eligble_female_relatives[!is.na(j_pair_id)]

########################
# Male partners of female relatives
eligble_male_relatives <-  rel_pairs[type == "full sibs" & sex != 0]

# Force male-id to be "i_id"
eligble_male_relatives[,rel_id := paste0(i_id, j_id)]
eligble_male_relatives <- rbind(eligble_male_relatives,eligble_male_relatives[,.(i_id = j_id, j_id = i_id, r, type, sex, rel_id)])
eligble_male_relatives[civil_status, j_pair_id := i.pair_id, on = c("j_id" = "male_id")]
eligble_male_relatives <- eligble_male_relatives[!is.na(j_pair_id)]


##############################
# Attach the relatives

# First, make temporary id (id of partner of twin)
RELmale[civil_status, temp2_id := i.female_id, on = c("c2id" = "pair_id")]
RELmale[civil_status, temp3_id := i.female_id, on = c("c3id" = "pair_id")]

RELfemale[civil_status, temp2_id := i.male_id, on = c("c2id" = "pair_id")]
RELfemale[civil_status, temp3_id := i.male_id, on = c("c3id" = "pair_id")]

RELoppos[civil_status, temp2_id := i.male_id,   on = c("c2id" = "pair_id")]
RELoppos[civil_status, temp3_id := i.female_id, on = c("c3id" = "pair_id")]


# Then, join final pair_id on the temporary id

RELmale <- merge(
  RELmale,
  eligble_female_relatives[, .(i_id, c1id = j_pair_id)],
  by.x = "temp2_id",
  by.y = "i_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)

RELmale <- merge(
  RELmale,
  eligble_female_relatives[, .(i_id, c4id = j_pair_id)],
  by.x = "temp3_id",
  by.y = "i_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)


RELfemale <- merge(
  RELfemale,
  eligble_male_relatives[, .(i_id, c1id = j_pair_id)],
  by.x = "temp2_id",
  by.y = "i_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)

RELfemale <- merge(
  RELfemale,
  eligble_male_relatives[, .(i_id, c4id = j_pair_id)],
  by.x = "temp3_id",
  by.y = "i_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)



RELoppos <- merge(
  RELoppos,
  eligble_male_relatives[, .(i_id, c1id = j_pair_id)],
  by.x = "temp2_id",
  by.y = "i_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)

RELoppos <- merge(
  RELoppos,
  eligble_female_relatives[, .(i_id, c4id = j_pair_id)],
  by.x = "temp3_id",
  by.y = "i_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)




## ---------------------------------------------------------------------------
## 7.  COMBINE INTO ONE DATA TABLE
## ---------------------------------------------------------------------------



all_units <- rbind(RELmale, RELfemale, RELoppos)

all_units[type == "twins" & r == 1 & sex == 1, type := "MZm"]
all_units[type == "twins" & r == 1 & sex == 2, type := "MZf"]
all_units[type == "twins" & r == .75 & sex == 1, type := "UZm"]
all_units[type == "twins" & r == .75 & sex == 2, type := "UZf"]
all_units[type == "twins" & r == .5 & sex == 1, type := "DZm"]
all_units[type == "twins" & r == .5 & sex == 2, type := "DZf"]
all_units[type == "twins" & r == .5 & sex == 0, type := "DZmf"]

all_units[type == "full sibs" & sex == 1, type := "FSm"]
all_units[type == "full sibs" & sex == 2, type := "FSf"]
all_units[type == "full sibs" & sex == 0, type := "FSmf"]


all_units <- all_units[,.(c1id,c2id,c3id, c4id, type, sex)]

# Clean up temporary objects
rm(RELmale, RELfemale, RELoppos,eligble_female_relatives, eligble_male_relatives)

table(all_units$type)







## ---------------------------------------------------------------------------
## 7.  Priority-aware de-duplication  (pair-id version ??? fast)
## ---------------------------------------------------------------------------
# Count completeness
all_units[, n_complete := rowSums(!is.na(.SD)), .SDcols = c("c1id", "c2id", "c3id", "c4id")]


  ## ---------- 7A.  build the priority score  (unchanged, vectorised) ----
all_units[, priority := (n_complete-2)]
  
all_units[type %in% c("MZm","MZf"), priority := priority+100]
all_units[type %in% c("DZm","DZf","DZmf"), priority := priority+50]
all_units[type %in% c("UZm","UZf"), priority := priority+25]



## ------------------------------------------------------------------
## Order by priority (highest first)
## ------------------------------------------------------------------
all_units[,randomizer := runif(nrow(all_units))] #Randomises priority within priority group
setorder(all_units, -priority, randomizer) # Set row order by priority
all_units[,randomizer:=NULL]
all_units[,n_complete:=NULL]


all_units <- all_units[(is.na(c1id) | c1id != c3id)] #Removes "double cousin creators" (where my wife's sister is my brother's wife)
all_units <- all_units[(is.na(c4id) | c4id != c2id)] #Removes "double cousin creators" (where my wife's sister is my brother's wife)
all_units <- all_units[(is.na(c1id) | is.na(c4id) | c4id != c1id)] #Removes other closed loops

# Converts a pair id to two numeric ids (for fast comparison)
split_half <- function(x) {as.numeric(c(substr(x, 4, 10), substr(x, 14, 20)))}

units <- copy(all_units)

## Pre-compute the three-marriage bundle per row  (no NA???s)
units[, bundle := Map(function(a,b,c,d) {split_half(na.omit(c(a,b,c,d)))},c1id, c2id, c3id, c4id)]

# Remove rows where the same person appears twice
units[, has_dup := vapply(bundle, function(x) anyDuplicated(x) > 0, logical(1))]
units <- units[has_dup == FALSE]
units[, has_dup := NULL]


## ------------------------------------------------------------------
## FIND PRE-USED
## ------------------------------------------------------------------


# units$bundle : list-column. Each row is a vector of IDs.
# Goal: go through rows in order, keep a row only if NONE of its IDs have been used
# in any previously kept row.

## ------------------------------------------------------------
## 1) Preprocess list-column into a dense integer representation
## ------------------------------------------------------------

# Pull all IDs from all rows into one long vector.
# use.names = FALSE keeps it flat and fast.
all_ids <- unlist(units$bundle, use.names = FALSE)

# Convert IDs to a compact integer code 1..K.
# match(x, unique(x)) returns 1 for the first unique ID, 2 for the next, etc.
# This avoids using the original IDs as indices if they are sparse or strings.
key <- match(all_ids, unique(all_ids))  # same length as all_ids

# Get the length (number of IDs) in each row of the original list-column.
# This tells us how many IDs belong to row 1, row 2, ...
lens <- lengths(units$bundle)

# Split the long encoded vector back into a list, one element per row.
# rep.int(seq_len(nrow(units)), lens) repeats row index i exactly lens[i] times,
# so split(...) knows which pieces belong to which row.
idx_list <- split(key, rep.int(seq_len(nrow(units)), lens))

## ------------------------------------------------------------
## 2) Greedy selection with O(1) membership checks
## ------------------------------------------------------------

# Make a logical vector "seen" of length = number of distinct IDs.
# seen[j] will be TRUE once we have used ID j in some accepted row.
seen <- rep(FALSE, length(unique(all_ids)))

# Output flag for each row: TRUE = keep, FALSE = drop.
keep <- logical(nrow(units))

system.time({
# Loop over rows in the given order (assumed already sorted by importance).
for (i in seq_along(idx_list)) {
  # v is the integer-coded IDs present in row i.
  v <- idx_list[[i]]
  
  # If ANY of these IDs have been seen before, we must skip this row.
  # seen[v] gives a logical vector for just the IDs in this row.
  if (any(seen[v])) next
  
  # Otherwise we accept the row.
  keep[i] <- TRUE
  
  # And we mark all IDs in this row as used so later rows cannot reuse them.
  seen[v] <- TRUE
}

})


units <- units[keep]         # Keep only accepted rows
units[, bundle := NULL]      # drop helper column
units[,priority:=NULL]


# Check number of marriages
units[, n_complete := rowSums(!is.na(.SD)), .SDcols = c("c1id", "c2id", "c3id", "c4id")]
units[,sum(n_complete)]



## ---------------------------------------------------------------------------
## 8.  RESULT
## ---------------------------------------------------------------------------
units


saveRDS(units, "01_data/ext_fam_units.rds")

t2 <- Sys.time() # (For measuring time of script)

t2-t1





