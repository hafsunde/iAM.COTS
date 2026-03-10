# Extended explanation of the iAM-COTS structural equation model

This note describes the model implemented in `R/iam.cots.fun.R` in plain language, with extra focus on **how the parent-offspring correlation is decomposed**.

## 1) What this model is trying to do

At a high level, the iAM-COTS model asks:

> When parents and offspring are similar on a trait, how much of that similarity is due to direct parent effects, shared genes, family environment, and assortative mating processes?

It is an extension of a children-of-twins style design. The model uses families where the parent generation includes related pairs (for example MZ twins, DZ twins, or siblings), their partners, and offspring in both nuclear families. This structure gives enough contrasts to separate several pathways that are hard to distinguish in ordinary parent-child data.

## 2) Who is in the model

The observed variables are ordered as:

1. Parent/sibling in family 1
2. Parent/sibling in family 2
3. Partner of person 1
4. Partner of person 2
5. Child 1 in family 1
6. Child 2 in family 1
7. Child 1 in family 2
8. Child 2 in family 2

So each row of data corresponds to a small kinship network (two related adults, their partners, and up to four children).

## 3) Latent components in the parent generation

For the parental phenotype, the model decomposes variance into:

- **A1**: additive genetic influences (`a1` path)
- **C1**: sibling-shared family environment (`c1` path)
- **T1**: twin-specific shared environment (`t1` path; omitted when not relevant)
- **E1**: non-shared environment (`e1` path)

A gene-environment correlation parameter **`w`** links A and C components (rather than forcing them to be independent). This allows genetic liability and family environmental liability to co-occur.

## 4) Intergenerational pathways (parent -> offspring)

Three main transmission/confounding paths are estimated:

- **`p`**: phenotypic transmission (the parental phenotype influencing offspring outcome)
- **`a1p`**: transmission/confounding through parental genetic factors
- **`c1p`**: transmission/confounding through parental shared-environment liability

In the code, the **raw parent-offspring covariance** is built from these pieces:

- phenotypic pathway: `var1 * p`
- genetic pathway: `a1 * 0.5 * a1p`
- environmental pathway: `c1 * c1p`
- two rGE cross-terms: `a1*w*c1p` and `c1*w*0.5*a1p`

These are combined in `covPO`.

## 5) Offspring-specific variance

The offspring phenotype also has variance not attributable to parent-generation liabilities:

- **A2** (`a2`)
- **C2** (`c2`)
- **E2** (`e2`)

So the offspring outcome is not treated as only a projection of parental factors; it has its own generation-specific residual architecture.

## 6) How assortative mating is represented

The model includes assortment using:

- a strength parameter **`mu`** (`mu_param` in parameter output), and
- a **sorting factor** `S` that partners assort on.

If `indirect.assortment = TRUE`, the sorting factor can have its own A/C/T/E loadings (`a1s`, `c1s`, `t1s`, `e1s`) instead of being identical to the focal parental phenotype.

Key intuition:

- partners become similar because they sort on `S`
- that similarity feeds back into intergenerational resemblance through extra covariance terms
- these extra terms are tracked separately from the "raw" (non-assortment) paths

## 7) The parent-offspring correlation decomposition (main focus)

The model computes a 3 x 4 decomposition object called `est_PO`, with rows:

- `raw` (without assortative-mating-induced inflation)
- `AMinf` (extra covariance induced by assortment)
- `total` (`raw + AMinf`)

and columns:

1. Parent Offspring Correlation (overall)
2. Phenotypic
3. Passive Genetic
4. Passive Environmental

All entries are standardized by dividing by `sqrt(var1)*sqrt(var2)`, so they are on the correlation scale.

### 7.1 Raw part (no AM inflation)

The raw decomposition corresponds to:

- **Overall raw correlation**: `covPO / sqrt(var1*var2)`
- **Phenotypic raw**: `(var1*p) / sqrt(var1*var2)`
- **Passive genetic raw**: `((a1 + c1*w)*0.5*a1p) / sqrt(var1*var2)`
- **Passive environmental raw**: `((c1 + a1*w)*c1p) / sqrt(var1*var2)`

The two "passive" components each contain a direct term and an rGE-weighted cross-term via `w`.

### 7.2 AM-induced part

The AM-induced increment is built from `covPS*mu*SortPO`, where `SortPO` is the covariance between the **parental sorting factor** and the offspring phenotype.

This increment is then split analogously into:

- **AM-induced phenotypic**: `covPS*mu*covPS*p`
- **AM-induced passive genetic**: `covPS*mu*(a1s + c1s*w)*0.5*a1p`
- **AM-induced passive environmental**: `covPS*mu*(c1s + a1s*w)*c1p`

again standardized to the correlation scale.

### 7.3 Total decomposition

The final reported decomposition is simply:

- `total = raw + AMinf`

component-wise and for the overall parent-offspring correlation.

This is useful because it distinguishes:

1. resemblance that would exist even without mate assortment in this parameterization, and
2. additional resemblance amplified by assortative mating.

## 8) Why twin/sibling structure is crucial

The model has multiple relatedness groups (for example MZ, DZ, FS, optionally HS/UZ). Their expected covariances are built using relatedness-specific genetic correlations (`1`, `FS_rG`, `UZ_rG`, `HS_rG`) and shared-environment assumptions. These cross-family contrasts are what help separate:

- parent-to-child phenotypic transmission (`p`),
- inherited/passive pathways (`a1p`, `c1p`), and
- AM-induced covariance inflation.

Without those contrasts, many of these terms would be statistically confounded.

## 9) Practical interpretation guidance

When reading estimates from `est_PO`:

- Start with **`total` / Parent Offspring Correlation** as the model-implied r(parent, child).
- Compare **`raw` vs `AMinf`** to gauge how much assortment contributes.
- Within each row, compare **Phenotypic**, **Passive Genetic**, and **Passive Environmental** columns.
- Remember that `w` means genetic and environmental liabilities are partly entangled; passive components include rGE cross-terms by design.

## 10) In one sentence

The iAM-COTS implementation decomposes parent-offspring similarity into direct phenotypic, passive genetic, and passive environmental pathways, and then explicitly splits each of those into a non-assortment part and an assortative-mating-induced increment.
