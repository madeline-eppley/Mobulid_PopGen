```R
> print(pairwise_results)
  Pop1 Pop2 Phi_ST_AMOVA p_value_AMOVA Fst_STACKS Phi_ST_STACKS
1 EAST  GAL     -0.00047        0.7219     0.0624      -0.00305
2 EAST  CTR     -0.00080        0.6977     0.0709      -0.00577
3  GAL  CTR     -0.00191        0.7961     0.0786      -0.00942
```

```bash
(base) [eppley.m@explorer-01 populations_gal]$ cat populations.phistats_summary.tsv
# Phi_st Means
	EAST	GAL	CTR	WST
EAST		-0.00304561	-0.00576946	-0.0852982
GAL			-0.00941724	-0.0832728
CTR				-0.0721796

# Fst' Means
	EAST	GAL	CTR	WST
EAST		-0.00239075	-0.00387333	-0.0218634
GAL			-0.00708876	-0.0316982
CTR				-0.0500769

# Dxy Means
	EAST	GAL	CTR	WST
EAST		0.00349326	0.00357786	0.00386717
GAL			0.00365231	0.00401476
CTR				0.00435992
(base) [eppley.m@explorer-01 populations_gal]$ cat populations.fst_summary.tsv
	EAST	GAL	CTR	WST
EAST		0.0623994	0.0709004	0.0997668
GAL			0.0785559	0.115853
CTR				0.150145
```
