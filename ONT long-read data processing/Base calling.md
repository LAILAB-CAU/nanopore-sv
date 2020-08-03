# 1.Base calling

Base_calling.sh
```
#!/bin/bash
ont-guppy/bin/guppy_basecaller	--device	cuda:0	--num_callers	8	--ipc_threads	16	--flowcell	FLO-PRO002	--kit	SQK-LSK109	--recursive	--qscore_filtering	--verbose_logs	--input_path	/fast5_path	--save_path	/xxx

```
