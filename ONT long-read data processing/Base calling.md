# 1.Base calling

Base_calling.sh
```
#!/bin/bash
ont-guppy/bin/guppy_basecaller --version           >>version-check.log
ont-guppy/bin/guppy_basecaller  --print_workflows  >>config-check.log
ont-guppy/bin/guppy_basecaller	--device cuda:0,1 --num_callers 10   --ipc_threads 20 --gpu_runners_per_device 4  --chunks_per_runner 1664	--flowcell	FLO-PRO002	--kit	SQK-LSK109	--recursive	--qscore_filtering	--verbose_logs	--input_path	/fast5_path	--save_path	/xxx

```









