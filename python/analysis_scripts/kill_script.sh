#!/bin/bash
job_num=$(squeue|grep 'mbaxter')
to_kill=${job_num:14:4}
$(scancel "$to_kill")
