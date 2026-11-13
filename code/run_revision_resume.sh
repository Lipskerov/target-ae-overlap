#!/bin/sh
# Amendment 2 — full re-run after internal review. Each step logs to logs/rev_*.log; stops at the first failure.
set -e
cd "$(dirname "$0")/.."
PY=python3; L=logs; FQ=data/derived/faers_q; S=code/python
TARGETS=all $PY $S/faers_q_profiles.py b > $L/rev_profiles_alltargets.log 2>&1
echo "[$(date +%H:%M)] H1";         : > $L/rev_H1.log
COHORT_DIR=$FQ/b MEASURES=shared_signals,jaccard $PY $S/primary_inference.py $FQ/b/profile_clean_prr_k20.csv R2_b_prr_clean >> $L/rev_H1.log 2>&1
COHORT_DIR=$FQ/b MEASURES=shared_signals,jaccard $PY $S/primary_inference.py $FQ/b/profile_clean2_ic025_k20.csv R2_b_ic025_clean2 >> $L/rev_H1.log 2>&1
for v in clean_ic025 clean2_prr clean_ror025 clean2_ror025; do COHORT_DIR=$FQ/b $PY $S/primary_inference.py $FQ/b/profile_${v}_k20.csv R2_b_$(echo $v | sed 's/\(.*\)_\(.*\)/\2_\1/') 999 1000 >> $L/rev_H1.log 2>&1; done
for e in c d e; do for v in clean_prr clean2_ic025; do COHORT_DIR=$FQ/$e $PY $S/primary_inference.py $FQ/$e/profile_${v}_k20.csv R2_${e}_$(echo $v | sed 's/\(.*\)_\(.*\)/\2_\1/') 999 1000 >> $L/rev_H1.log 2>&1; done; done
for v in clean_prr clean2_ic025; do COHORT_DIR=$FQ/b_alltargets $PY $S/primary_inference.py $FQ/b_alltargets/profile_${v}_k20.csv R2_balltargets_$(echo $v | sed 's/\(.*\)_\(.*\)/\2_\1/') 999 1000 >> $L/rev_H1.log 2>&1; done
echo "[$(date +%H:%M)] pairs_all";  $PY $S/pairs_all.py $FQ/b/profile_clean_prr_k20.csv R2_prr_clean > $L/rev_pairs_prr.log 2>&1
$PY $S/pairs_all.py $FQ/b/profile_clean2_ic025_k20.csv R2_ic025_clean2 > $L/rev_pairs_ic025.log 2>&1
echo "[$(date +%H:%M)] comention";  MEASURE=prr FLAG=clean $PY $S/comention.py > $L/rev_comention_prr.log 2>&1
MEASURE=ic025 FLAG=clean2 $PY $S/comention.py > $L/rev_comention_ic025.log 2>&1
echo "[$(date +%H:%M)] asof";       MEASURE=prr FLAG=clean $PY $S/asof_q.py > $L/rev_asof_prr.log 2>&1
MEASURE=ic025 FLAG=clean2 $PY $S/asof_q.py > $L/rev_asof_ic025.log 2>&1
echo "[$(date +%H:%M)] trials";     rm -f $FQ/pt_vocabulary.csv; $PY $S/ctgov_profiles.py > $L/rev_ctgov_profiles.log 2>&1
C=data/derived/ctgov/profile; : > $L/rev_H4.log
MEASURES=shared_signals COHORT_DIR=$C $PY $S/primary_inference.py $C/profile_clean_k20.csv R2_ctgov_all >> $L/rev_H4.log 2>&1
MEASURES=shared_signals COHORT_DIR=$C EXCLUDE_PAIRS=$C/shared_condition_pairs.csv $PY $S/primary_inference.py $C/profile_clean_k20.csv R2_ctgov_no_shared_condition >> $L/rev_H4.log 2>&1
MEASURES=shared_signals COHORT_DIR=$C EXCLUDE_PAIRS=$C/shared_trial_pairs.csv $PY $S/primary_inference.py $C/profile_clean_k20.csv R2_ctgov_no_shared_trial >> $L/rev_H4.log 2>&1
$PY $S/ctgov_within_trial.py > $L/rev_ctgov_within.log 2>&1
$PY $S/positive_control.py $FQ/b/profile_clean_prr_k20.csv R2_prr_clean >> $L/rev_H4.log 2>&1
$PY $S/positive_control.py $FQ/b/profile_clean2_ic025_k20.csv R2_ic025_clean2 >> $L/rev_H4.log 2>&1
echo "[$(date +%H:%M)] DONE"
