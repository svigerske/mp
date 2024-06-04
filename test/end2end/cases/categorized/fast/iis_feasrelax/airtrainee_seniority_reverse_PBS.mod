
# Reverse seniority constraints,
# sublimated as post-processing objectives.
# Optimize for all trainees with E[i]==1, ranked by reverse seniority,
# then for all with E[i]==2.
# In this special case the objectives can be aggregated for each value of e
# (see Solve with aggregated preferences.)
maximize ReverseSeniority {e in 1..2, i in I: E[i]==e}:
  sum {t in V[i]: Pr[i, t]==0}
    S[i] * x[i, t]
    suffix objpriority (2-e)*S_range + 1 + S[i] - min {j in I} S[j];
