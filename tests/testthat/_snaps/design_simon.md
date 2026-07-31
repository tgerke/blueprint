# print and protocol text output are stable

    Code
      print(design)
    Output
      <Simon two-stage design: optimal> 
        Hypotheses: H0 p = 0.2 vs. H1 p = 0.4 (one-sided)
        Stage 1: enroll 13; stop for futility if <= 3 responses
        Stage 2: enroll 30 more (43 total); reject H0 if >= 13 responses
        Exact operating characteristics:
          type I error 0.0496 | power 0.8002
          PET(H0) 0.747 | E[N | H0] 20.6

---

    Code
      cat(draft_protocol_text(design))
    Output
      A Simon optimal two-stage design (Simon, 1989) will be used. The null hypothesis that the true response rate is 20% will be tested against a one-sided alternative. In the first stage, 13 patients will be accrued. If there are 3 or fewer responses in these 13 patients, the study will be stopped. Otherwise, 30 additional patients will be accrued for a total of 43. The null hypothesis will be rejected if 13 or more responses are observed in 43 patients. This design yields a type I error rate of 0.050 and power of 0.800 when the true response rate is 40%. Under the null hypothesis, the probability of early termination is 0.75 and the expected sample size is 20.6 patients.

