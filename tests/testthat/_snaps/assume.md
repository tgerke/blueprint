# scenario print methods are stable

    Code
      print(assume_response(0.3))
    Output
      <scenario: response>
        true response rate(s): 0.3 

---

    Code
      print(assume_survival("weibull", shape = 1.5, scale = 4.2))
    Output
      <scenario: survival>
        distribution: weibull (shape = 1.5, scale = 4.2)

