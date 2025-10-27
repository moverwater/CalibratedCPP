birthRates = runif(200,0,10)
deathRates = runif(200,0,10)

birthAndDeathRates = data.frame(birthRates = birthRates,
                                deathRates = deathRates)

write.csv(birthAndDeathRates,"~/code/beast_and_friends/CalibratedCoalescentPointProcess/validation/bd_rates.csv")
