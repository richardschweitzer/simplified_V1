estimate_TRFs_for_SFs
================
Richard Schweitzer
2023-05-26

Libraries…

``` r
library(data.table)
library(ggplot2)
library(scales)
library(viridis)
```

    ## Loading required package: viridisLite

    ## 
    ## Attaching package: 'viridis'

    ## The following object is masked from 'package:scales':
    ## 
    ##     viridis_pal

``` r
library(cowplot)
library(minpack.lm)
library(pracma)
library(mgcv)
```

    ## Loading required package: nlme

    ## This is mgcv 1.9-1. For overview type 'help("mgcv-package")'.

    ## 
    ## Attaching package: 'mgcv'

    ## The following object is masked from 'package:pracma':
    ## 
    ##     magic

Shall we save the data during this run?

``` r
do_save <- FALSE
```

# 1. Replicate the findings by Burr & Morrone (1993)

The model and convolution operation in Fig. 2.

``` r
source("exp_damp_sin.R")

# create temporal response function
framedur <- 1000/120
t_res <- framedur/12
trf <- exp_damp_sin(x_step = t_res)
trf_time <- trf[[2]]
trf <- trf[[1]]
plot(trf_time, trf)
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# create a sequence of SOA=60 ms and two positive pulses, and pos/neg pulses
t_stim_on <- 10
soa <- 60
stim_time <- seq(1, 500, by = t_res)
# ... two positive pulses
stim_1 <- rep(0, times = length(stim_time))
stim_1[stim_time >= t_stim_on & stim_time <= t_stim_on+framedur] <- 1
stim_1[stim_time >= t_stim_on+soa & stim_time <= t_stim_on+soa+framedur] <- 1
# ... one positive, one negative
stim_2 <- rep(0, times = length(stim_time))
stim_2[stim_time >= t_stim_on & stim_time <= t_stim_on+framedur] <- 1
stim_2[stim_time >= t_stim_on+soa & stim_time <= t_stim_on+soa+framedur] <- (-1)
plot(stim_time, stim_2, col = "red", type = "l")
lines(stim_time, stim_1, col = "black")
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
# now convolute -- we should get Fig. 2 of the paper
resp_1 <- zapsmall(convolve(stim_1, rev(trf), type = "open"))[1:length(stim_1)]
resp_2 <- zapsmall(convolve(stim_2, rev(trf), type = "open"))[1:length(stim_2)]
plot(stim_time, resp_2, type = "l", col = "red")
lines(stim_time, resp_1, col = "black")
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

This above solution uses convolution, but let’s try the solution
proposed in the paper… Turns out that to reproduce the absolute values
in Fig. 1, we *a0* parameter must be *1.8e3*, not *1.8e5* as stated in
the paper.

``` r
two_pulse_response <- function(t_ms, soa_1 = 0, soa_2, polarity, K = 1,
                               a0 = 1.8e3, a1 = 13, a2 = 5, a3 = 27.3 ) {
  res <- K * ( exp_damp_sin_fun(t_ms = t_ms - soa_1, a0 = a0, a1 = a1, a2 = a2, a3 = a3) + 
                polarity * exp_damp_sin_fun(t_ms = t_ms - soa_2, a0 = a0, a1 = a1, a2 = a2, a3 = a3) )
  return(res)
}

# again, replication of Fig. 2
soa_60ms_pospos <- two_pulse_response(t_ms = stim_time, soa_1 = 10, soa_2 = 70, polarity = 1)
soa_60ms_posneg <- two_pulse_response(t_ms = stim_time, soa_1 = 10, soa_2 = 70, polarity = -1)
plot(stim_time, soa_60ms_posneg, type = "l", col = "red")
lines(stim_time, soa_60ms_pospos, col = "black")
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# probability summation with integrate
beta_now <- 4
integrate(f = function(x) { 
  abs(two_pulse_response(t_ms = x, soa_1 = 0, soa_2 = 60, polarity = 1))^beta_now }, lower = 0, upper = Inf)$value^(1/beta_now)
```

    ## [1] 52.15366

And now the model results from Figure 1B - should be the result of the
probability summation. We can use the code above and iterate over SOAs

``` r
# beta parameter for probability summation:
beta <- 4
# not preallocate
prob_sum_results <- NULL
# loop over SOAs
for (soa in seq(8.3, 300, by = 8.3)) {
  ## the convolution solution
  trf_now <- exp_damp_sin(x_step = t_res, a0 = 1.8e3)
  trf_now <- trf_now[[1]]
  # the stimuli (see above)
  stim_time <- seq(1, 1000, by = t_res)
  stim_1 <- rep(0, times = length(stim_time))
  stim_1[stim_time >= t_stim_on & stim_time <= t_stim_on+framedur] <- 1 / 12
  stim_1[stim_time >= t_stim_on+soa & stim_time <= t_stim_on+soa+framedur] <- 1 / 12
  stim_2 <- rep(0, times = length(stim_time))
  stim_2[stim_time >= t_stim_on & stim_time <= t_stim_on+framedur] <- 1 / 12
  stim_2[stim_time >= t_stim_on+soa & stim_time <= t_stim_on+soa+framedur] <- (-1) / 12
  # response
  resp_1 <- zapsmall(convolve(stim_1, rev(trf_now), type = "open"))[1:length(stim_1)]
  resp_2 <- zapsmall(convolve(stim_2, rev(trf_now), type = "open"))[1:length(stim_2)]
  # plot(stim_time, resp_2, type = "l", col = "red")
  # lines(stim_time, resp_1, col = "black")
  # probability summation
  sum_1 <- trapz(x = stim_time, y = abs(resp_1)^beta)^(1/beta)
  sum_2 <- trapz(x = stim_time, y = abs(resp_2)^beta)^(1/beta)
  ## the integration-based solution:
  inta_1 <- integrate(f = function(x) { 
    abs(two_pulse_response(t_ms = x, soa_1 = 0, soa_2 = soa, polarity = 1))^beta }, lower = 0, upper = Inf)$value^(1/beta)
  inta_2 <- integrate(f = function(x) { 
    abs(two_pulse_response(t_ms = x, soa_1 = 0, soa_2 = soa, polarity = -1))^beta }, lower = 0, upper = Inf)$value^(1/beta)
  # save results
  prob_sum_results <- rbind(prob_sum_results, 
                            data.table(soa = soa, sum_pospos = sum_1, sum_posneg = sum_2, 
                                       inta_pospos = inta_1, inta_posneg = inta_2))
}
# now this should match Fig. 1B - DB's luminance data
# ... sum
plot(prob_sum_results$soa, prob_sum_results$sum_posneg, col = "red", type = "l", 
     ylim = c(min(prob_sum_results[ ,c(2,3)]), max(prob_sum_results[ ,c(2,3)]) ))
lines(prob_sum_results$soa, prob_sum_results$sum_pospos, col = "black", type = "l")
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# ... integrate
plot(prob_sum_results$soa, prob_sum_results$inta_posneg, col = "red", type = "l", 
     ylim = c(min(prob_sum_results[ ,c(4,5)]), max(prob_sum_results[ ,c(4,5)]) ))
lines(prob_sum_results$soa, prob_sum_results$inta_pospos, col = "black", type = "l")
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
# ... integrate in loglog coordinates
ggplot(data = prob_sum_results, aes(x = soa, y = inta_pospos, color = "pos-pos")) + 
  geom_line() + 
  geom_line(data = prob_sum_results, aes(x = soa, y = inta_posneg, color = "pos-neg")) + 
  scale_y_log10() + annotation_logticks(sides = "l") + 
  theme_classic()
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

Now, finally, let’s try to predict the flicker sensitivity data. Burr &
Morrone used flicker from 0.5 to 50 Hz with a Gaussian temporal envelope
of SD=230 ms. We’ll create that data, too.

``` r
source('create_flicker_stim.R')

# parameters
(temporal_gaussian_sd <- 28 * framedur)
```

    ## [1] 233.3333

``` r
pres_dur <- 10*temporal_gaussian_sd
# not preallocate
flicker_sens_results <- NULL
# loop over SOAs
for (tfreq in exp(seq(log(1), log(40), length.out = 40)) ) {
  ## solution 1
  trf_now <- exp_damp_sin(x_step = t_res, a0 = 1.8e3)
  trf_now <- trf_now[[1]] #/ sum(trf_now[[1]])
  # gaussian envelope
  stim_time <- seq(0, pres_dur, by = t_res)
  gaussian <- dnorm(x = stim_time, mean = pres_dur/2, sd = temporal_gaussian_sd)
  gaussian <- gaussian / max(gaussian)
  x <- create_flicker_stim(pres_dur, t_res, tfreq) # flicker stimulus
  stim_tfreq <- x #* gaussian
  #plot(stim_time, stim_tfreq, type = "l")
  # response
  resp_tfreq <- zapsmall(convolve(stim_tfreq, rev(trf_now), type = "open"))[1:length(stim_tfreq)]
  #plot(stim_time, resp_tfreq, type = "l", col = "red")
  # probability summation
  sum_tfreq <- trapz(x = stim_time, y = abs(resp_tfreq)^beta)^(1/beta)
  # save results
  flicker_sens_results <- rbind(flicker_sens_results, 
                                data.table(tfreq = tfreq, sum_tfreq = sum_tfreq, inta_tfreq = NaN))
}
# sum approach
ggplot(data = flicker_sens_results, aes(x = tfreq, y = sum_tfreq)) + 
  geom_line() + 
  annotation_logticks() + 
  scale_x_log10() + 
  scale_y_log10(#breaks = trans_breaks("log10", function(x) 10^x), 
                #labels = trans_format("log10", math_format(10^.x))
    ) + 
  theme_classic()
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

That worked so so, but sufficiently well.

# 2. The IRF used by Bergen & Wilson (1984)

This function here may be better suited, actually:

``` r
source('bergen_wilson.R')

# looks as predicted with standard parameters:
plot(bergen_wilson(x_step = t_res)[[2]], bergen_wilson(x_step = t_res)[[1]])
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# 3. Fit flicker sensitivity data to estimate TRFs

Let’s combine the the two functions for the IRF and the stimulus
generation.

``` r
source("flicker_to_sens.R")

# test this (Burr & Morrone):
flicker_to_sens(tfreqs = c(1, 40), 
                a0 = 1.8e3, a1 = 13, a2 = 5, a3 = 27.3, 
                t_res = 1000/1440, pres_dur = 5000, ramp_sd = 28 * (1000/120))
```

    ## [1] 1597.5961  101.2486

``` r
# ... and Bergen & Wilson:
flicker_to_sens(tfreqs = c(1, 40), 
                a0 = 1.8e3, a1 = 0.92, a2 = 9.5, a3 = 4.0, a4 = 1.5, 
                t_res = 1000/1440, pres_dur = 5000, ramp_sd = 28 * (1000/120))
```

    ## [1] 9076.0688  731.7458

Load digitized flicker sensitivity data from Kelly (1969) and try to fit
the data with the above function

``` r
Kelly69 <- read.csv(file = "Kelly_1969_data_Fig8.csv")
str(Kelly69)
```

    ## 'data.frame':    36 obs. of  3 variables:
    ##  $ tfreq     : num  1.05 1.54 3.1 4.03 6.18 ...
    ##  $ modulation: num  0.00528 0.00455 0.00351 0.0031 0.00271 ...
    ##  $ SF        : num  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...

``` r
Kelly69$sens <- 1 / Kelly69$modulation
setDT(Kelly69)

# have a look:
ggplot(data = Kelly69, aes(x = tfreq, y = sens, color = as.factor(SF) )) + 
  geom_point() + geom_line() + 
  scale_x_log10() + scale_y_log10() + theme_classic()
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
# let's try to fit
Kelly69[ , fit_nls_burr := NaN]
Kelly69[ , fit_nls_bergen := NaN]
for (sf_now in sort(unique(Kelly69$SF))) {
  # Burr & Morrone
  print(paste(sf_now, "Burr & Morrone"))
  fit1 <- nlsLM(data = Kelly69, subset = SF==sf_now, 
                formula = log(sens) ~ log(flicker_to_sens(tfreq, 
                                                          a0 = a0, a1 = a1, a2 = a2, a3 = a3, 
                                                          t_res = 1000/1440, 
                                                          pres_dur = 8 * 28 * (1000/120), 
                                                          ramp_sd = 28 * (1000/120)
                                                          ) ), 
                start = c(a0 = 1.8e3, a1 = 13, a2 = 5, a3 = 27.3), trace = TRUE, 
                lower = c(a0 = 0, a1 = 0, a2 = 0, a3 = 0) )
  coef_fit1 <- coefficients(fit1)
  # Bergen & Wilson
  print(paste(sf_now, "Bergen & Wilson"))
  fit2 <- nlsLM(data = Kelly69, subset = SF==sf_now, 
                formula = log(sens) ~ log(flicker_to_sens(tfreq, 
                                                          a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4,
                                                          t_res = 1000/1440, 
                                                          pres_dur = 8 * 28 * (1000/120), 
                                                          ramp_sd = 28 * (1000/120)
                                                          ) ), 
                start = c(a0 = 400, a1 = 1, a2 = 10, a3 = 4.0, a4 = 1.5), trace = TRUE, 
                lower = c(a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0) )
  coef_fit2 <- coefficients(fit2)
  # plot both functions
  plot(exp_damp_sin(a0 = coef_fit1[1], a1 = coef_fit1[2], a2 = coef_fit1[3], a3 = coef_fit1[4])[[1]], main = sf_now)
  points(bergen_wilson(A = coef_fit2[1], B = coef_fit2[2], tau = coef_fit2[3], n = coef_fit2[4], k = coef_fit2[5])[[1]], col = "red")
  # model comparison
  print(anova(fit1, fit2))
  # predict based on fits fits
  Kelly69[SF==sf_now, fit_nls_burr := exp(predict(fit1))]
  Kelly69[SF==sf_now, fit_nls_bergen := exp(predict(fit2))]
}
```

    ## [1] "0.25 Burr & Morrone"
    ## It.    0, RSS =    76.5608, Par. =       1800         13          5       27.3
    ## It.    1, RSS =    14.8863, Par. =    674.219    13.6347    4.33258    34.4454
    ## It.    2, RSS =    10.1939, Par. =    122.137    15.9179    6.08229    48.1541
    ## It.    3, RSS =    1.60497, Par. =    242.939    16.6184    6.72034    56.8018
    ## It.    4, RSS =  0.0817568, Par. =    315.017    16.0202    6.50137    52.0973
    ## It.    5, RSS =   0.068165, Par. =     320.69    16.9866    7.98198    53.5619
    ## It.    6, RSS =  0.0655953, Par. =    322.135    16.9368    7.78534    53.1485
    ## It.    7, RSS =  0.0655953, Par. =    322.135    16.9368    7.78534    53.1485
    ## [1] "0.25 Bergen & Wilson"
    ## It.    0, RSS =    97.1779, Par. =        400          1         10          4        1.5
    ## It.    1, RSS =    96.7042, Par. =    393.146   0.923224    9.98563     4.0171    1.46931
    ## It.    2, RSS =    95.3624, Par. =    389.836   0.925764     9.9808    4.02544      1.454
    ## It.    3, RSS =    92.7226, Par. =    383.319   0.930753    9.96747    4.04163    1.42435
    ## It.    4, RSS =    87.5671, Par. =    370.909    0.94039    9.94116    4.07528    1.36671
    ## It.    5, RSS =    77.3916, Par. =    347.522   0.959124     9.8978    4.13437    1.24742
    ## It.    6, RSS =    59.5347, Par. =    305.015   0.987279    9.78198    4.24845    1.04924
    ## It.    7, RSS =    29.4076, Par. =      225.1    1.01376    9.52135    4.47034   0.728234
    ## It.    8, RSS =    2.13898, Par. =     110.99   0.958685    8.62811    4.72783   0.302935
    ## It.    9, RSS =    1.83475, Par. =    129.388   0.992629    9.76979    3.14171   0.114149
    ## It.   10, RSS =  0.0834931, Par. =     142.18   0.984816    9.92705    3.09483   0.127079
    ## It.   11, RSS =  0.0679615, Par. =     134.38   0.983042    10.0295    3.03879   0.134638
    ## It.   12, RSS =  0.0679396, Par. =    134.275   0.982913    10.0303    3.04098   0.134431
    ## It.   13, RSS =  0.0675854, Par. =     134.39   0.982902    10.0319    3.04272   0.134609
    ## It.   14, RSS =    0.06735, Par. =    134.523    0.98289    10.0328     3.0456   0.134656
    ## It.   15, RSS =  0.0669433, Par. =    134.784   0.982919     10.032    3.05051   0.134873
    ## It.   16, RSS =  0.0663935, Par. =    135.189   0.983025    10.0456    3.05855   0.135394
    ## It.   17, RSS =  0.0663914, Par. =    135.039   0.983044    10.0413    3.05899   0.135385
    ## It.   18, RSS =   0.066372, Par. =    135.025    0.98306    10.0424     3.0578   0.135431
    ## It.   19, RSS =   0.066372, Par. =    135.025    0.98306    10.0424     3.0578   0.135431

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreq, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 8 * 28 * (1000/120), ramp_sd = 28 * (1000/120)))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreq, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 8 * 28 * (1000/120), ramp_sd = 28 * (1000/120)))
    ##   Res.Df Res.Sum Sq Df      Sum Sq F value Pr(>F)
    ## 1      8   0.065595                              
    ## 2      7   0.066372  1 -0.00077667 -0.0819      1
    ## [1] "0.5 Burr & Morrone"
    ## It.    0, RSS =    57.4267, Par. =       1800         13          5       27.3
    ## It.    1, RSS =    12.4882, Par. =    809.477    13.4102    4.51317    34.0531
    ## It.    2, RSS =   0.560017, Par. =    347.034    13.5129    4.41685    39.9052
    ## It.    3, RSS =   0.190335, Par. =     399.81    13.7175    5.01311    48.1894
    ## It.    4, RSS =   0.154081, Par. =    410.167    13.8118    4.72873    50.5703
    ## It.    5, RSS =   0.148728, Par. =    421.888    14.0054    4.83655    51.7637
    ## It.    6, RSS =   0.148669, Par. =    423.373    14.0285    4.84537    51.8239
    ## It.    7, RSS =   0.148608, Par. =     422.97    14.0301    4.85233      51.78
    ## It.    8, RSS =   0.148606, Par. =    423.378    14.0352    4.85663    51.7818
    ## It.    9, RSS =   0.148596, Par. =    423.067    14.0328    4.85641      51.78
    ## It.   10, RSS =   0.148596, Par. =    422.907    14.0318    4.85637    51.7798
    ## It.   11, RSS =   0.148595, Par. =    422.989    14.0323    4.85638    51.7793
    ## It.   12, RSS =   0.148595, Par. =    422.992    14.0328    4.85628    51.7786
    ## It.   13, RSS =   0.148595, Par. =     422.97    14.0328    4.85631    51.7784
    ## It.   14, RSS =   0.148594, Par. =    422.959    14.0327    4.85631    51.7783
    ## It.   15, RSS =   0.148594, Par. =    422.957    14.0327    4.85633    51.7783
    ## It.   16, RSS =   0.148594, Par. =    422.959    14.0327    4.85634    51.7784
    ## It.   17, RSS =   0.148594, Par. =     422.96    14.0327    4.85634    51.7783
    ## It.   18, RSS =   0.148594, Par. =    422.959    14.0327    4.85634    51.7783
    ## It.   19, RSS =   0.148594, Par. =    422.959    14.0327    4.85634    51.7783
    ## It.   20, RSS =   0.148594, Par. =    422.959    14.0327    4.85634    51.7783
    ## It.   21, RSS =   0.148594, Par. =    422.959    14.0327    4.85634    51.7783
    ## [1] "0.5 Bergen & Wilson"
    ## It.    0, RSS =    75.7318, Par. =        400          1         10          4        1.5
    ## It.    1, RSS =    74.9839, Par. =    387.902   0.896849    9.97976    4.03186    1.44649
    ## It.    2, RSS =    72.8759, Par. =    381.784   0.901982    9.97132    4.04747    1.42273
    ## It.    3, RSS =    68.6879, Par. =    369.884   0.911585    9.95334     4.0792    1.37387
    ## It.    4, RSS =    60.6605, Par. =    347.732   0.928157    9.91788    4.13747    1.27427
    ## It.    5, RSS =    45.9112, Par. =    310.673   0.960036    9.83238    4.24557    1.08188
    ## It.    6, RSS =    21.2996, Par. =    239.719    1.01637    9.61948    4.43645   0.754885
    ## It.    7, RSS =    1.73247, Par. =    121.798    1.02902      8.305    4.74851   0.268708
    ## It.    8, RSS =   0.579641, Par. =     77.437    1.05075    8.27906    4.11537    0.40709
    ## It.    9, RSS =  0.0799342, Par. =    73.2634    1.06545    8.37464    3.97881   0.490901
    ## It.   10, RSS =  0.0706245, Par. =     68.429    1.07367    8.15797    4.07811   0.558166
    ## It.   11, RSS =  0.0700086, Par. =    60.7023    1.08777    8.08077    4.07123   0.619567
    ## It.   12, RSS =  0.0697997, Par. =    59.5129    1.08899    8.04402    4.10449    0.65107
    ## It.   13, RSS =   0.068404, Par. =    56.7023    1.08922    8.34046    3.90723   0.643743
    ## It.   14, RSS =   0.068404, Par. =    56.7023    1.08922    8.34046    3.90723   0.643743

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->

    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreq, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 8 * 28 * (1000/120), ramp_sd = 28 * (1000/120)))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreq, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 8 * 28 * (1000/120), ramp_sd = 28 * (1000/120)))
    ##   Res.Df Res.Sum Sq Df  Sum Sq F value  Pr(>F)  
    ## 1      8   0.148594                             
    ## 2      7   0.068404  1 0.08019  8.2061 0.02418 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "2 Burr & Morrone"
    ## It.    0, RSS =    50.9869, Par. =       1800         13          5       27.3
    ## It.    1, RSS =    10.5201, Par. =    759.527     12.801    4.63794    33.4657
    ## It.    2, RSS =    5.58118, Par. =    190.906    14.4051    7.53384    42.5507
    ## It.    3, RSS =    0.25443, Par. =    295.811    16.4975    11.8386    40.0278
    ## It.    4, RSS =   0.236568, Par. =    438.096    12.1349    6.99021    57.9091
    ## It.    5, RSS =   0.135555, Par. =    367.786    14.3893     8.6202    46.6606
    ## It.    6, RSS =   0.103897, Par. =    360.042    14.7142    9.31493     43.883
    ## It.    7, RSS =   0.103844, Par. =    358.425    14.6769    9.30244    43.8845
    ## It.    8, RSS =   0.103712, Par. =     358.12    14.7079    9.30888    43.8759
    ## It.    9, RSS =   0.103704, Par. =    358.544    14.7169    9.31121     43.879
    ## It.   10, RSS =   0.103688, Par. =     358.66    14.7194    9.31439    43.8725
    ## It.   11, RSS =   0.103688, Par. =    358.709    14.7184    9.31448    43.8694
    ## It.   12, RSS =   0.103685, Par. =    358.692    14.7178    9.31462    43.8681
    ## It.   13, RSS =   0.103681, Par. =    358.681     14.717    9.31486    43.8681
    ## It.   14, RSS =   0.103681, Par. =    358.681     14.717    9.31486    43.8681
    ## [1] "2 Bergen & Wilson"
    ## It.    0, RSS =    69.6271, Par. =        400          1         10          4        1.5
    ## It.    1, RSS =    68.1507, Par. =    373.704    0.83294    9.99231    4.07597    1.38217
    ## It.    2, RSS =    63.8932, Par. =    361.104   0.844256    9.99193    4.11349    1.33344
    ## It.    3, RSS =     56.003, Par. =    338.336   0.863682    9.97792    4.17813    1.23521
    ## It.    4, RSS =    41.9843, Par. =    301.566   0.896111    9.96561    4.29078    1.04032
    ## It.    5, RSS =    17.7413, Par. =    238.319   0.945058     9.8843    4.49792   0.679563
    ## It.    6, RSS =    3.54977, Par. =    139.071   0.953075    9.27517     4.5545   0.146994
    ## It.    7, RSS =   0.294307, Par. =    143.544   0.940734    8.96005    3.63797    0.14232
    ## It.    8, RSS =   0.235669, Par. =    141.151   0.942472    9.62012    3.32595   0.137184
    ## It.    9, RSS =   0.163047, Par. =    128.717   0.947315    11.2049     2.8114   0.129399
    ## It.   10, RSS =  0.0962123, Par. =    126.323   0.955712    12.9507    2.53576   0.126074
    ## It.   11, RSS =  0.0541821, Par. =     108.55   0.957826    14.6827    2.29758   0.139082
    ## It.   12, RSS =  0.0462121, Par. =    109.675   0.959924    15.0843    2.29511     0.1414
    ## It.   13, RSS =  0.0460585, Par. =    110.342   0.962662    15.4498    2.28467   0.144906
    ## It.   14, RSS =  0.0403744, Par. =    109.193   0.961861    15.5267    2.24321   0.140421
    ## It.   15, RSS =  0.0353763, Par. =    110.669   0.963412    15.8912    2.20058   0.134343
    ## It.   16, RSS =  0.0342937, Par. =     116.72   0.966231    16.0843     2.1923   0.128244
    ## It.   17, RSS =  0.0312748, Par. =    114.205   0.965239     16.346    2.12238   0.123056
    ## It.   18, RSS =  0.0286178, Par. =    118.098    0.96781     16.675    2.11046   0.120256
    ## It.   19, RSS =   0.028395, Par. =    129.198   0.972726    17.1119    2.08405   0.108578
    ## It.   20, RSS =  0.0265103, Par. =    135.693   0.973242    17.2218    2.06361   0.101162
    ## It.   21, RSS =   0.026413, Par. =    135.429   0.973173    17.2144    2.05454   0.100723
    ## It.   22, RSS =   0.026312, Par. =    135.459   0.973237    17.2125    2.05759    0.10117
    ## It.   23, RSS =  0.0262765, Par. =    135.676   0.973303     17.215    2.05982   0.101576
    ## It.   24, RSS =  0.0262438, Par. =    135.614   0.973357    17.2194    2.06151   0.101684
    ## It.   25, RSS =   0.026243, Par. =    135.577   0.973359    17.2195    2.06124   0.101638
    ## It.   26, RSS =   0.026243, Par. =    135.577   0.973359    17.2195    2.06124   0.101638

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->

    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreq, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 8 * 28 * (1000/120), ramp_sd = 28 * (1000/120)))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreq, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 8 * 28 * (1000/120), ramp_sd = 28 * (1000/120)))
    ##   Res.Df Res.Sum Sq Df   Sum Sq F value   Pr(>F)   
    ## 1      8   0.103681                                
    ## 2      7   0.026243  1 0.077438  20.656 0.002652 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# look at fits
ggplot(data = Kelly69, aes(x = tfreq, y = sens, color = as.factor(SF) )) + 
  geom_point() + 
  geom_line(data = Kelly69, aes(x = tfreq, y = fit_nls_burr, color = as.factor(SF), linetype = "Burr & Morrone")) + 
  geom_line(data = Kelly69, aes(x = tfreq, y = fit_nls_bergen, color = as.factor(SF), linetype = "Bergen & Wilson")) + 
  scale_x_log10() + scale_y_log10() + theme_classic() + 
  labs(x = "Temporal frequency [Hz]", y = "Contrast sensitivity", color = "SF", linetype = "Fit")
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-9-5.png)<!-- -->

Both of the two formulations of the TRF seem to fit the data well.

Now we can try to fit Kelly’s (1979) stabilized data that can be
approximated by a function.

``` r
source('kelly_vel.R') # this is Kelly's function
source('get_best_NLS.R') # for rigorous fitting

Kelly_stabilized <- expand.grid(list(TF = 10^(linspace(log10(0.5), log10(50), 33)), 
                                     SF = 10^(linspace(log10(0.05), log10(10), 31)) 
                                     ))
setDT(Kelly_stabilized)
# apply the factor of 2 for flicker thresholds (see Kelly, 1979)
Kelly_stabilized[ , sens := kelly_vel(sf = SF, v = TF / SF) / 2] 

# this is Kelly's surface:
p_kelly_surface <- ggplot(data = Kelly_stabilized, aes(x = SF, y = TF, z = log10(sens) )) + 
  #geom_raster(interpolate = TRUE) + 
  geom_contour() + geom_contour_filled(bins = 30) + 
  coord_cartesian(expand = FALSE) + 
  scale_x_log10() + scale_y_log10() + 
  annotation_logticks() + 
  #scale_fill_viridis_c(trans = "log10") +
  theme_classic(base_size = 12.5) + theme() + 
  labs(x = "Spatial frequency [cpd]", y = "Temporal frequency [Hz]", fill = "Contrast\nsensitivity")
p_kelly_surface
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# parameters for fitting
do_plot <- FALSE
pres_dur_stabilized <- 4 * 8 * 28 * (1000/120) # 4 * 8 * 28 * (1000/120)
ramp_sd_stabilized <- pres_dur_stabilized / 8
beta_stabilized <- 1.5
update_start <- FALSE # use starting parameters of previous fit? 
start_burr <- c(a0 = 1/4, a1 = 10, a2 = 3, a3 = 5)  # , G = 1.4
start_bergen <- c(a0 = 1/4, a1 = 0.92, a2 = 6, a3 = 10, a4 = 1.5)
start_params_n <- 3 # for rigorous fitting set >1, how many steps between lowest and highest starting value?
# preallocate
Kelly_stabilized[ , fit_nls_burr := NaN]
Kelly_stabilized[ , fit_nls_bergen := NaN]
Kelly_stabilized[ , fit_nls_burr_rigorous := NaN]
Kelly_stabilized[ , fit_nls_bergen_rigorous := NaN]
Kelly_IRFs <- NULL
Kelly_IRF_coefs <- NULL
# start cluster if necessary
if (start_params_n > 1) {
  library(parallel)
  (n_cores <- detectCores()-2) # how many cores?
  library(doParallel)
  par_cluster <- makeCluster(n_cores) # , outfile="" # for verbose parallel output
  registerDoParallel(par_cluster)
  # for some reason, we must manually export our functions...
  clusterExport(cl = par_cluster, 
                varlist = c("exp_damp_sin", "exp_damp_sin_fun", "bergen_wilson", "bergen_wilson_fun", 
                            "create_flicker_stim", "flicker_to_sens") 
  )
} else {
  par_cluster <- NULL
}
```

    ## Loading required package: foreach

    ## Loading required package: iterators

``` r
# fit a IRF to each SF
for (sf_now in sort(unique(Kelly_stabilized$SF))) {
  # Burr & Morrone
  print(paste(sf_now, "Burr & Morrone"))
  if (start_params_n > 1) {
    fit1_stab <- get_best_NLS_foreach(cl = par_cluster, 
                                      df = Kelly_stabilized[SF==sf_now], 
                                      model_form = paste0("log(sens) ~ log(flicker_to_sens(tfreqs = TF, 
                                                           a0 = a0, a1 = a1, a2 = a2, a3 = a3, 
                                                           t_res = 1000/1440, 
                                                           pres_dur = ", pres_dur_stabilized, ", ",
                                                           "ramp_sd = ", ramp_sd_stabilized, ", ",
                                                           "beta = ", beta_stabilized, ", ",
                                                           "G = 1))"), 
                                      use_robust = FALSE, 
                                      model_control = nls.lm.control(maxiter = 200), 
                                      absolute_lowest = c(a0 = 0, a1 = 0, a2 = 0, a3 = 0),
                                      start_params_low = start_burr - start_burr*0.5,
                                      start_params_high = start_burr + start_burr*0.5,
                                      start_params_n = start_params_n, 
                                      debug_mode = FALSE
                                      )
  } 
  fit1_stab_s <- nlsLM(data = Kelly_stabilized, subset = SF==sf_now, 
                       formula = log(sens) ~ log(flicker_to_sens(tfreqs = TF, 
                                                                 a0 = a0, a1 = a1, a2 = a2, a3 = a3, 
                                                                 t_res = 1000/1440, 
                                                                 pres_dur = pres_dur_stabilized, 
                                                                 ramp_sd = ramp_sd_stabilized, 
                                                                 beta = beta_stabilized, 
                                                                 G = 1) ), 
                       #weights = 1/sens, # use weights to account for lowest sensitivities?
                       start = start_burr, 
                       lower = c(a0 = 0, a1 = 0, a2 = 0, a3 = 0),  # , G = 1
                       trace = FALSE, control = nls.lm.control(maxiter = 200)
  )
  coef_fit1_stab <- coefficients(fit1_stab)
  coef_fit1_stab_s <- coefficients(fit1_stab_s)
  if (update_start) {
    start_burr <- coef_fit1_stab_s
  }
  # Bergen & Wilson
  print(paste(sf_now, "Bergen & Wilson"))
  if (start_params_n > 1) {
    fit2_stab <- get_best_NLS_foreach(cl = par_cluster, 
                                      df = Kelly_stabilized[SF==sf_now], 
                                      model_form = paste0("log(sens) ~ log(flicker_to_sens(tfreqs = TF, 
                                                           a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4,
                                                           t_res = 1000/1440, 
                                                           pres_dur = ", pres_dur_stabilized, ", ",
                                                           "ramp_sd = ", ramp_sd_stabilized, ", ",
                                                           "beta = ", beta_stabilized, ", ",
                                                           "G = 1))"), 
                                      use_robust = FALSE, 
                                      model_control = nls.lm.control(maxiter = 200), 
                                      absolute_lowest = c(a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0),
                                      start_params_low = start_bergen - start_bergen*0.5,
                                      start_params_high = start_bergen + start_bergen*0.5,
                                      start_params_n = start_params_n, 
                                      debug_mode = FALSE
                                      )
  } 
  fit2_stab_s <- nlsLM(data = Kelly_stabilized, subset = SF==sf_now, 
                       formula = log(sens) ~ log(flicker_to_sens(tfreqs = TF, 
                                                                 a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4,
                                                                 t_res = 1000/1440, 
                                                                 pres_dur = pres_dur_stabilized, 
                                                                 ramp_sd = ramp_sd_stabilized,
                                                                 beta = beta_stabilized, 
                                                                 G = 1
                       ) ), 
                       #weights = 1/sens, # use weights to account for lowest sensitivities?
                       start = start_bergen, 
                       lower = c(a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0), 
                       trace = FALSE, control = nls.lm.control(maxiter = 200) )
#  plot(exp(predict(fit2_stab)), type = "l")
#  points(Kelly_stabilized[SF==sf_now, sens], col = "red")
  coef_fit2_stab <- coefficients(fit2_stab)
  coef_fit2_stab_s <- coefficients(fit2_stab_s)
  if (update_start) {
    start_bergen <- coef_fit2_stab_s
  }
  # model comparison
  if (start_params_n > 1) {
    print(anova(fit1_stab, fit2_stab))
  }
  print(anova(fit1_stab_s, fit2_stab_s))
  # save parameters
  if (start_params_n > 1) {
    Kelly_IRF_coefs <- rbind(Kelly_IRF_coefs, 
                             data.table(SF = sf_now, model = "Burr & Morrone", fit_type = "rigorous",
                                        a0 = coef_fit1_stab[1], a1 = coef_fit1_stab[2], a2 = coef_fit1_stab[3], a3 = coef_fit1_stab[4], 
                                        a4 = NaN, G = coef_fit1_stab[5], rss = summary(fit1_stab)$sigma))
    Kelly_IRF_coefs <- rbind(Kelly_IRF_coefs, 
                             data.table(SF = sf_now, model = "Bergen & Wilson", fit_type = "rigorous",
                                        a0 = coef_fit2_stab[1], a1 = coef_fit2_stab[2], a2 = coef_fit2_stab[3], a3 = coef_fit2_stab[4], 
                                        a4 = coef_fit2_stab[5], G = coef_fit2_stab[6], rss = summary(fit2_stab)$sigma))
  }
  Kelly_IRF_coefs <- rbind(Kelly_IRF_coefs, 
                           data.table(SF = sf_now, model = "Burr & Morrone", fit_type = "simple",
                                      a0 = coef_fit1_stab_s[1], a1 = coef_fit1_stab_s[2], a2 = coef_fit1_stab_s[3], a3 = coef_fit1_stab_s[4], 
                                      a4 = NaN, G = coef_fit1_stab_s[5], rss = summary(fit1_stab_s)$sigma))
  Kelly_IRF_coefs <- rbind(Kelly_IRF_coefs, 
                           data.table(SF = sf_now, model = "Bergen & Wilson", fit_type = "simple",
                                      a0 = coef_fit2_stab_s[1], a1 = coef_fit2_stab_s[2], a2 = coef_fit2_stab_s[3], a3 = coef_fit2_stab_s[4], 
                                      a4 = coef_fit2_stab_s[5], G = coef_fit2_stab_s[6], rss = summary(fit2_stab_s)$sigma))
  # extract IRF waveforms
  if (start_params_n > 1) {
    irf_stab_fit1 <- exp_damp_sin(x_dur = 300, x_step = t_res, 
                                  a0 = coef_fit1_stab[1], a1 = coef_fit1_stab[2], a2 = coef_fit1_stab[3], a3 = coef_fit1_stab[4])
    irf_stab_fit2 <- bergen_wilson(x_dur = 300, x_step = t_res, 
                                   A = coef_fit2_stab[1], B = coef_fit2_stab[2], tau = coef_fit2_stab[3], 
                                   n = coef_fit2_stab[4], k = coef_fit2_stab[5])
  }
  irf_stab_fit1_s <- exp_damp_sin(x_dur = 300, x_step = t_res, 
                                  a0 = coef_fit1_stab_s[1], a1 = coef_fit1_stab_s[2], a2 = coef_fit1_stab_s[3], a3 = coef_fit1_stab_s[4])
  irf_stab_fit2_s <- bergen_wilson(x_dur = 300, x_step = t_res, 
                                   A = coef_fit2_stab_s[1], B = coef_fit2_stab_s[2], tau = coef_fit2_stab_s[3], 
                                   n = coef_fit2_stab_s[4], k = coef_fit2_stab_s[5])
  if (do_plot) {
    plot(irf_stab_fit1[[2]], irf_stab_fit1[[1]], main = sf_now)
    points(irf_stab_fit2[[2]], irf_stab_fit2[[1]], col = "red")
  }
  # add
  Kelly_IRFs <- rbind(Kelly_IRFs, data.table(irf_time = irf_stab_fit1_s[[2]], 
                                             irf_burr = irf_stab_fit1_s[[1]], 
                                             irf_bergen = irf_stab_fit2_s[[1]],
                                             SF = sf_now, 
                                             fit_type = "simple"))
  if (start_params_n > 1) {
    Kelly_IRFs <- rbind(Kelly_IRFs, data.table(irf_time = irf_stab_fit1[[2]], 
                                               irf_burr = irf_stab_fit1[[1]], 
                                               irf_bergen = irf_stab_fit2[[1]],
                                               SF = sf_now, 
                                               fit_type = "rigorous"))
  }
  # compare fit
  Kelly_stabilized[SF==sf_now, fit_nls_burr := exp(predict(fit1_stab_s))]
  Kelly_stabilized[SF==sf_now, fit_nls_bergen := exp(predict(fit2_stab_s))]
  if (start_params_n > 1) {
    Kelly_stabilized[SF==sf_now, fit_nls_burr_rigorous := exp(predict(fit1_stab))]
    Kelly_stabilized[SF==sf_now, fit_nls_bergen_rigorous := exp(predict(fit2_stab))]
  } else {
    Kelly_stabilized[SF==sf_now, fit_nls_burr_rigorous := NaN]
    Kelly_stabilized[SF==sf_now, fit_nls_bergen_rigorous := NaN]
  }
}
```

    ## [1] "0.05 Burr & Morrone"

    ## Loading required package: assertthat

    ## [1] "0.05 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    10.1598                                
    ## 2     28     0.8084  1 9.3514   323.9 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     17.783                                
    ## 2     28      0.817  1 16.966  581.42 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.0596583179235247 Burr & Morrone"
    ## [1] "0.0596583179235247 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     9.6387                                
    ## 2     28     0.6677  1  8.971  376.19 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    17.6214                                
    ## 2     28     0.6833  1 16.938  694.13 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.0711822979492869 Burr & Morrone"
    ## [1] "0.0711822979492869 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     9.0252                                
    ## 2     28     0.5147  1 8.5105  462.98 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     17.463                                
    ## 2     28      0.522  1 16.941  908.66 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.0849323232317123 Burr & Morrone"
    ## [1] "0.0849323232317123 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     8.4771                                
    ## 2     28     0.3684  1 8.1087  616.25 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    17.3022                                
    ## 2     28     0.3737  1 16.929  1268.5 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.101338390826821 Burr & Morrone"
    ## [1] "0.101338390826821 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     7.8852                                
    ## 2     28     0.2545  1 7.6307  839.61 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    17.1454                                
    ## 2     28     0.2626  1 16.883  1799.8 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.120913558756098 Burr & Morrone"
    ## [1] "0.120913558756098 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     7.4021                                
    ## 2     28     0.1964  1 7.2058  1027.5 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    16.9915                                
    ## 2     28     0.2075  1 16.784  2264.3 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.144269990590721 Burr & Morrone"
    ## [1] "0.144269990590721 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     6.6873                                
    ## 2     28     0.1998  1 6.4875  908.97 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     16.840                                
    ## 2     28      0.206  1 16.634  2261.4 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.172138099309703 Burr & Morrone"
    ## [1] "0.172138099309703 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     6.2287                                
    ## 2     28     0.2299  1 5.9988  730.55 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    16.7075                                
    ## 2     28     0.2349  1 16.473  1963.4 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.205389389107391 Burr & Morrone"
    ## [1] "0.205389389107391 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     5.8739                                
    ## 2     28     0.2861  1 5.5878  546.77 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    16.5996                                
    ## 2     28     0.2924  1 16.307  1561.5 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.245063709469745 Burr & Morrone"
    ## [1] "0.245063709469745 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     5.6178                                
    ## 2     28     0.3694  1 5.2484  397.85 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    10.3656                                
    ## 2     28     0.3897  1  9.976  716.84 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.292401773821287 Burr & Morrone"
    ## [1] "0.292401773821287 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     5.2773                                
    ## 2     28     0.4385  1 4.8388  308.99 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     9.9769                                
    ## 2     28     0.4941  1 9.4827  537.35 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.348883959680657 Burr & Morrone"
    ## [1] "0.348883959680657 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     5.3232                                
    ## 2     28     0.4922  1 4.8311  274.84 5.216e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     9.6437                                
    ## 2     28     0.5980  1 9.0457  423.54 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.416276603700937 Burr & Morrone"
    ## [1] "0.416276603700937 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     5.3956                                
    ## 2     28     0.5244  1 4.8712  260.11 1.051e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     10.850                                
    ## 2     28      0.711  1 10.139  399.28 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.496687239354311 Burr & Morrone"
    ## [1] "0.496687239354311 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     5.9051                                
    ## 2     28     0.5203  1 5.3848  289.81 2.649e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     10.796                                
    ## 2     28      0.813  1  9.983  343.82 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.592630504679146 Burr & Morrone"
    ## [1] "0.592630504679146 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     5.8146                                
    ## 2     28     0.4729  1 5.3417  316.31 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    11.4898                                
    ## 2     28     0.9122  1 10.578  324.69 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.707106781186547 Burr & Morrone"
    ## [1] "0.707106781186547 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29      5.097                                
    ## 2     28      0.389  1  4.708  338.89 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     9.4157                                
    ## 2     28     0.9133  1 8.5024  260.66 1.023e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "0.843696023158145 Burr & Morrone"
    ## [1] "0.843696023158145 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     4.8285                                
    ## 2     28     0.2908  1 4.5377  436.93 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    10.6101                                
    ## 2     28     0.9005  1 9.7096   301.9 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "1.00666971160764 Burr & Morrone"
    ## [1] "1.00666971160764 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     4.9530                                
    ## 2     28     0.2022  1 4.7508  657.78 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    10.6007                                
    ## 2     28     0.8346  1 9.7661  327.64 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "1.20112443398143 Burr & Morrone"
    ## [1] "1.20112443398143 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     6.6725                                
    ## 2     28     0.1375  1 6.5351  1330.9 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    10.7226                                
    ## 2     28     0.7374  1 9.9852  379.14 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "1.43314126696356 Burr & Morrone"
    ## [1] "1.43314126696356 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     6.1097                                
    ## 2     28     0.0987  1 6.0111  1705.8 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    10.9172                                
    ## 2     28     0.6059  1 10.311  476.49 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "1.7099759466767 Burr & Morrone"
    ## [1] "1.7099759466767 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     6.0876                                
    ## 2     28     0.0812  1 6.0065  2072.4 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    11.2170                                
    ## 2     28     0.4633  1 10.754  649.94 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "2.04028577336837 Burr & Morrone"
    ## [1] "2.04028577336837 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     6.6687                                
    ## 2     28     0.0789  1 6.5899  2339.8 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    11.6365                                
    ## 2     28     0.3017  1 11.335  1051.8 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "2.43440034644909 Burr & Morrone"
    ## [1] "2.43440034644909 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     6.6651                                
    ## 2     28     0.1648  1 6.5003  1104.6 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     12.812                                
    ## 2     28      0.165  1 12.647  2145.7 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "2.90464459643197 Burr & Morrone"
    ## [1] "2.90464459643197 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     7.6667                                
    ## 2     28     0.1159  1 7.5508  1824.6 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     13.341                                
    ## 2     28      0.117  1 13.224  3165.6 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "3.46572421577573 Burr & Morrone"
    ## [1] "3.46572421577573 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     8.7465                                
    ## 2     28     0.1498  1 8.5968  1607.2 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    14.5021                                
    ## 2     28     0.1513  1 14.351  2655.4 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "4.13518554200014 Burr & Morrone"
    ## [1] "4.13518554200014 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     9.4411                                
    ## 2     28     0.2571  1 9.1839  1000.1 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    14.9723                                
    ## 2     28     0.2581  1 14.714  1596.1 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "4.93396427474814 Burr & Morrone"
    ## [1] "4.93396427474814 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29     10.053                                
    ## 2     28      0.422  1 9.6305  638.94 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    16.0019                                
    ## 2     28     0.4315  1  15.57  1010.3 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "5.88704018652475 Burr & Morrone"
    ## [1] "5.88704018652475 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    10.4434                                
    ## 2     28     0.6196  1 9.8238  443.94 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    16.1064                                
    ## 2     28     0.6298  1 15.477  688.04 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "7.02421830152518 Burr & Morrone"
    ## [1] "7.02421830152518 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    12.2724                                
    ## 2     28     0.8203  1 11.452  390.93 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    17.1629                                
    ## 2     28     0.8243  1 16.339  554.98 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "8.3810609719326 Burr & Morrone"
    ## [1] "8.3810609719326 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    12.5924                                
    ## 2     28     0.9921  1   11.6  327.39 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    17.4371                                
    ## 2     28     1.0041  1 16.433  458.23 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "10 Burr & Morrone"
    ## [1] "10 Bergen & Wilson"
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = 7466.66666666667, ramp_sd = 933.333333333333, beta = 1.5, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    14.9115                                
    ## 2     28     1.1136  1 13.798  346.91 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Analysis of Variance Table
    ## 
    ## Model 1: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ## Model 2: log(sens) ~ log(flicker_to_sens(tfreqs = TF, a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, t_res = 1000/1440, pres_dur = pres_dur_stabilized, ramp_sd = ramp_sd_stabilized, beta = beta_stabilized, G = 1))
    ##   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
    ## 1     29    17.9775                                
    ## 2     28     1.1227  1 16.855  420.37 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# shutdown cluster
if (!is.null(par_cluster)) {
  stopCluster(par_cluster)
}
# goodness of fit: Burr & Morrone
p_burrmorrone <- ggplot(data = Kelly_stabilized, aes(x = TF, y = sens, color = SF, group = SF )) + 
  #geom_point(size = 2.0) + 
  geom_line(data = Kelly_stabilized, aes(x = TF, y = fit_nls_burr_rigorous, color = SF, group = SF ),
            size = 1.5, alpha = 1) +
  # geom_line(data = Kelly_stabilized, aes(x = TF, y = fit_nls_burr, color = SF, group = SF ), 
  #           size = 1.5, alpha = 1) + 
  scale_x_log10() + 
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) + 
  annotation_logticks() + 
  theme_classic(base_size = 12.5) + theme(legend.position = "bottom") + 
  coord_cartesian(ylim = c(0.01, 200), expand = FALSE) + 
  scale_color_viridis_c(trans = "log10") + 
  labs(x = "Temporal frequency [Hz]", y = "Contrast sensitivity", color = "Spatial frequency [cpd]", 
       title = "IRF model (Burr & Morrone)")
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
p_burrmorrone
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
# goodness of fit: Bergen & Wilson
p_bergenwilson <- ggplot(data = Kelly_stabilized, aes(x = TF, y = sens, color = SF, group = SF )) + 
  #geom_point(size = 2.0) + 
  geom_line(data = Kelly_stabilized, aes(x = TF, y = fit_nls_bergen_rigorous, color = SF, group = SF ),
            size = 1.5, alpha = 1) +
  # geom_line(data = Kelly_stabilized, aes(x = TF, y = fit_nls_bergen, color = SF, group = SF ), 
  #           size = 1.5, alpha = 1) + 
  scale_x_log10() + 
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) + 
  annotation_logticks() + 
  theme_classic(base_size = 12.5) + theme(legend.position = "bottom") + 
  coord_cartesian(ylim = c(0.01, 200), expand = FALSE) + 
  scale_color_viridis_c(trans = "log10") + 
  labs(x = "Temporal frequency [Hz]", y = "Contrast sensitivity", color = "Spatial frequency [cpd]", 
       title = "IRF model (Bergen & Wilson)")
p_bergenwilson
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

``` r
# ground truth
p_groundtruth <- ggplot(data = Kelly_stabilized, aes(x = TF, y = sens, color = SF, group = SF )) + 
  #geom_point(size = 2.0) + 
  geom_line(size = 1.5, alpha = 1) + 
  scale_x_log10() + 
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) + 
  annotation_logticks() + 
  theme_classic(base_size = 12.5) + theme(legend.position = "bottom") + 
  coord_cartesian(ylim = c(0.01, 200), expand = FALSE) + 
  scale_color_viridis_c(trans = "log10") + 
  labs(x = "Temporal frequency [Hz]", y = "Contrast sensitivity", color = "Spatial frequency [cpd]", 
       title = "Kelly's (1979) function")
p_groundtruth
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->

``` r
# combine
p_groundtruth_bergenwilson <- plot_grid(p_groundtruth, 
                                        p_bergenwilson, 
                                        p_burrmorrone, 
                                        nrow = 1)
p_groundtruth_bergenwilson
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-10-5.png)<!-- -->

``` r
## look at parameters
Kelly_IRF_coefs[model=="Bergen & Wilson"]
```

    ##              SF           model fit_type         a0         a1       a2
    ##           <num>          <char>   <char>      <num>      <num>    <num>
    ##  1:  0.05000000 Bergen & Wilson rigorous 0.32658727 0.99999336 6.004004
    ##  2:  0.05000000 Bergen & Wilson   simple 0.13537809 0.99999590 6.011701
    ##  3:  0.05965832 Bergen & Wilson rigorous 0.34503753 1.00000541 6.000671
    ##  4:  0.05965832 Bergen & Wilson   simple 0.14036052 0.99978475 6.064025
    ##  5:  0.07118230 Bergen & Wilson rigorous 0.24226582 0.99999711 5.965182
    ##  6:  0.07118230 Bergen & Wilson   simple 0.14844887 1.00008772 6.000636
    ##  7:  0.08493232 Bergen & Wilson rigorous 0.53964622 0.99999061 6.059688
    ##  8:  0.08493232 Bergen & Wilson   simple 0.15914469 1.00003592 5.996157
    ##  9:  0.10133839 Bergen & Wilson rigorous 0.49370529 0.99998107 6.141780
    ## 10:  0.10133839 Bergen & Wilson   simple 0.14866566 0.99990396 6.044725
    ## 11:  0.12091356 Bergen & Wilson rigorous 1.32068015 0.99998601 6.229332
    ## 12:  0.12091356 Bergen & Wilson   simple 0.14820896 1.00054971 6.090169
    ## 13:  0.14426999 Bergen & Wilson rigorous 1.09187665 0.99691648 6.324813
    ## 14:  0.14426999 Bergen & Wilson   simple 0.19930807 0.98274614 6.235294
    ## 15:  0.17213810 Bergen & Wilson rigorous 0.64203931 1.00920824 6.379588
    ## 16:  0.17213810 Bergen & Wilson   simple 0.22427842 0.97321817 6.310482
    ## 17:  0.20538939 Bergen & Wilson rigorous 1.03057760 0.99185137 6.481696
    ## 18:  0.20538939 Bergen & Wilson   simple 0.22912881 0.96314012 6.427339
    ## 19:  0.24506371 Bergen & Wilson rigorous 1.22241612 0.99088381 6.578476
    ## 20:  0.24506371 Bergen & Wilson   simple 0.16237847 0.92815599 6.349180
    ## 21:  0.29240177 Bergen & Wilson rigorous 0.11584877 0.75744530 3.422642
    ## 22:  0.29240177 Bergen & Wilson   simple 0.17546801 0.91504442 6.435725
    ## 23:  0.34888396 Bergen & Wilson rigorous 0.11989675 0.71095376 3.462611
    ## 24:  0.34888396 Bergen & Wilson   simple 0.26539262 0.93204524 6.718235
    ## 25:  0.41627660 Bergen & Wilson rigorous 0.12389588 0.66073532 3.501481
    ## 26:  0.41627660 Bergen & Wilson   simple 0.32191657 0.93171146 6.884831
    ## 27:  0.49668724 Bergen & Wilson rigorous 0.12704186 0.60936987 3.587207
    ## 28:  0.49668724 Bergen & Wilson   simple 0.27519273 0.90262762 6.994317
    ## 29:  0.59263050 Bergen & Wilson rigorous 0.13046121 0.56058269 3.698843
    ## 30:  0.59263050 Bergen & Wilson   simple 0.16739934 0.79655560 6.841066
    ## 31:  0.70710678 Bergen & Wilson rigorous 0.13324880 0.51052826 3.842140
    ## 32:  0.70710678 Bergen & Wilson   simple 0.26007717 0.85203459 7.247023
    ## 33:  0.84369602 Bergen & Wilson rigorous 0.13533584 0.45757971 4.050357
    ## 34:  0.84369602 Bergen & Wilson   simple 0.23700113 0.80488853 7.290715
    ## 35:  1.00666971 Bergen & Wilson rigorous 0.13840800 0.40724108 4.268582
    ## 36:  1.00666971 Bergen & Wilson   simple 0.27493690 0.79843584 7.285638
    ## 37:  1.20112443 Bergen & Wilson rigorous 0.14003319 0.35401023 4.549172
    ## 38:  1.20112443 Bergen & Wilson   simple 0.25512694 0.74004248 7.184132
    ## 39:  1.43314127 Bergen & Wilson rigorous 0.14167594 0.29631678 4.816802
    ## 40:  1.43314127 Bergen & Wilson   simple 0.35480513 0.77887921 7.065068
    ## 41:  1.70997595 Bergen & Wilson rigorous 0.14176732 0.23421487 5.081253
    ## 42:  1.70997595 Bergen & Wilson   simple 0.24686585 0.60211882 6.474807
    ## 43:  2.04028577 Bergen & Wilson rigorous 0.14017340 0.16755304 5.312550
    ## 44:  2.04028577 Bergen & Wilson   simple 0.35200730 0.65668312 5.813591
    ## 45:  2.43440035 Bergen & Wilson rigorous 0.31911527 1.41993521 5.440962
    ## 46:  2.43440035 Bergen & Wilson   simple 0.29381522 0.54456094 5.452535
    ## 47:  2.90464460 Bergen & Wilson rigorous 0.12688978 0.01201489 5.750513
    ## 48:  2.90464460 Bergen & Wilson   simple 0.32057549 0.60490842 5.731701
    ## 49:  3.46572422 Bergen & Wilson rigorous 0.28590356 1.40161831 6.044159
    ## 50:  3.46572422 Bergen & Wilson   simple 0.27454393 0.58278935 6.067781
    ## 51:  4.13518554 Bergen & Wilson rigorous 0.28009937 1.35236336 6.377757
    ## 52:  4.13518554 Bergen & Wilson   simple 0.25005665 0.60394639 6.356407
    ## 53:  4.93396427 Bergen & Wilson rigorous 0.75617007 1.10694487 6.615564
    ## 54:  4.93396427 Bergen & Wilson   simple 0.19288310 0.58293923 6.617445
    ## 55:  5.88704019 Bergen & Wilson rigorous 0.43538215 1.13823740 7.018521
    ## 56:  5.88704019 Bergen & Wilson   simple 0.18439479 0.68312574 7.271228
    ## 57:  7.02421830 Bergen & Wilson rigorous 0.16820121 1.24434563 7.368181
    ## 58:  7.02421830 Bergen & Wilson   simple 0.04692713 0.11689176 7.294484
    ## 59:  8.38106097 Bergen & Wilson rigorous 0.86744301 1.02898395 7.705786
    ## 60:  8.38106097 Bergen & Wilson   simple 0.09936827 0.75475969 7.985724
    ## 61: 10.00000000 Bergen & Wilson rigorous 0.94189023 1.01411830 8.110762
    ## 62: 10.00000000 Bergen & Wilson   simple 0.03255722 0.58473554 7.909026
    ##              SF           model fit_type         a0         a1       a2
    ##            a3           a4     G        rss
    ##         <num>        <num> <num>      <num>
    ##  1: 12.483355  0.753361116    NA 0.16991528
    ##  2: 11.817349  1.819611262    NA 0.17082065
    ##  3: 12.455652  0.757386850    NA 0.15442511
    ##  4: 11.625446  1.833709132    NA 0.15621110
    ##  5: 12.301100  1.160094513    NA 0.13558048
    ##  6: 11.752737  1.878487125    NA 0.13654196
    ##  7: 12.383837  0.539848603    NA 0.11470873
    ##  8: 11.763206  1.872048568    NA 0.11552021
    ##  9: 12.131575  0.617649150    NA 0.09533333
    ## 10: 11.481767  2.120102903    NA 0.09685217
    ## 11: 12.135860  0.242271681    NA 0.08374252
    ## 12: 11.297596  2.257994412    NA 0.08609498
    ## 13: 11.872356  0.301332381    NA 0.08448159
    ## 14: 11.282791  1.709009751    NA 0.08576512
    ## 15: 11.617210  0.534557077    NA 0.09061620
    ## 16: 11.185678  1.563155360    NA 0.09159669
    ## 17: 11.492982  0.336216338    NA 0.10109214
    ## 18: 10.921553  1.561583631    NA 0.10219068
    ## 19: 11.317396  0.288704946    NA 0.11485638
    ## 20: 10.603117  2.395566149    NA 0.11796866
    ## 21: 20.892272 13.746281468    NA 0.12514061
    ## 22: 10.512398  2.273523001    NA 0.13284259
    ## 23: 20.585224 14.362100537    NA 0.13258102
    ## 24: 10.444733  1.435369475    NA 0.14614198
    ## 25: 20.325222 15.291123555    NA 0.13684856
    ## 26: 10.308401  1.180643316    NA 0.15935598
    ## 27: 19.715624 16.230599431    NA 0.13631106
    ## 28: 10.023470  1.411195489    NA 0.17039760
    ## 29: 19.010284 17.229949605    NA 0.12995315
    ## 30:  9.669929  2.607177475    NA 0.18049191
    ## 31: 18.172361 18.417470464    NA 0.11786647
    ## 32:  9.643738  1.514461987    NA 0.18060771
    ## 33: 17.071657 18.961334353    NA 0.10190885
    ## 34:  9.528558  1.697637220    NA 0.17933572
    ## 35: 16.122299 19.172242723    NA 0.08498549
    ## 36:  9.724603  1.425868085    NA 0.17264735
    ## 37: 15.039235 18.697549435    NA 0.07007367
    ## 38:  9.873408  1.550155404    NA 0.16228418
    ## 39: 14.210561 17.970579265    NA 0.05936271
    ## 40: 10.354918  0.993737087    NA 0.14710549
    ## 41: 13.521056 17.299125898    NA 0.05383649
    ## 42: 11.136842  1.413392752    NA 0.12863055
    ## 43: 13.038262 16.919358092    NA 0.05306939
    ## 44: 12.643951  0.608023142    NA 0.10381028
    ## 45: 13.169277  0.497358402    NA 0.07671302
    ## 46: 12.989870  0.121891743    NA 0.07677388
    ## 47: 12.291709 11.187854017    NA 0.06432940
    ## 48: 12.471458  0.081455753    NA 0.06463389
    ## 49: 12.131949  0.481827417    NA 0.07313517
    ## 50: 11.760647  0.009023445    NA 0.07351417
    ## 51: 11.688665  0.442687670    NA 0.09582860
    ## 52: 11.520097  0.101717764    NA 0.09601457
    ## 53: 11.542983  0.143309377    NA 0.12277040
    ## 54: 11.487541  0.281132726    NA 0.12414392
    ## 55: 11.000106  0.185657341    NA 0.14875787
    ## 56: 10.524370  0.194151759    NA 0.14997967
    ## 57: 10.570984  0.319287259    NA 0.17115746
    ## 58: 10.296601  0.325535745    NA 0.17158073
    ## 59: 10.370947  0.040983920    NA 0.18823481
    ## 60:  9.878284  0.145063544    NA 0.18937137
    ## 61: 10.001531  0.020167184    NA 0.19943157
    ## 62:  9.723625  0.019712916    NA 0.20023839
    ##            a3           a4     G        rss

``` r
Kelly_IRF_coefs[model=="Burr & Morrone"]
```

    ##              SF          model fit_type          a0        a1       a2
    ##           <num>         <char>   <char>       <num>     <num>    <num>
    ##  1:  0.05000000 Burr & Morrone rigorous 0.083084711 15.307892 2.384464
    ##  2:  0.05000000 Burr & Morrone   simple 0.237027640  8.697440 2.652366
    ##  3:  0.05965832 Burr & Morrone rigorous 0.045646397 19.431653 2.206225
    ##  4:  0.05965832 Burr & Morrone   simple 0.254789640  8.726802 2.650678
    ##  5:  0.07118230 Burr & Morrone rigorous 0.049300048 19.279234 2.186592
    ##  6:  0.07118230 Burr & Morrone   simple 0.273125594  8.725043 2.617888
    ##  7:  0.08493232 Burr & Morrone rigorous 0.058185552 18.937479 2.143584
    ##  8:  0.08493232 Burr & Morrone   simple 0.290637085  8.668266 2.564319
    ##  9:  0.10133839 Burr & Morrone rigorous 0.054365000 19.599037 2.227059
    ## 10:  0.10133839 Burr & Morrone   simple 0.312166961  8.599832 2.467226
    ## 11:  0.12091356 Burr & Morrone rigorous 0.053715287 20.056533 2.283389
    ## 12:  0.12091356 Burr & Morrone   simple 0.329800747  8.436747 2.309690
    ## 13:  0.14426999 Burr & Morrone rigorous 0.059112431 19.452248 2.206575
    ## 14:  0.14426999 Burr & Morrone   simple 0.352767394  8.390336 2.269064
    ## 15:  0.17213810 Burr & Morrone rigorous 0.059196897 19.858120 2.256927
    ## 16:  0.17213810 Burr & Morrone   simple 0.374575035  8.242764 2.158668
    ## 17:  0.20538939 Burr & Morrone rigorous 0.063054036 19.602697 2.223983
    ## 18:  0.20538939 Burr & Morrone   simple 0.399485084  8.122389 2.074788
    ## 19:  0.24506371 Burr & Morrone rigorous 0.058110614 20.874941 2.378931
    ## 20:  0.24506371 Burr & Morrone   simple 0.169417519 13.296131 3.321731
    ## 21:  0.29240177 Burr & Morrone rigorous 0.057445614 20.909499 2.381137
    ## 22:  0.29240177 Burr & Morrone   simple 0.173188093 13.599469 3.404421
    ## 23:  0.34888396 Burr & Morrone rigorous 0.063978290 19.979251 2.267589
    ## 24:  0.34888396 Burr & Morrone   simple 0.177007785 13.777435 3.445765
    ## 25:  0.41627660 Burr & Morrone rigorous 0.068000801 20.149240 2.288356
    ## 26:  0.41627660 Burr & Morrone   simple 0.222713963 11.667616 3.084133
    ## 27:  0.49668724 Burr & Morrone rigorous 0.074647554 20.014561 2.271421
    ## 28:  0.49668724 Burr & Morrone   simple 0.230919818 11.540661 3.093363
    ## 29:  0.59263050 Burr & Morrone rigorous 0.048366604 24.170104 2.734956
    ## 30:  0.59263050 Burr & Morrone   simple 0.259115347 10.645117 2.917383
    ## 31:  0.70710678 Burr & Morrone rigorous 0.054922908 22.767741 2.586605
    ## 32:  0.70710678 Burr & Morrone   simple 0.190237891 13.884428 3.218539
    ## 33:  0.84369602 Burr & Morrone rigorous 0.044991945 23.264993 2.637815
    ## 34:  0.84369602 Burr & Morrone   simple 0.236221628 11.408302 3.165174
    ## 35:  1.00666971 Burr & Morrone rigorous 0.041476233 23.811988 2.694236
    ## 36:  1.00666971 Burr & Morrone   simple 0.228804203 11.520901 3.218430
    ## 37:  1.20112443 Burr & Morrone rigorous 0.092775685 17.917675 3.136992
    ## 38:  1.20112443 Burr & Morrone   simple 0.219022236 11.622047 3.265260
    ## 39:  1.43314127 Burr & Morrone rigorous 0.043673426 23.212673 2.630869
    ## 40:  1.43314127 Burr & Morrone   simple 0.205698359 11.824950 3.329052
    ## 41:  1.70997595 Burr & Morrone rigorous 0.042090065 23.896986 3.025660
    ## 42:  1.70997595 Burr & Morrone   simple 0.186454861 12.060123 3.395858
    ## 43:  2.04028577 Burr & Morrone rigorous 0.043541797 22.812734 2.911364
    ## 44:  2.04028577 Burr & Morrone   simple 0.163858896 12.353134 3.468048
    ## 45:  2.43440035 Burr & Morrone rigorous 0.036609450 22.156993 2.840270
    ## 46:  2.43440035 Burr & Morrone   simple 0.171368978 11.352028 3.256118
    ## 47:  2.90464460 Burr & Morrone rigorous 0.031621228 21.344595 2.748821
    ## 48:  2.90464460 Burr & Morrone   simple 0.132404116 12.743465 3.554236
    ## 49:  3.46572422 Burr & Morrone rigorous 0.027947852 21.081345 2.718527
    ## 50:  3.46572422 Burr & Morrone   simple 0.134197137 10.860770 3.156294
    ## 51:  4.13518554 Burr & Morrone rigorous 0.044950853 16.690423 2.970322
    ## 52:  4.13518554 Burr & Morrone   simple 0.105890281 11.069115 3.207372
    ## 53:  4.93396427 Burr & Morrone rigorous 0.034584808 16.742232 2.977992
    ## 54:  4.93396427 Burr & Morrone   simple 0.090874195 10.655448 3.113091
    ## 55:  5.88704019 Burr & Morrone rigorous 0.022367809 16.312097 2.914938
    ## 56:  5.88704019 Burr & Morrone   simple 0.059561278 11.274835 3.256565
    ## 57:  7.02421830 Burr & Morrone rigorous 0.017611606 15.716755 2.822969
    ## 58:  7.02421830 Burr & Morrone   simple 0.041713592 10.079440 2.980627
    ## 59:  8.38106097 Burr & Morrone rigorous 0.009425348 15.670867 2.816283
    ## 60:  8.38106097 Burr & Morrone   simple 0.023927412 10.070033 2.979932
    ## 61: 10.00000000 Burr & Morrone rigorous 0.006690821 14.408183 3.162353
    ## 62: 10.00000000 Burr & Morrone   simple 0.012914814 10.068148 2.979436
    ##              SF          model fit_type          a0        a1       a2
    ##            a3    a4     G       rss
    ##         <num> <num> <num>     <num>
    ##  1:  9.725416   NaN    NA 0.5918945
    ##  2: 19.811016   NaN    NA 0.7830663
    ##  3:  6.789983   NaN    NA 0.5765158
    ##  4: 20.045816   NaN    NA 0.7795097
    ##  5:  6.788581   NaN    NA 0.5578666
    ##  6: 20.262285   NaN    NA 0.7759936
    ##  7:  6.956530   NaN    NA 0.5406605
    ##  8: 20.387114   NaN    NA 0.7724175
    ##  9:  6.740092   NaN    NA 0.5214443
    ## 10: 20.565488   NaN    NA 0.7689082
    ## 11:  6.605944   NaN    NA 0.5052189
    ## 12: 20.476986   NaN    NA 0.7654497
    ## 13:  6.506880   NaN    NA 0.4802044
    ## 14: 20.853224   NaN    NA 0.7620256
    ## 15:  6.372049   NaN    NA 0.4634465
    ## 16: 20.964405   NaN    NA 0.7590258
    ## 17:  6.255231   NaN    NA 0.4500552
    ## 18: 21.256707   NaN    NA 0.7565703
    ## 19:  6.033481   NaN    NA 0.4401330
    ## 20: 12.009567   NaN    NA 0.5978596
    ## 21:  5.729690   NaN    NA 0.4265873
    ## 22: 11.807474   NaN    NA 0.5865403
    ## 23:  5.602059   NaN    NA 0.4284388
    ## 24: 11.587991   NaN    NA 0.5766654
    ## 25:  5.560992   NaN    NA 0.4313405
    ## 26: 11.709021   NaN    NA 0.6116801
    ## 27:  5.543954   NaN    NA 0.4512461
    ## 28: 11.512826   NaN    NA 0.6101449
    ## 29:  4.522129   NaN    NA 0.4477757
    ## 30: 11.568665   NaN    NA 0.6294437
    ## 31:  4.442191   NaN    NA 0.4192351
    ## 32: 10.604232   NaN    NA 0.5698066
    ## 33:  3.599176   NaN    NA 0.4080456
    ## 34: 10.539452   NaN    NA 0.6048691
    ## 35:  3.111958   NaN    NA 0.4132730
    ## 36: 10.015824   NaN    NA 0.6045999
    ## 37:  5.376894   NaN    NA 0.4796745
    ## 38:  9.419019   NaN    NA 0.6080660
    ## 39:  2.891601   NaN    NA 0.4589990
    ## 40:  8.717136   NaN    NA 0.6135596
    ## 41:  2.212191   NaN    NA 0.4581680
    ## 42:  7.832594   NaN    NA 0.6219274
    ## 43:  2.395020   NaN    NA 0.4795370
    ## 44:  6.839410   NaN    NA 0.6334491
    ## 45:  1.786335   NaN    NA 0.4794058
    ## 46:  7.358158   NaN    NA 0.6646793
    ## 47:  1.440453   NaN    NA 0.5141681
    ## 48:  5.614881   NaN    NA 0.6782668
    ## 49:  1.315550   NaN    NA 0.5491851
    ## 50:  6.588323   NaN    NA 0.7071591
    ## 51:  2.814833   NaN    NA 0.5705736
    ## 52:  5.928316   NaN    NA 0.7185300
    ## 53:  2.515210   NaN    NA 0.5887596
    ## 54:  6.087239   NaN    NA 0.7428257
    ## 55:  1.846599   NaN    NA 0.6000988
    ## 56:  5.120855   NaN    NA 0.7452479
    ## 57:  2.492885   NaN    NA 0.6505283
    ## 58:  5.408744   NaN    NA 0.7693013
    ## 59:  1.898843   NaN    NA 0.6589544
    ## 60:  4.978939   NaN    NA 0.7754218
    ## 61:  2.396094   NaN    NA 0.7170691
    ## 62:  4.974296   NaN    NA 0.7873461
    ##            a3    a4     G       rss

``` r
# model comparison?
t.test(y = Kelly_IRF_coefs[model=="Bergen & Wilson" & fit_type=="rigorous", rss], 
       x = Kelly_IRF_coefs[model=="Burr & Morrone" & fit_type=="rigorous", rss], paired = TRUE)
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  Kelly_IRF_coefs[model == "Burr & Morrone" & fit_type == "rigorous", rss] and Kelly_IRF_coefs[model == "Bergen & Wilson" & fit_type == "rigorous", rss]
    ## t = 33.76, df = 30, p-value < 2.2e-16
    ## alternative hypothesis: true mean difference is not equal to 0
    ## 95 percent confidence interval:
    ##  0.3719070 0.4197998
    ## sample estimates:
    ## mean difference 
    ##       0.3958534

# 4. Create a continuous TRF space and save it for later

Have a look at the individual IRFs for each SF, then fit a Time x SF GAM
so that we can interpolate an IRF for every possible SF. Then we’ll save
the model for later…

``` r
## IRFs
p_irf_burrmorrone <- ggplot(data = Kelly_IRFs, aes(x = irf_time, y = irf_burr, color = SF, group = SF )) + 
  geom_line(size = 1.5) + theme_classic() + 
  scale_color_viridis_c(trans = "log10") + 
  labs(color = "SF [cpd]", x = "Time [ms]", y = "H(t)") + 
  facet_wrap(~fit_type)
p_irf_burrmorrone
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# Bergen & Wilson
p_irf_bergenwilson <- ggplot(data = Kelly_IRFs, aes(x = irf_time, y = irf_bergen, color = SF, group = SF )) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_classic(base_size = 12.5) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_color_viridis_c(trans = "log10") + 
  labs(color = "SF [cpd]", x = "Time [ms]", y = "H(t)") + 
  facet_wrap(~fit_type)
p_irf_bergenwilson
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
# IRF surface
p_irf_surface <- ggplot(data = Kelly_IRFs, 
                        aes(x = irf_time, y = SF, fill = irf_bergen)) + 
  geom_raster(interpolate = TRUE) + 
  coord_cartesian(expand = FALSE) + 
  scale_x_continuous(limits = c(0, 250), expand = c(0,0)) + 
  scale_y_log10() + 
  annotation_logticks(sides = "l") + 
  scale_fill_gradient2() + 
  theme_minimal(base_size = 12.5) + theme() + 
  labs(y = "Spatial frequency [cpd]", x = "Time [ms]", fill = "H(t)") + 
  facet_wrap(~fit_type)
p_irf_surface
```

    ## Warning: Removed 4588 rows containing missing values or values outside the scale range
    ## (`geom_raster()`).

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
# on the side, a few exemplary IRFs
p_irf_examples <- ggplot(data = Kelly_IRFs[fit_type=="simple"], 
                        aes(x = irf_time, y = irf_bergen)) + 
  scale_x_continuous(limits = c(0, 250), expand = c(0,0)) + 
  geom_line(size = 1.5) + 
  theme_minimal(base_size = 12.5) + theme() + 
  labs(x = "Time [ms]", y = "H(t)") + 
  facet_wrap(~round(SF,2))
p_irf_examples
```

    ## Warning: Removed 73 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

``` r
# approximate the IRFs with a GAM ~ time x SF
gam_irf <- bam(data = Kelly_IRFs, subset = fit_type=="simple", 
               formula = irf_bergen ~ 
                 s(log(SF), k = 10, bs = "cr") + # initially 16
                 s(irf_time, k = 50, bs = "cr") + # initially 40
                 ti(log(SF), irf_time, k = c(10, 50), bs = c("cr", "cr") )
)
summary(gam_irf)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## irf_bergen ~ s(log(SF), k = 10, bs = "cr") + s(irf_time, k = 50, 
    ##     bs = "cr") + ti(log(SF), irf_time, k = c(10, 50), bs = c("cr", 
    ##     "cr"))
    ## 
    ## Parametric coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 9.293e-04  7.642e-07    1216   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                          edf Ref.df      F p-value    
    ## s(log(SF))             8.997    9.0 138117  <2e-16 ***
    ## s(irf_time)           48.816   49.0 268159  <2e-16 ***
    ## ti(log(SF),irf_time) 406.122  433.4  13936  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.999   Deviance explained = 99.9%
    ## fREML = -1.0503e+05  Scale est. = 7.8253e-09  n = 13423

``` r
Kelly_IRFs[fit_type=="simple", irf_bergen_gam_pred := predict(gam_irf)]

# IRF surface, fitted by GAM
p_irf_surface_gam <- ggplot(data = Kelly_IRFs[fit_type=="simple"], 
                        aes(x = irf_time, y = SF, fill = irf_bergen_gam_pred)) + 
  geom_raster(interpolate = TRUE) + 
  coord_cartesian(expand = FALSE) + 
  scale_x_continuous(limits = c(0, 250), expand = c(0,0)) + 
  scale_y_log10() +  # limits = c(0.01, 10)
  annotation_logticks(sides = "l") + 
  scale_fill_gradient2() + 
  theme_minimal(base_size = 12.5) + theme() + 
  labs(y = "Spatial frequency [cpd]", x = "Time [ms]", fill = "H(t)")
p_irf_surface_gam
```

    ## Warning: Removed 2294 rows containing missing values or values outside the scale range
    ## (`geom_raster()`).

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->

``` r
# save the models
if (do_save) {
  save(gam_irf, file = "gam_irf.rda", compress = "xz")
}
```

# 5. Fourier analysis of TRFs

Look at the Fourier domain of those IRFs. The amplitude spectrum should
look very much like Kelly’s sensitivity curve!

``` r
# pick a single SF
unique_SFs <- unique(Kelly_IRFs$SF)
IRF_ffts <- NULL
for (unique_SF in unique_SFs) {
  # run fft
  single_IRF <- Kelly_IRFs[fit_type=="simple" & SF==unique_SF]
  # zero padding?
  n = 2^nextpow2(nrow(single_IRF)*10)+1;
  zeros = rep(0, times = n - nrow(single_IRF))
  # from https://de.mathworks.com/help/matlab/ref/fft.html
  Y <- fft(c(single_IRF$irf_bergen_gam_pred, zeros))
  n <- length(Y) # length of signal 
  Fs <- unique(round(1/diff(single_IRF$irf_time/1000))) # in Hz, artifical sampling rate
  # amplitude
  P2 = abs(Y) #abs(Y/n) 
  P1 = P2[1:floor(n/2+1)]
  P1[2:(length(P1)-1)] = 2*P1[2:(length(P1)-1)]
  # phase
  pha2 = Imag(Y) 
  pha1 = pha2[1:floor(n/2+1)]
  pha1[2:(length(pha1)-1)] = 2*pha1[2:(length(pha1)-1)]
  # save:
  single_IRF_fft <- data.table(SF = unique_SF, 
                               f = Fs*(0:(n/2))/n, 
                               Y = P1, 
                               pha = pha1)
  IRF_ffts <- rbind(IRF_ffts, single_IRF_fft)
}
# plot the fourier domain
p_fft_amp <- ggplot(IRF_ffts, aes(x = f, y = Y, color = SF, group = SF)) + 
  #coord_cartesian(expand = FALSE, xlim = c(0.5, 50), ylim = c(1e-7, 1e-3/2)) + 
  coord_cartesian(expand = FALSE, xlim = c(0.5, 50), ylim = c(0.01, 200)) + 
  geom_line(size = 1.5, alpha = 1) + 
  scale_x_log10() + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) + 
  annotation_logticks() + 
  theme_classic(base_size = 12.5) + theme(legend.position = "bottom") + 
  scale_color_viridis_c(trans = "log10") + 
  labs(x = "Temporal frequency [Hz]", y = "Amplitude", color = "Spatial frequency [cpd]", 
       title = "IRF model Fourier domain - amplitude")
# plot the phase
p_fft_phase <- ggplot(IRF_ffts, aes(x = f, y = pha, color = SF, group = SF)) + 
  coord_cartesian(expand = FALSE, xlim = c(0.5, 50)) + 
  geom_line(size = 1.5, alpha = 1) + 
  scale_x_log10() + 
  annotation_logticks(sides = "b") + 
  theme_classic(base_size = 12.5) + theme(legend.position = "bottom") + 
  scale_color_viridis_c(trans = "log10") + 
  labs(x = "Temporal frequency [Hz]", y = "Phase", color = "Spatial frequency [cpd]", 
       title = "IRF model Fourier domain - phase")

# combine
p_fft_combined <- plot_grid(p_groundtruth, 
                            p_fft_amp, p_fft_phase, nrow = 1)
```

    ## Warning in scale_x_log10(): log-10 transformation introduced infinite values.
    ## log-10 transformation introduced infinite values.

``` r
p_fft_combined
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

The Fourier spectrum seems to match - how about some more explicit model
predictions?

``` r
# predict Kelly's SF-TF surface using the model
gam_and_flicker_to_sens <- function(tfreqs,  # can be numeric vector
                                    SF_now, # single value, for now
                                    irf_gam, # is the SF x time GAM that encodes the IRFs
                                    gam_df=NULL, # if you want to not specify a GAM, but a DF containing the IRFs
                                    irf_dur=500, # 
                                    t_res, 
                                    pres_dur = 4 * 8 * 28 * (1000/120), 
                                    ramp_sd = 4 * 28 * (1000/120), # those are constants
                                    beta=1.5, G=1 # slope of prob summation, negative-lobe asymmetry
) {
  require(data.table)
  require(mgcv)
  require(pracma)
  # create temporal reference
  stim_time <- seq(0, pres_dur, by = t_res) 
  # iterate over tfreqs
  sums_tfreqs <- vector(mode = "numeric", length = length(tfreqs))
  for (tfreq_i in (1:length(tfreqs))) {
    # create flicker stimulus and gaussian temporal ramp
    x <- create_flicker_stim(dur = pres_dur, t_res = t_res, 
                             tfreq = tfreqs[tfreq_i]) 
    if (!is.na(ramp_sd)) {
      gaussian <- dnorm(x = stim_time, mean = pres_dur/2, sd = ramp_sd) # gaussian ramp
      gaussian <- gaussian / max(gaussian)
    } else {
      gaussian <- rep(1, length.out = length(stim_time))
    }
    stim_tfreq <- x * gaussian # combined
    # extract IRF from GAM
    if (is.null(gam_df)) {
      gam_df <- data.table(irf_time = seq(0, irf_dur, by = t_res), 
                         SF = SF_now)
      gam_df[ , irf_predict := predict.bam(object = irf_gam, newdata = gam_df)]
      trf_time <- gam_df$irf_time
      trf <- gam_df$irf_predict
    } else {
      trf_time <- as.numeric(t(gam_df[,1]))
      trf <- as.numeric(t(gam_df[,2]))
    }
    # convolution
    resp_tfreq <- zapsmall(convolve(stim_tfreq, rev(trf), type = "open"))[1:length(stim_tfreq)]
    # add the systematic asymmetry? (as in Kelly & Savoie, 1978 or Bergen & Wilson, 1984)
    if (G != 1) {
      resp_tfreq[resp_tfreq<0] <- resp_tfreq[resp_tfreq<0] * G
    }
    # probability summation
    sums_tfreqs[tfreq_i] <- trapz(x = stim_time, y = abs(resp_tfreq)^beta)^(1/beta)
  }
  return(sums_tfreqs)
}

# in case we do not want to read out directly from the GAM
# e.g., to see the effect of fitting by the GAM
not_gam_alternative <- Kelly_IRFs[fit_type=="simple", c("irf_time", "irf_bergen", "SF")]
are_equal(sort(unique(Kelly_stabilized$SF)), sort(unique(not_gam_alternative$SF)))
```

    ## [1] TRUE

``` r
# let the TRFs predict
Kelly_stabilized[ , fit_ultimate := NaN]
for (r in 1:nrow(Kelly_stabilized)) {
  Kelly_stabilized$fit_ultimate[r] <- gam_and_flicker_to_sens(tfreqs = Kelly_stabilized$TF[r], 
                                                              SF = Kelly_stabilized$SF[r], 
                                                              irf_gam = gam_irf, 
                                                              gam_df = NULL,
                                                              #gam_df = not_gam_alternative[SF==Kelly_stabilized$SF[r]], 
                                                              t_res = t_res)
}



# have a look:
p_kelly_surface_predicted <- ggplot(data = Kelly_stabilized, aes(x = SF, y = TF, z = log10(fit_ultimate) )) + 
  #geom_raster(interpolate = TRUE) + 
  geom_contour() + geom_contour_filled(bins = 30) + 
  coord_cartesian(expand = FALSE) + 
  scale_x_log10() + scale_y_log10() + 
  annotation_logticks() + 
  #scale_fill_viridis_c(trans = "log10") +
  theme_classic(base_size = 12.5) + theme() + 
  labs(x = "Spatial frequency [cpd]", y = "Temporal frequency [Hz]", fill = "Contrast\nsensitivity")

p_kelly_surface_original_and_predicted <- 
  plot_grid(p_kelly_surface + scale_fill_viridis_d(option = "mako") + 
              theme(legend.position = "none") + 
              ggtitle("Kelly's original function"), 
            p_kelly_surface_predicted + scale_fill_viridis_d(option = "mako") + 
              theme(legend.position = "none") + 
              ggtitle("IRF model prediction"), 
            nrow = 1)
p_kelly_surface_original_and_predicted
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

… this is not exactly like Kelly’s surface, but good enough. The TRFs
are able to explain a very good deal of the data!

Finally, save this workspace:

``` r
if (do_save) {
  save.image(file = "estimate_TRFs_for_SFs.RData", compress = "xz")
}
```

Make a plot for the paper:

``` r
SDECTheme <- function(base_size=15, base_family="Helvetica") { 
  theme_classic(base_size=base_size, base_family=base_family) %+replace% 
    theme(
      # size of text
      axis.text = element_text(size = base_size), # 0.9*base_size
      legend.text = element_text(size = base_size),
      # remove grid horizontal and vertical lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background  = element_blank(),
      # panel boxes
      panel.border = element_blank(),
      axis.line = element_line(colour = "grey20"), 
      # facet strips
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill="transparent", colour = "transparent"),
      strip.placement = "outside",
      # legend
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}

#### VERSION 1: surface

# for labeled contours
require(metR) # needs: sudo apt install libudunits2-dev
```

    ## Loading required package: metR

    ## 
    ## Attaching package: 'metR'

    ## The following object is masked from 'package:pracma':
    ## 
    ##     cross

``` r
# define breaks manually
breaks_here <- (c(0, exp(seq(log(0.12), log(120), length.out = 100))))
breaks_here_2 <- c(0, 0.1, 1, 7.5, 25, 50, 75, 100, 120)
# original kelly
p_kelly_surface_paper <- ggplot(data = Kelly_stabilized, aes(x = SF, y = TF, z = sens )) + 
  geom_contour_fill(bins = 100, breaks = breaks_here) + 
  geom_contour(colour = "black", alpha = 0.5, breaks = breaks_here_2) + 
  geom_text_contour(breaks = breaks_here_2,
                      label.placer = label_placer_fraction(frac = 0.9)) + 
  #geom_contour_filled(bins = 30) + geom_contour(bins = 30) + 
  coord_cartesian(expand = FALSE) + 
  scale_x_log10() + scale_y_log10() + 
  annotation_logticks() + 
  scale_fill_viridis_c(trans = "log10", option = "mako", breaks = c(0.1, 1, 10, 100)#, 
                       # limits = c(0.1, 
                       #            max(c(Kelly_stabilized$sens, Kelly_stabilized$fit_ultimate)))
                       ) +
  theme_classic(base_size = 12.5) + SDECTheme() + 
  labs(x = "Spatial frequency [cpd]", y = "Temporal frequency [Hz]", fill = "Contrast\nsensitivity")

p_kelly_surface_predicted_paper <- ggplot(data = Kelly_stabilized, aes(x = SF, y = TF, z = fit_ultimate )) + 
  geom_contour_fill(bins = 100, breaks = breaks_here) + 
  geom_contour(colour = "black", alpha = 0.5, breaks = breaks_here_2) + 
  geom_text_contour(breaks = breaks_here_2,
                      label.placer = label_placer_fraction(frac = 0.9)) + 
  #geom_contour_filled(bins = 30) + geom_contour(bins = 30) + 
  coord_cartesian(expand = FALSE) + 
  scale_x_log10() + scale_y_log10() + 
  annotation_logticks() + 
  scale_fill_viridis_c(trans = "log10", option = "mako", breaks = c(0.1, 1, 10, 100)#, 
                       # limits = c(0.1, 
                       #            max(c(Kelly_stabilized$sens, Kelly_stabilized$fit_ultimate)))
                       ) +
  theme_classic(base_size = 12.5) + SDECTheme() + 
  labs(x = "Spatial frequency [cpd]", y = "Temporal frequency [Hz]", fill = "Contrast\nsensitivity")

plot_grid(p_kelly_surface_paper, p_kelly_surface_predicted_paper)
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
##### VERSION 2: contrast sensitivity plot
p_kelly_surface_paper_2 <- ggplot(data = Kelly_stabilized, aes(x = TF, y = sens, 
                                                               color = SF, group = SF)) + 
  geom_line(size = 1.5) + 
  #geom_contour_filled(bins = 30) + geom_contour(bins = 30) + 
  coord_cartesian(expand = FALSE, ylim = c(0.1, 200)) + 
  scale_y_log10(breaks = c(0.1, 1, 10, 100)) + 
  scale_x_log10() + 
  annotation_logticks() + 
  scale_color_viridis_c(trans = "log10") +
  theme_classic(base_size = 12.5) + SDECTheme() + 
  labs(color = "SF [cpd]", x = "Temporal frequency [Hz]", 
       y = "Contrast sensitivity (Kelly, 1979)")
p_kelly_surface_paper_2
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
p_kelly_surface_predicted_paper_2 <- ggplot(data = Kelly_stabilized, aes(x = TF, y = fit_ultimate, 
                                                               color = SF, group = SF)) + 
  geom_line(size = 1.5) + 
  #geom_contour_filled(bins = 30) + geom_contour(bins = 30) + 
  coord_cartesian(expand = FALSE, ylim = c(0.1, 200)) + 
  scale_y_log10(breaks = c(0.1, 1, 10, 100)) + 
  scale_x_log10() + 
  annotation_logticks() + 
  scale_color_viridis_c(trans = "log10") +
  theme_classic(base_size = 12.5) + SDECTheme() + 
  labs(color = "SF [cpd]", x = "Temporal frequency [Hz]", 
       y = "Predicted contrast sensitivity")
p_kelly_surface_predicted_paper_2
```

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

``` r
# show IRFs
p_irf_bergenwilson_paper <- ggplot(data = Kelly_IRFs[fit_type=="simple"], 
                                   aes(x = irf_time, y = irf_bergen, 
                                       color = SF, group = SF )) + 
  geom_line(size = 1.5, alpha = 0.8) + 
  theme_classic(base_size = 12.5) + SDECTheme() + 
  scale_x_continuous(expand = c(0,0), limits = c(0, 250)) + 
  scale_color_viridis_c(trans = "log10") + 
  labs(color = "SF [cpd]", x = "Time [ms]", y = "H(t)") 
p_irf_bergenwilson_paper
```

    ## Warning: Removed 2263 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

``` r
p_irf_surface_gam_paper <- p_irf_surface_gam + SDECTheme()

# combine all
plot_grid(p_kelly_surface_paper_2, p_irf_bergenwilson_paper, 
          p_irf_surface_gam_paper, p_kelly_surface_predicted_paper_2, 
          nrow = 2, ncol = 2, align = "hv") # export as 5 x 9
```

    ## Warning: Removed 2263 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

    ## Warning: Removed 2294 rows containing missing values or values outside the scale range
    ## (`geom_raster()`).

![](estimate_TRFs_for_SFs_files/figure-gfm/unnamed-chunk-15-5.png)<!-- -->
