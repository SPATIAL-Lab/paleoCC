source("fernSpike/helpers.R")

# 1 kyr injection
## Baseline for 1 kyr injection
base1k = sim(list(inj = 1,
                  fb_bio = 0,
                  fb_oc = 0,
                  fb_fpoc = 0,
                  fb_ow = 0,
                  duration = 1000.0,
                  injmass = 0.3,
                  assfb = 2.0e4,
                  casename = "base1k"))

## Feedbacks for 1 kyr injection
fb1k = sim(list(inj = 1,
                  fb_bio = 1,
                  fb_oc = 1,
                  fb_fpoc = 0,
                  fb_ow = 1,
                  duration = 1000.0,
                  injmass = 0.3,
                  assfb = 2.0e4,
                  casename = "fb1k"))

plot.case(fb1k, base1k, "fb1k")

## Scaled feedbacks for 1 kyr injection
sfb1k = sim(list(inj = 1,
                fb_bio = 1,
                fb_oc = 1,
                fb_fpoc = 0,
                fb_ow = 1,
                duration = 1000.0,
                injmass = 0.3,
                assfb = 2.0e3,
                casename = "sfb1k"))

plot.case(sfb1k, base1k, "sfb1k")

# 
## Baseline for 5 kyr injection
base5k = sim(list(inj = 1,
                  fb_bio = 0,
                  fb_oc = 0,
                  fb_fpoc = 0,
                  fb_ow = 0,
                  duration = 5000.0,
                  injmass = 0.3,
                  assfb = 2.0e4,
                  casename = "base1k"))

## Feedbacks for 5 kyr injection
fb5k = sim(list(inj = 1,
                fb_bio = 1,
                fb_oc = 1,
                fb_fpoc = 0,
                fb_ow = 1,
                duration = 5000.0,
                injmass = 0.3,
                assfb = 2.0e4,
                casename = "fb1k"))

plot.case(fb5k, base5k, "fb5k")

## Scaled feedbacks for 5 kyr injection
sfb5k = sim(list(inj = 1,
                 fb_bio = 1,
                 fb_oc = 1,
                 fb_fpoc = 0,
                 fb_ow = 1,
                 duration = 5000.0,
                 injmass = 0.3,
                 assfb = 1.0e4,
                 casename = "sfb5k"))

plot.case(sfb5k, base5k, "sfb5k")
