source("fernSpike/helpers.R")

## 3 kyr baseline
b3k = sim(list(inj = 1,
               fb_bio = 0,
               fb_oc = 0,
               fb_fpoc = 0,
               fb_ow = 0,
               duration = 3000.0,
               injmass = 0.242,
               assfb = 0.5,
               rlim = 1.000003,
               casename = "base3k"))

## 3 kyr feedback
f3k = sim(list(inj = 1,
               fb_bio = 1,
               fb_oc = 1,
               fb_fpoc = 0,
               fb_ow = 2,
               duration = 3000.0,
               injmass = 0.242,
               assfb = 0.5,
               rlim = 1.000003,
               casename = "fb3k"))

plot.case(f3k, b3k, "fb3k")

write.csv(b3k, "fernSpike/control.csv")
write.csv(f3k, "fernSpike/feedbacks.csv")

## Sensitivity tests - assfb
f3k.weak = sim(list(inj = 1,
               fb_bio = 1,
               fb_oc = 1,
               fb_fpoc = 0,
               fb_ow = 2,
               duration = 3000.0,
               injmass = 0.242,
               assfb = 0.3,
               rlim = 1.000003,
               casename = "weak3k"))

f3k.strong = sim(list(inj = 1,
               fb_bio = 1,
               fb_oc = 1,
               fb_fpoc = 0,
               fb_ow = 2,
               duration = 3000.0,
               injmass = 0.242,
               assfb = 0.7,
               rlim = 1.000003,
               casename = "strong3k"))

write.csv(f3k.weak, "fernSpike/weak.csv")
write.csv(f3k.strong, "fernSpike/strong.csv")

plot.case(f3k.strong, f3k, "strong")
plot.case(f3k.weak, f3k, "weak")

## Sensitivity tests - recovery
f3k.fast = sim(list(inj = 1,
                    fb_bio = 1,
                    fb_oc = 1,
                    fb_fpoc = 0,
                    fb_ow = 2,
                    duration = 3000.0,
                    injmass = 0.242,
                    assfb = 0.3,
                    rlim = 1.00001,
                    casename = "fast3k"))

f3k.slow = sim(list(inj = 1,
                      fb_bio = 1,
                      fb_oc = 1,
                      fb_fpoc = 0,
                      fb_ow = 2,
                      duration = 3000.0,
                      injmass = 0.242,
                      assfb = 0.7,
                      rlim = 1.000001,
                      casename = "slow3k"))

write.csv(f3k.fast, "fernSpike/fast.csv")
write.csv(f3k.slow, "fernSpike/slow.csv")

plot.case(f3k.fast, f3k, "fast")
plot.case(f3k.slow, f3k, "slow")
