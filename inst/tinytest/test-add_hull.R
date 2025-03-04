# Test added on 20250304
info <- "add_hull() returns expected values"
dat <- iris[, 1:2]
# Start NULL graphics device
pdf(NULL)
plot(dat)
ihull <- add_hull(dat)
dev.off()
expect_equal(range(ihull), c(14, 132), info = info)
