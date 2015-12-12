######################################
#Convert Sculptor data into form for abc
######################################

###################
#Functions to convert
#equatorial to 'standard'
#coordinates
###################
conv.xi <- function(alpha,delta,alpha0,delta0){
  num <- cos(delta)*sin(alpha-alpha0)
  denom <- sin(delta0)*sin(delta) + cos(delta0)*cos(delta)*cos(alpha-alpha0)
  return(num/denom)
}

conv.eta <- function(alpha,delta,alpha0,delta0){
  num <- cos(delta0)*sin(delta)-sin(delta0)*cos(delta)*cos(alpha-alpha0)
  denom <- sin(delta0)*sin(delta) + cos(delta0)*cos(delta)*cos(alpha-alpha0)
  return(num/denom)
}

###################
#Convert Sculptor data
###################

sculptor.ra.rad <- (1 + 0 / 60 + 9.4 / 60^2)*2*pi / 24
sculptor.dec.rad <- (-33 - 42 / 60 - 33 / 60^2) * 2 * pi / 360
sculptor.dist.kpc <- 86
sculptor.velocity.kms <- 110

stars <- read.table("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/data/b08mg_members.res")
names(stars) <- c("r.ascension.deg", "declination.deg", "velocity.kms",
                  "ve.error", "metallicity", "me.error", "mg.index", 
                  "mg.error", "p.membership","class.m")

stars$r.ascension.rad <- stars$r.ascension.deg * 2 * pi / 360
stars$declination.rad <- stars$declination.deg * 2 * pi / 360

stars$xi <- conv.xi(stars$r.ascension.rad, stars$declination.rad, sculptor.ra.rad, sculptor.dec.rad)
stars$eta <- conv.eta(stars$r.ascension.rad, stars$declination.rad, sculptor.ra.rad, sculptor.dec.rad)

stars.sculptor <- stars[stars$p.membership > .5,]
stars.sculptor$x <- tan(stars.sculptor$xi)*sculptor.dist.kpc
stars.sculptor$y <- tan(stars.sculptor$eta)*sculptor.dist.kpc
stars.sculptor$vz <- stars.sculptor$velocity.kms - sculptor.velocity.kms

stars.sculptor.mr <- stars.sculptor[stars.sculptor$class.m == 1,c("x","y","vz")]
stars.sculptor.mp <- stars.sculptor[stars.sculptor$class.m == 2,c("x","y","vz")]

save(stars.sculptor.mr, file = "/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/data/sculptor_mr_abc.Rdata")
save(stars.sculptor.mp, file = "/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/data/sculptor_mp_abc.Rdata")